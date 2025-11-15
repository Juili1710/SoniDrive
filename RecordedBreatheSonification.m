% RECORDED FILE SONIFICATION APP
% Reads a recorded audio file (.wav, .mp3) block-by-block to simulate live input,
% bypassing the need for a live microphone.
% Displays FFT, Envelope, Spectrogram, and Karaoke plots with 3 Mapping choices.

% --- GLOBAL VARIABLES AND INITIAL CONFIGURATION ---
global timerHandle;
global h_fig h_ax_karaoke h_bpm_text h_status_text;
global h_map_type_dd h_waveform_dd;
global h_ax_spectrogram h_ax_envelope h_ax_fft h_ax_raw;

% Core data storage
global audio_file_data Fs_audio current_block_index; 

% CORE CONFIGURATION PARAMETERS
global Fs Update_Interval_Sec Target_Pitch_Hz Target_BPM Target_Range;
global prev_breathing_data buffer_size lpf_breathing;

Fs = 44100;                      % Target internal processing sample rate
Update_Interval_Sec = 0.5;       % How often to process a block of data (500 ms)
Target_Pitch_Hz = 440;           % A4 note for target BPM (14 BPM)
Target_BPM = 14;                 % Center of the Normal range
Target_Range = [12 15];          % Normal (Green) BPM range
buffer_size = Fs * 10;           % 10 seconds of history for stable BPM calculation

% DSP Filter Setup (Bandpass 0.1 to 1 Hz for Breathing Envelope)
low_freq = 0.1 / (Fs/2);
high_freq = 1 / (Fs/2);
[b_breathing, a_breathing] = butter(4, [low_freq high_freq], 'bandpass');
lpf_breathing = struct('b', b_breathing, 'a', a_breathing);


% --- 1. SETUP AND GUI CREATION ---
function recorded_file_sonification
    global h_fig h_ax_karaoke h_bpm_text h_status_text;
    global h_map_type_dd h_waveform_dd h_ax_spectrogram h_ax_envelope h_ax_fft h_ax_raw Target_Range;

    % --- Create GUI Figure and Layout Adjustments ---
    h_fig = figure('Name', 'Recorded File Sonification (File Input)', 'Position', [100 100 1200 800], 'CloseRequestFcn', @cleanUp);
    
    % --- Control Panel (UI elements) ---
    control_x = 20;
    control_width = 150;
    
    uicontrol('Style', 'text', 'String', 'Mapping Type:', 'Position', [control_x 720 control_width 20], 'HorizontalAlignment', 'left');
    h_map_type_dd = uicontrol('Style', 'popupmenu', 'String', {'Linear', 'Non-Linear', 'Fibonacci'}, ...
                              'Position', [control_x 690 control_width 25], 'Value', 1);

    uicontrol('Style', 'text', 'String', 'Waveform:', 'Position', [control_x 660 control_width 20], 'HorizontalAlignment', 'left');
    h_waveform_dd = uicontrol('Style', 'popupmenu', 'String', {'Sine', 'Chirp'}, ...
                              'Position', [control_x 630 control_width 25], 'Value', 1);
    
    % BPM Value Display
    h_bpm_text = uicontrol('Style', 'text', 'String', 'BPM: --', ...
                           'Position', [control_x 570 control_width 30], 'FontSize', 16, 'FontWeight', 'bold', 'ForegroundColor', [0 0 0.8], 'HorizontalAlignment', 'left');
    
    % Start/Stop Buttons
    uicontrol('Style', 'pushbutton', 'String', 'LOAD FILE & START', 'Position', [control_x 480 control_width 50], 'Callback', @start_file_processing, 'FontSize', 14, 'BackgroundColor', [0.7 1 0.7]);
    uicontrol('Style', 'pushbutton', 'String', 'STOP', 'Position', [control_x 420 control_width 50], 'Callback', @stop_processing, 'FontSize', 14, 'BackgroundColor', [1 0.7 0.7]);
    
    h_status_text = uicontrol('Style', 'text', 'String', 'Status: Stopped', 'Position', [control_x 390 control_width 20], 'FontSize', 10, 'ForegroundColor', [0.8 0 0], 'HorizontalAlignment', 'center');

    % --- Plotting Axes (Manual Positions for Clean Layout) ---
    plot_start_x = 0.2; 
    plot_width = 0.75;
    plot_width_tertiary = plot_width / 3;
    plot_height_top = 0.25;
    
    % TOP ROW (DSP CHAIN)
    h_ax_raw = axes('Parent', h_fig, 'Position', [plot_start_x 0.7 0.24 plot_height_top]); 
    title(h_ax_raw, '1. Raw Audio Block (0.5s)');
    xlabel(h_ax_raw, 'Time (s)');
    
    h_ax_fft = axes('Parent', h_fig, 'Position', [plot_start_x + plot_width_tertiary 0.7 0.24 plot_height_top]); 
    title(h_ax_fft, '2. Filtered FFT (Breathing Freq)');
    xlabel(h_ax_fft, 'Frequency (Hz)');
    
    h_ax_envelope = axes('Parent', h_fig, 'Position', [plot_start_x + 2*plot_width_tertiary 0.7 0.24 plot_height_top]); 
    title(h_ax_envelope, '3. Filtered Breathing Envelope (10s)');
    xlabel(h_ax_envelope, 'Time (s)');

    % MIDDLE/BOTTOM ROW (FEEDBACK)
    h_ax_spectrogram = axes('Parent', h_fig, 'Position', [plot_start_x 0.37 plot_width 0.3]); 
    title(h_ax_spectrogram, '4. Sonification Spectrogram (Output Pitch)');
    
    h_ax_karaoke = axes('Parent', h_fig, 'Position', [plot_start_x 0.05 plot_width 0.3]); 
    title(h_ax_karaoke, ['5. Karaoke Feedback | Target: ' num2str(Target_Range(1)) '-' num2str(Target_Range(2)) ' BPM']);
    ylabel(h_ax_karaoke, 'BPM Level');
    ylim(h_ax_karaoke, [5 25]);
    xlabel(h_ax_karaoke, 'Time (s)'); 
    grid(h_ax_karaoke, 'on');

end

% --- 2. START/STOP FUNCTIONS ---

function start_file_processing(~, ~)
    global timerHandle h_status_text;
    
    if ~isempty(timerHandle) && isvalid(timerHandle) && strcmp(timerHandle.Running, 'on')
        set(h_status_text, 'String', 'Status: Already Running', 'ForegroundColor', [0.8 0.8 0]);
        return;
    end
    
    % --- STEP 1: FILE SELECTION AND LOAD ---
    [filename, pathname] = uigetfile({'*.wav'; '*.mp3'}, 'Select Recorded Breathing Audio File');
    
    if isequal(filename, 0)
        set(h_status_text, 'String', 'Status: File selection cancelled', 'ForegroundColor', [0.8 0 0]);
        return;
    end
    
    filepath = fullfile(pathname, filename);
    set(h_status_text, 'String', 'Status: Loading audio file...', 'ForegroundColor', [0 0 0.8]);
    drawnow;

    try
        global audio_file_data Fs_audio current_block_index prev_breathing_data buffer_size Fs;
        
        [y_raw, Fs_raw] = audioread(filepath);
        
        % Check if stereo and convert to mono (first channel)
        if size(y_raw, 2) > 1
            y_mono = y_raw(:, 1);
        else
            y_mono = y_raw;
        end
        
        % Resample the audio to the target internal processing rate (Fs=44100)
        audio_file_data = resample(y_mono, Fs, Fs_raw);
        Fs_audio = Fs; % Fs_audio now holds the resampled rate (44100)

        % 1. INITIALIZE DATA POINTERS
        tic;
        current_block_index = 1;
        prev_breathing_data = zeros(buffer_size, 1);
        
        % 2. SETUP TIMER
        global Update_Interval_Sec;
        timerHandle = timer('ExecutionMode', 'fixedRate', ...
                            'Period', Update_Interval_Sec, ...
                            'TimerFcn', @process_audio_block_file, ...
                            'StopFcn', @stop_cleanup);
        
        % 3. START PROCESSING
        start(timerHandle);
        set(h_status_text, 'String', ['Status: Processing STARTED (' filename ')'], 'ForegroundColor', [0 0.5 0]);

    catch ME
        errordlg(['Error loading/processing file: ' ME.message], 'File Error');
        set(h_status_text, 'String', 'Status: Failed to load file', 'ForegroundColor', [0.8 0 0]);
    end
end

function stop_processing(~, ~)
    global timerHandle h_status_text;
    
    if ~isempty(timerHandle) && isvalid(timerHandle)
        stop(timerHandle);
        delete(timerHandle);
        timerHandle = [];
    end
    
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Processing STOPPED', 'ForegroundColor', [0.8 0 0]);
    end
end

function stop_cleanup(~, ~)
    global h_status_text;
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Processing Finished (File End)', 'ForegroundColor', [0 0 0.8]);
    end
end

function cleanUp(~, ~)
    stop_processing();
    delete(gcf);
end

% --- 3. CORE PROCESSING LOOP (File Input) ---

function process_audio_block_file(~, ~)
    global Fs_audio Fs Update_Interval_Sec Target_Pitch_Hz Target_BPM Target_Range;
    global h_ax_karaoke h_bpm_text h_map_type_dd h_waveform_dd h_ax_spectrogram h_ax_envelope h_ax_fft h_ax_raw;
    global prev_breathing_data buffer_size lpf_breathing;
    global timerHandle audio_file_data current_block_index;

    samples_to_process = Fs * Update_Interval_Sec;
    
    % --- 1. DATA ACQUISITION (Read block from file buffer) ---
    
    end_index = current_block_index + samples_to_process - 1;
    
    if end_index > length(audio_file_data)
        stop(timerHandle); % Stop when file ends
        return;
    end
    
    new_block = audio_file_data(current_block_index : end_index);
    current_block_index = end_index + 1;
    
    % --- 2. DSP & BPM EXTRACTION ---
    
    envelope_full = abs(new_block); 
    breathing_envelope = filter(lpf_breathing.b, lpf_breathing.a, envelope_full);
    prev_breathing_data = [prev_breathing_data(length(new_block)+1:end); breathing_envelope];
    
    Y = fft(prev_breathing_data);
    L = length(prev_breathing_data);
    P1 = abs(Y/L); P1 = P1(1:floor(L/2)+1); P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;

    breathing_range_indices = find(f >= 0.1 & f <= 1);
    
    if isempty(breathing_range_indices)
        current_BPM = Target_BPM;
        max_amplitude = 0.5;
    else
        [~, idx] = max(P1(breathing_range_indices));
        breathing_freq = f(breathing_range_indices(idx));
        current_BPM = breathing_freq * 60;
        max_amplitude = max(P1(breathing_range_indices));
    end

    if isgraphics(h_bpm_text)
        set(h_bpm_text, 'String', sprintf('BPM: %.1f', current_BPM));
    end

    % --- 3. SONIFICATION ENGINE ---
    
    map_type_idx = get(h_map_type_dd, 'Value'); map_type_strings = get(h_map_type_dd, 'String'); map_type = map_type_strings{map_type_idx};
    waveform_idx = get(h_waveform_dd, 'Value'); waveform_type_strings = get(h_waveform_dd, 'String'); waveform_type = waveform_type_strings{waveform_idx};
    
    Pitch_Shift_Semitones = calculate_pitch_shift(current_BPM, Target_BPM, map_type);
    current_Pitch_Hz = Target_Pitch_Hz * 2^(Pitch_Shift_Semitones / 12);
    volume_scale = min(1.0, 0.2 + max_amplitude * 30); 
    
    tone_duration = Update_Interval_Sec;
    t_audio = 0:1/Fs:tone_duration - 1/Fs;

    if strcmp(waveform_type, 'Sine')
        audio_output = sin(2*pi*current_Pitch_Hz*t_audio) * 0.5;
    else 
        start_freq = current_Pitch_Hz * 0.95; end_freq = current_Pitch_Hz * 1.05;
        audio_output = chirp(t_audio, start_freq, tone_duration, end_freq, 'linear') * 0.3; 
    end
    
    player = audioplayer(audio_output * volume_scale, Fs);
    play(player); 

    % --- 4. VISUAL FEEDBACK & PLOTS ---
    
    current_time = toc;
    
    % 4a. Raw Audio Block Plot (Top Left)
    if isgraphics(h_ax_raw)
        cla(h_ax_raw);
        t_raw_block = linspace(current_time - Update_Interval_Sec, current_time, samples_to_process);
        plot(h_ax_raw, t_raw_block, new_block);
        xlim(h_ax_raw, [current_time - Update_Interval_Sec, current_time]);
        title(h_ax_raw, '1. Raw Audio Block (0.5s)');
    end


    % 4b. Filtered FFT Plot (Top Center)
    if isgraphics(h_ax_fft)
        cla(h_ax_fft);
        plot(h_ax_fft, f(2:end), P1(2:end)); 
        xlim(h_ax_fft, [0 1]); 
        
        if ~isempty(breathing_range_indices)
            max_fft_mag = max(P1(breathing_range_indices)) * 1.2;
            ylim(h_ax_fft, [0 max_fft_mag]);
        end

        title(h_ax_fft, '2. Filtered FFT (Breathing Freq)');
    end

    % 4c. Breathing Envelope Plot (Top Right)
    if isgraphics(h_ax_envelope)
        cla(h_ax_envelope);
        time_vector = linspace(current_time - buffer_size/Fs, current_time, buffer_size);
        plot(h_ax_envelope, time_vector, prev_breathing_data);
        xlim(h_ax_envelope, [current_time - 10, current_time]);
        ylim(h_ax_envelope, [0 max(prev_breathing_data)*1.2 + 0.01]);
        title(h_ax_envelope, '3. Filtered Breathing Envelope (10s History)');
    end


    % 4d. Spectrogram Plot (Middle)
    persistent audio_buffer;
    if isempty(audio_buffer)
        audio_buffer = [];
    end
    audio_buffer = [audio_buffer; audio_output'];

    if isgraphics(h_ax_spectrogram) && ~isempty(audio_buffer) && rem(current_time, 2) < Update_Interval_Sec
        spectrogram(audio_buffer, hamming(512), 128, 512, Fs, 'yaxis', 'Parent', h_ax_spectrogram);
        title(h_ax_spectrogram, ['4. Sonification Spectrogram | Time: ' num2str(current_time, '%.1f') 's | BPM: ' num2str(current_BPM, '%.1f')]);
        ylim(h_ax_spectrogram, [0 1]); 
        colormap(h_ax_spectrogram, 'jet');
        colorbar('peer', h_ax_spectrogram);
    end
    
    % 4e. Karaoke Bar Plot (Bottom)
    if isgraphics(h_ax_karaoke)
        cla(h_ax_karaoke);
        
        rect_height = Target_Range(2) - Target_Range(1);
        rectangle('Parent', h_ax_karaoke, 'Position', [0, Target_Range(1), current_time + 5, rect_height], 'FaceColor', [0 0.8 0 0.2], 'EdgeColor', 'none'); 

        if current_BPM > Target_Range(2) || current_BPM < Target_Range(1)
            color = 'r'; 
        else
            color = 'g'; 
        end
        
        bar_height = 0.5;
        rectangle('Parent', h_ax_karaoke, 'Position', [current_time - Update_Interval_Sec, current_BPM - bar_height/2, Update_Interval_Sec, bar_height], 'FaceColor', color, 'EdgeColor', color);
        
        title(h_ax_karaoke, ['5. Karaoke Feedback | Map: ' map_type ' | Wave: ' waveform_type]);
        xlim(h_ax_karaoke, [max(0, current_time - 10), current_time + 2]); 
    end
    
    drawnow limitrate;
end

% --- PITCH MAPPING FUNCTIONS (Identical to previous versions) ---

function shift = calculate_pitch_shift(current_BPM, Target_BPM, map_type)
    
    BPM_deviation = current_BPM - Target_BPM;
    
    if strcmp(map_type, 'Linear')
        shift = BPM_deviation * 1; 
        
    elseif strcmp(map_type, 'Non-Linear')
        max_dev = 10;
        normalized_dev = BPM_deviation / max_dev;
        
        if normalized_dev >= 0
            shift = 12 * (normalized_dev^2); 
        else
            shift = -12 * abs(normalized_dev)^0.5; 
        end
        
    elseif strcmp(map_type, 'Fibonacci')
        
        fib_semitones = [0 1 2 3 5 8 13]; 
        abs_dev = abs(BPM_deviation);
        bpm_steps = [0.5 1.5 2.5 4.5 7.5 12.5]; 
        
        if abs_dev < bpm_steps(1)
            shift_idx = 1;
        else
            shift_idx = find(abs_dev >= bpm_steps, 1, 'last') + 1;
        end
        
        shift_idx = min(shift_idx, length(fib_semitones));
        shift = fib_semitones(shift_idx) * sign(BPM_deviation);
        
    else
        shift = 0;
    end
end

% Function that runs the app upon execution
recorded_file_sonification();