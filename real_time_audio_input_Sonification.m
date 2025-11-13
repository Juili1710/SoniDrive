% REAL-TIME BREATHING SONIFICATION APP (LIVE MICROPHONE INPUT)
% Features: Live FFT, Envelope Plotting, 3 Mappings (Linear, Non-Linear, Fibonacci), Karaoke Feedback.

% --- GLOBAL VARIABLES AND INITIAL CONFIGURATION ---
global audioRecorder timerHandle;
global h_fig h_ax_karaoke h_bpm_text h_status_text;
global h_map_type_dd h_waveform_dd;
global h_ax_spectrogram h_ax_envelope h_ax_fft; % Added h_ax_fft
global Fs Update_Interval_Sec Target_Pitch_Hz Target_BPM Target_Range;
global prev_breathing_data buffer_size lpf_breathing;

Fs = 44100;                      % Audio Sampling Rate
Update_Interval_Sec = 0.5;       % How often to process a block of audio (500 ms)
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
function live_sonification_app
    global h_fig h_ax_karaoke h_bpm_text h_status_text;
    global h_map_type_dd h_waveform_dd h_ax_spectrogram h_ax_envelope h_ax_fft Target_Range;

    % --- Create GUI Figure and Layout Adjustments ---
    h_fig = figure('Name', 'Real-Time Sonification (Live Mic Input)', 'Position', [100 100 1000 750], 'CloseRequestFcn', @cleanUp); 
    
    % --- Control Panel (UI elements) ---
    control_x = 20;
    control_width = 150;
    
    uicontrol('Style', 'text', 'String', 'Mapping Type:', 'Position', [control_x 670 control_width 20], 'HorizontalAlignment', 'left');
    h_map_type_dd = uicontrol('Style', 'popupmenu', 'String', {'Linear', 'Non-Linear', 'Fibonacci'}, ...
                              'Position', [control_x 640 control_width 25], 'Value', 1);

    uicontrol('Style', 'text', 'String', 'Waveform:', 'Position', [control_x 610 control_width 20], 'HorizontalAlignment', 'left');
    h_waveform_dd = uicontrol('Style', 'popupmenu', 'String', {'Sine', 'Chirp'}, ...
                              'Position', [control_x 580 control_width 25], 'Value', 1);
    
    % BPM Value Display
    h_bpm_text = uicontrol('Style', 'text', 'String', 'BPM: --', ...
                           'Position', [control_x 520 control_width 30], 'FontSize', 16, 'FontWeight', 'bold', 'ForegroundColor', [0 0 0.8], 'HorizontalAlignment', 'left');
    
    % Start/Stop Buttons
    uicontrol('Style', 'pushbutton', 'String', 'START RECORDING', 'Position', [control_x 430 control_width 50], 'Callback', @start_sonification, 'FontSize', 14, 'BackgroundColor', [0.7 1 0.7]);
    uicontrol('Style', 'pushbutton', 'String', 'STOP', 'Position', [control_x 370 control_width 50], 'Callback', @stop_sonification, 'FontSize', 14, 'BackgroundColor', [1 0.7 0.7]);
    
    h_status_text = uicontrol('Style', 'text', 'String', 'Status: Stopped', 'Position', [control_x 340 control_width 20], 'FontSize', 10, 'ForegroundColor', [0.8 0 0], 'HorizontalAlignment', 'center');

    % --- Plotting Axes (Manual Positions for Clean Layout) ---
    plot_start_x = 0.2; 
    plot_width_half = 0.36;
    plot_width_full = 0.75;
    plot_height = 0.25;

    % 1. Raw Audio FFT (Top Left)
    h_ax_fft = axes('Parent', h_fig, 'Position', [plot_start_x 0.7 0.35 plot_height]); 
    title(h_ax_fft, '1. Raw Audio FFT (0.5s Block)');
    xlabel(h_ax_fft, 'Frequency (Hz)');
    
    % 2. Breathing Envelope (Top Right)
    h_ax_envelope = axes('Parent', h_fig, 'Position', [plot_start_x + plot_width_half 0.7 0.35 plot_height]); 
    title(h_ax_envelope, '2. Filtered Breathing Envelope (10s History)');
    xlabel(h_ax_envelope, 'Time (s)');

    % 3. Sonification Spectrogram (Middle, Full Width)
    h_ax_spectrogram = axes('Parent', h_fig, 'Position', [plot_start_x 0.38 plot_width_full 0.28]); 
    title(h_ax_spectrogram, '3. Sonification Spectrogram (Output Pitch)');
    
    % 4. Karaoke Feedback (Bottom, Full Width)
    h_ax_karaoke = axes('Parent', h_fig, 'Position', [plot_start_x 0.05 plot_width_full 0.28]); 
    title(h_ax_karaoke, ['4. Karaoke Feedback | Target: ' num2str(Target_Range(1)) '-' num2str(Target_Range(2)) ' BPM']);
    ylabel(h_ax_karaoke, 'BPM Level');
    ylim(h_ax_karaoke, [5 25]);
    xlabel(h_ax_karaoke, 'Time (s)'); 
    grid(h_ax_karaoke, 'on');

end

% --- 2. START AND STOP FUNCTIONS ---

function start_sonification(~, ~)
    global audioRecorder timerHandle h_status_text prev_breathing_data;
    global Fs Update_Interval_Sec buffer_size;
    
    % Initialization Guard
    if ~isempty(timerHandle) && isvalid(timerHandle) && strcmp(timerHandle.Running, 'on')
        set(h_status_text, 'String', 'Status: Already Recording', 'ForegroundColor', [0.8 0.8 0]);
        return;
    end

    try
        % 1. INITIALIZE AUDIO RECORDER
        tic; % Start tic before starting the timer loop
        audioRecorder = audiorecorder(Fs, 16, 1);
        
        % Pre-allocate buffer for breathing envelope data (10s)
        prev_breathing_data = zeros(buffer_size, 1);
        
        % 2. SETUP TIMER
        timerHandle = timer('ExecutionMode', 'fixedRate', ...
                            'Period', Update_Interval_Sec, ...
                            'TimerFcn', @process_audio_block, ...
                            'StopFcn', @stop_cleanup);
        
        % Start recording and timer
        record(audioRecorder);
        start(timerHandle);

        set(h_status_text, 'String', 'Status: Recording & Sonifying STARTED', 'ForegroundColor', [0 0.5 0]);
        
    catch ME
        errordlg(['Error initializing audiorecorder. Details: ' ME.message], 'Audio Error');
        if isgraphics(h_status_text)
            set(h_status_text, 'String', 'Status: Failed to Start', 'ForegroundColor', [0.8 0 0]);
        end
    end
end

function stop_sonification(~, ~)
    global audioRecorder timerHandle h_status_text;
    
    if ~isempty(audioRecorder) && isvalid(audioRecorder)
        stop(audioRecorder);
        delete(audioRecorder);
        audioRecorder = [];
    end
    
    if ~isempty(timerHandle) && isvalid(timerHandle)
        stop(timerHandle);
        delete(timerHandle);
        timerHandle = [];
    end
    
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Recording STOPPED', 'ForegroundColor', [0.8 0 0]);
    end
end

function stop_cleanup(~, ~)
    global h_status_text;
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Finished', 'ForegroundColor', [0 0 0.8]);
    end
end

function cleanUp(~, ~)
    stop_sonification();
    delete(gcf);
end

% --- 3. CORE PROCESSING LOOP ---

function process_audio_block(~, ~)
    global audioRecorder Fs Update_Interval_Sec Target_Pitch_Hz Target_BPM Target_Range;
    global h_ax_karaoke h_bpm_text h_map_type_dd h_waveform_dd h_ax_spectrogram h_ax_envelope h_ax_fft;
    global prev_breathing_data buffer_size lpf_breathing;
    
    % --- 1. DATA ACQUISITION & DSP ---
    
    % Read the new block of audio (data length = Fs * Update_Interval_Sec)
    audio_data_raw = getaudiodata(audioRecorder);
    stop_sample = length(audio_data_raw);
    start_sample = max(1, stop_sample - Fs * Update_Interval_Sec + 1);
    new_block = audio_data_raw(start_sample:stop_sample);
    
    if isempty(new_block) || length(new_block) < Fs * Update_Interval_Sec
        return; % Not enough data yet
    end
    
    % Calculate Envelope (proxy for breathing depth)
    envelope_full = abs(new_block); 
    
    % Filter the envelope (Bandpass 0.1Hz - 1Hz)
    breathing_envelope = filter(lpf_breathing.b, lpf_breathing.a, envelope_full);
    
    % Update the history buffer for stable BPM calculation
    prev_breathing_data = [prev_breathing_data(length(new_block)+1:end); breathing_envelope];
    
    % --- 2. BPM EXTRACTION (using FFT on the history buffer) ---
    
    Y = fft(prev_breathing_data);
    L = length(prev_breathing_data);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L; % Frequency vector
    
    % Find the peak in the breathing frequency range (0.1 Hz to 1 Hz)
    breathing_range_indices = find(f >= 0.1 & f <= 1);
    
    if isempty(breathing_range_indices)
        current_BPM = Target_BPM; % Default if no peak found
        max_amplitude = 0.5;
    else
        [max_power, idx] = max(P1(breathing_range_indices));
        breathing_freq = f(breathing_range_indices(idx));
        current_BPM = breathing_freq * 60; % Convert Hz to BPM (x 60)
        max_amplitude = max(P1(breathing_range_indices)); % Max FFT amplitude for Volume control
    end

    if isgraphics(h_bpm_text)
        set(h_bpm_text, 'String', sprintf('BPM: %.1f', current_BPM));
    end

    % --- 3. SONIFICATION ENGINE (PITCH, VOLUME, WAVEFORM) ---
    
    map_type_idx = get(h_map_type_dd, 'Value');
    map_type_strings = get(h_map_type_dd, 'String');
    map_type = map_type_strings{map_type_idx};
    
    waveform_idx = get(h_waveform_dd, 'Value');
    waveform_type_strings = get(h_waveform_dd, 'String');
    waveform_type = waveform_type_strings{waveform_idx};
    
    % PITCH CALCULATION
    Pitch_Shift_Semitones = calculate_pitch_shift(current_BPM, Target_BPM, map_type);
    current_Pitch_Hz = Target_Pitch_Hz * 2^(Pitch_Shift_Semitones / 12);
    
    % VOLUME CALCULATION (Based on Breathing Depth/Amplitude)
    volume_scale = min(1.0, 0.2 + max_amplitude * 30); 
    
    % WAVEFORM GENERATION 
    tone_duration = Update_Interval_Sec;
    t_audio = 0:1/Fs:tone_duration - 1/Fs;

    if strcmp(waveform_type, 'Sine')
        audio_output = sin(2*pi*current_Pitch_Hz*t_audio) * 0.5;
    else % Chirp
        start_freq = current_Pitch_Hz * 0.95;
        end_freq = current_Pitch_Hz * 1.05;
        audio_output = chirp(t_audio, start_freq, tone_duration, end_freq, 'linear') * 0.3; 
    end
    
    player = audioplayer(audio_output * volume_scale, Fs);
    play(player); 

    % --- 4. VISUAL FEEDBACK & PLOTS ---
    
    current_time = toc;
    
    % 4a. Raw Audio FFT Plot (Top Left)
    if isgraphics(h_ax_fft)
        cla(h_ax_fft);
        % FFT of just the new block (raw audio)
        Y_raw = fft(new_block);
        L_raw = length(new_block);
        P2_raw = abs(Y_raw/L_raw);
        P1_raw = P2_raw(1:floor(L_raw/2)+1);
        P1_raw(2:end-1) = 2*P1_raw(2:end-1);
        f_raw = Fs*(0:(L_raw/2))/L_raw;

        plot(h_ax_fft, f_raw, P1_raw);
        xlim(h_ax_fft, [0 Fs/2]);
        ylim(h_ax_fft, [0 max(P1_raw) * 1.2]);
        title(h_ax_fft, '1. Raw Audio FFT (0.5s Block)');
    end

    % 4b. Breathing Envelope Plot (Top Right)
    if isgraphics(h_ax_envelope)
        cla(h_ax_envelope);
        time_vector = linspace(current_time - buffer_size/Fs, current_time, buffer_size);
        plot(h_ax_envelope, time_vector, prev_breathing_data);
        xlim(h_ax_envelope, [current_time - 10, current_time]);
        ylim(h_ax_envelope, [0 max(prev_breathing_data)*1.2 + 0.01]);
        title(h_ax_envelope, '2. Filtered Breathing Envelope (10s History)');
    end


    % 4c. Spectrogram Plot (Middle)
    persistent audio_buffer;
    if isempty(audio_buffer)
        audio_buffer = [];
    end
    audio_buffer = [audio_buffer; audio_output'];

    if isgraphics(h_ax_spectrogram) && ~isempty(audio_buffer) && rem(current_time, 2) < Update_Interval_Sec
        spectrogram(audio_buffer, hamming(512), 128, 512, Fs, 'yaxis', 'Parent', h_ax_spectrogram);
        title(h_ax_spectrogram, ['3. Sonification Spectrogram | Time: ' num2str(current_time, '%.1f') 's | BPM: ' num2str(current_BPM, '%.1f')]);
        ylim(h_ax_spectrogram, [0 1]); 
        colormap(h_ax_spectrogram, 'jet');
        colorbar('peer', h_ax_spectrogram);
    end
    
    % 4d. Karaoke Bar Plot (Bottom)
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
        
        title(h_ax_karaoke, ['4. Karaoke Feedback | Map: ' map_type ' | Wave: ' waveform_type]);
        xlim(h_ax_karaoke, [max(0, current_time - 10), current_time + 2]); 
    end
    
    drawnow limitrate;
end

% --- PITCH MAPPING FUNCTIONS ---
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
live_sonification_app();