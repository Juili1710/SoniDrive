% SIMULATED_SONIFICATION_APP.M
% Simulated Comparative Sonification Study App for Breathing Rate (BPM)
%
% This script simulates a driver's breathing signal transitioning through
% Normal, Fatigue, and Stress states. It allows the user to compare
% Linear vs. Non-Linear vs. Fibonacci pitch mapping and Sine vs. Chirp waveforms in real-time.

% --- GLOBAL VARIABLES AND INITIAL CONFIGURATION ---
global timerHandle;
global h_fig h_ax_karaoke h_bpm_text h_status_text;
global h_map_type_dd h_waveform_dd;
global h_ax_spectrogram;

% CORE CONFIGURATION PARAMETERS - Declared globally and initialized immediately
global Fs Update_Interval_Sec Target_Pitch_Hz Target_BPM total_sim_time simulation_time;

Fs = 44100;                 % Audio Sampling Rate
Update_Interval_Sec = 0.5;  % How often to update the simulation
Target_Pitch_Hz = 440;      % A4 note for target BPM (14 BPM)
Target_BPM = 14;            % Center of the Normal range
total_sim_time = 30;        % Total length of simulation
simulation_time = 0;        % Current time tracker

% --- 1. SETUP AND GUI CREATION ---
function simulated_sonification_app
    global total_sim_time; 
    global h_fig h_ax_karaoke h_bpm_text h_status_text;
    global h_map_type_dd h_waveform_dd;
    global h_ax_spectrogram;
    
    % --- Create GUI Figure and Layout Adjustments ---
    h_fig = figure('Name', 'Simulated Comparative Sonification Study', 'Position', [100 100 1000 700], 'CloseRequestFcn', @cleanUp);
    
    % --- Control Panel (UI elements) ---
    control_x = 20;
    control_width = 150;
    
    uicontrol('Style', 'text', 'String', 'Mapping Type:', 'Position', [control_x 630 control_width 20], 'HorizontalAlignment', 'left');
    % ADDED TRUE FIBONACCI OPTION
    h_map_type_dd = uicontrol('Style', 'popupmenu', 'String', {'Linear', 'Non-Linear', 'Fibonacci'}, ...
                              'Position', [control_x 600 control_width 25], 'Value', 1);
    uicontrol('Style', 'text', 'String', 'Waveform:', 'Position', [control_x 570 control_width 20], 'HorizontalAlignment', 'left');
    h_waveform_dd = uicontrol('Style', 'popupmenu', 'String', {'Sine', 'Chirp'}, ...
                              'Position', [control_x 540 control_width 25], 'Value', 1);
    uicontrol('Style', 'text', 'String', 'BPM:', 'Position', [control_x 500 control_width 20], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    h_bpm_text = uicontrol('Style', 'text', 'String', '--', 'Position', [control_x 470 control_width 20], 'FontSize', 14, 'FontWeight', 'bold', 'ForegroundColor', [0 0 0.8], 'HorizontalAlignment', 'left');
    
    % Start Button
    uicontrol('Style', 'pushbutton', 'String', 'START SIMULATION', 'Position', [control_x 360 control_width 50], 'Callback', @start_simulation, 'FontSize', 14, 'BackgroundColor', [0.7 1 0.7]);
    
    % Stop Button
    uicontrol('Style', 'pushbutton', 'String', 'STOP', 'Position', [control_x 300 control_width 50], 'Callback', @stop_simulation, 'FontSize', 14, 'BackgroundColor', [1 0.7 0.7]);
    
    h_status_text = uicontrol('Style', 'text', 'String', 'Status: Stopped', 'Position', [control_x 270 control_width 20], 'FontSize', 10, 'ForegroundColor', [0.8 0 0], 'HorizontalAlignment', 'center');
    
    % --- Plotting Axes ---
    plot_start_x = 0.2; % Normalized position for plots
    
    % 1. Spectrogram (Top Plot)
    h_ax_spectrogram = subplot(2, 1, 1, 'Parent', h_fig, 'Position', [plot_start_x 0.5 0.75 0.45]);
    title(h_ax_spectrogram, 'Real-Time Sonification Spectrogram');
    
    % 2. Karaoke Feedback (Bottom Plot)
    h_ax_karaoke = subplot(2, 1, 2, 'Parent', h_fig, 'Position', [plot_start_x 0.05 0.75 0.4]);
    title(h_ax_karaoke, 'Karaoke Feedback | Target: 12-15 BPM');
    ylabel(h_ax_karaoke, 'BPM Level');
    ylim(h_ax_karaoke, [5 25]);
    xlabel(h_ax_karaoke, 'Time (s)');
    grid(h_ax_karaoke, 'on');
    
    cla(h_ax_spectrogram);
    cla(h_ax_karaoke);
end

% --- 2. START AND STOP FUNCTIONS ---
function start_simulation(~, ~)
    global timerHandle simulation_time h_status_text;
    global Update_Interval_Sec;
    if ~isempty(timerHandle) && isvalid(timerHandle) && strcmp(timerHandle.Running, 'on')
        set(h_status_text, 'String', 'Status: Already Running', 'ForegroundColor', [0.8 0.8 0]);
        return;
    end
    
    global h_ax_karaoke h_ax_spectrogram;
    
    cla(h_ax_karaoke);
    cla(h_ax_spectrogram);
    
    simulation_time = 0;
    tic;
    
    timerHandle = timer('ExecutionMode', 'fixedRate', ...
                        'Period', Update_Interval_Sec, ...
                        'TimerFcn', @process_simulation_block, ...
                        'StopFcn', @stop_cleanup);
    
    start(timerHandle);
    set(h_status_text, 'String', 'Status: Running Simulation', 'ForegroundColor', [0 0.5 0]);
end

function stop_simulation(~, ~)
    global timerHandle h_status_text;
    
    if ~isempty(timerHandle) && isvalid(timerHandle)
        stop(timerHandle);
        delete(timerHandle);
        timerHandle = [];
    end
    
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Stopped', 'ForegroundColor', [0.8 0 0]);
    end
end

function stop_cleanup(~, ~)
    global h_status_text;
    if isgraphics(h_status_text)
        set(h_status_text, 'String', 'Status: Simulation Complete', 'ForegroundColor', [0 0 0.8]);
    end
end

function cleanUp(~, ~)
    stop_simulation();
    delete(gcf);
end

% --- 3. CORE SIMULATION AND PROCESSING LOOP ---
function process_simulation_block(~, ~)
    global Fs Update_Interval_Sec Target_Pitch_Hz total_sim_time simulation_time;
    global h_ax_karaoke h_bpm_text h_map_type_dd h_waveform_dd h_ax_spectrogram;
    global timerHandle Target_BPM;
    
    % End condition
    if simulation_time >= total_sim_time
        stop(timerHandle);
        return;
    end
    
    t_end = simulation_time + Update_Interval_Sec;
    
    % 1. Simulate Breathing Signal (Breathing Rate/BPM)
    
    if t_end <= 10
        current_BPM = 14 + 1.5*sin(t_end*0.5); 
        state = 'Normal (13-15 BPM)';
    elseif t_end <= 20
        current_BPM = 8 + 1*rand(); 
        state = 'Fatigue (7-9 BPM)';
    else
        current_BPM = 22 + 2*sin(t_end*2); 
        state = 'Stress (20-24 BPM)';
    end
    
    set(h_bpm_text, 'String', sprintf('%.1f', current_BPM));
    
    % 2. SONIFICATION ENGINE (Mapping & Waveform Selection)
    
    map_type_idx = get(h_map_type_dd, 'Value');
    map_type_strings = get(h_map_type_dd, 'String');
    map_type = map_type_strings{map_type_idx}; % Assumes cell array structure from uicontrol
    
    waveform_idx = get(h_waveform_dd, 'Value');
    waveform_type_strings = get(h_waveform_dd, 'String');
    waveform_type = waveform_type_strings{waveform_idx}; % Assumes cell array structure from uicontrol
    
    % PITCH CALCULATION - Uses selected map_type
    Pitch_Shift_Semitones = calculate_pitch_shift(current_BPM, Target_BPM, map_type);
    Pitch_Shift_Semitones = max(-12, min(12, Pitch_Shift_Semitones));
    current_Pitch_Hz = Target_Pitch_Hz * 2^(Pitch_Shift_Semitones / 12);
    
    % Generate the sonified tone
    tone_duration = Update_Interval_Sec;
    t_audio = 0:1/Fs:tone_duration - 1/Fs;
    
    % WAVEFORM GENERATION - Uses selected waveform_type
    if strcmp(waveform_type, 'Sine')
        audio_output = sin(2*pi*current_Pitch_Hz*t_audio) * 0.5;
    elseif strcmp(waveform_type, 'Chirp')
        start_freq = current_Pitch_Hz * 0.9;
        end_freq = current_Pitch_Hz * 1.1;
        audio_output = chirp(t_audio, start_freq, tone_duration, end_freq, 'linear') * 0.3;
    else
        audio_output = zeros(size(t_audio));
    end
    
    player = audioplayer(audio_output, Fs);
    play(player);
    
    % --- 3. VISUAL FEEDBACK & PLOTS ---
    
    current_time = toc;
    
    % 3a. Karaoke Bar Plot (Bottom Plot)
    
    if isgraphics(h_ax_karaoke)
        cla(h_ax_karaoke);
        
        % 1. Draw Static Target Box (12-15 BPM)
        rectangle('Parent', h_ax_karaoke, 'Position', [0, 12, total_sim_time, 3], 'FaceColor', [0 0.8 0 0.2], 'EdgeColor', 'none');
        
        % 2. Determine color and position for Dynamic Bar
        if current_BPM > 15 || current_BPM < 12
            color = 'r'; % Red for deviation
        else
            color = 'g'; % Green for on-target
        end
        
        bar_height = 0.5;
        
        % 3. Draw Dynamic Bar (The current BPM marker)
        rectangle('Parent', h_ax_karaoke, 'Position', [current_time - Update_Interval_Sec, current_BPM - bar_height/2, Update_Interval_Sec, bar_height], 'FaceColor', color, 'EdgeColor', color);
        
        % 4. Set Limits and Titles
        title(h_ax_karaoke, ['Karaoke Feedback | Mapping: ' map_type ' | Waveform: ' waveform_type]);
        xlabel(h_ax_karaoke, 'Time (s)');
        xlim(h_ax_karaoke, [max(0, current_time - 10), current_time + 2]);
        ylim(h_ax_karaoke, [5 25]);
        grid(h_ax_karaoke, 'on');
    end
    
    % 3b. Spectrogram Plot (Top Plot)
    persistent audio_buffer;
    if isempty(audio_buffer) || simulation_time == 0
        audio_buffer = [];
    end
    audio_buffer = [audio_buffer; audio_output'];
    if isgraphics(h_ax_spectrogram) && ~isempty(audio_buffer) && rem(current_time, 2) < Update_Interval_Sec
        spectrogram(audio_buffer, hamming(512), 128, 512, Fs, 'yaxis', 'Parent', h_ax_spectrogram);
        title(h_ax_spectrogram, ['Sonification Spectrogram | Time: ' num2str(current_time, '%.1f') 's | State: ' state]);
        ylim(h_ax_spectrogram, [0 5]);
        colormap(h_ax_spectrogram, 'jet');
        colorbar('peer', h_ax_spectrogram);
    end
    
    drawnow limitrate;
    
    % Advance time
    simulation_time = t_end;
end

% --- PITCH MAPPING FUNCTIONS ---
function shift = calculate_pitch_shift(current_BPM, Target_BPM, map_type)
    
    BPM_deviation = current_BPM - Target_BPM;
    
    if strcmp(map_type, 'Linear')
        % Linear: Direct proportionality (1 semitone per 1 BPM deviation)
        shift = BPM_deviation * 1; 
        
    elseif strcmp(map_type, 'Non-Linear')
        % Non-Linear (Quadratic/Root): Accelerating deviation at extremes
        max_dev = 10;
        normalized_dev = BPM_deviation / max_dev;
        
        if normalized_dev >= 0
            % Stress (High BPM): Quadratic shift (accelerates rapidly)
            shift = 12 * (normalized_dev^2); 
        else
            % Fatigue (Low BPM): Square Root shift (accelerates less sharply)
            shift = -12 * abs(normalized_dev)^0.5; 
        end
        
    elseif strcmp(map_type, 'Fibonacci')
        % True Fibonacci Series: Discrete steps based on harmonic intervals
        
        % Fibonacci intervals (in semitones): 0, 1, 2, 3, 5, 8, 13
        fib_semitones = [0 1 2 3 5 8 13]; 
        
        % Determine the absolute deviation from the target
        abs_dev = abs(BPM_deviation);
        
        % Define the BPM thresholds that trigger a jump to the next Fibonacci interval
        % e.g., 0.5 BPM deviation = 0 semitones; 1.5 BPM dev = 2 semitones
        bpm_thresholds = [0.5 1.5 2.5 4.5 7.5 12.5]; 
        
        % Find the step index
        if abs_dev < bpm_thresholds(1)
            shift_idx = 1; % Index 1 is 0 semitone shift
        else
            % Find the largest threshold that the current deviation exceeds
            shift_idx = find(abs_dev >= bpm_thresholds, 1, 'last') + 1;
        end
        
        % Clamp index to prevent array out of bounds
        shift_idx = min(shift_idx, length(fib_semitones));
        
        % Apply the shift and restore the sign (+ for stress, - for fatigue)
        shift = fib_semitones(shift_idx) * sign(BPM_deviation);
        
    else
        shift = 0; % Default
    end
end

% Function that runs the app upon execution
simulated_sonification_app();