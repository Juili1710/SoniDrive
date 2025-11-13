% SIMULATED BREATHING DATA GENERATOR
% Generates complex, realistic breathing signal data that can be used for
% offline testing, comparing against the live app, or creating static plots
% for the conference paper.

% --- CONFIGURATION ---
Fs = 44100;         % Sampling rate of the final output
Duration_Sec = 60;  % Total simulation length (60 seconds)
t = 0:1/Fs:Duration_Sec - 1/Fs;

% --- BREATHING STATES (BPM) ---
% Define BPM over time to simulate a driver transitioning through states.
% Normal (14 BPM) -> Fatigue (8 BPM) -> Stress (20 BPM)

% 0-20s: Normal Breathing (14 BPM, clean)
t_normal = t(t <= 20);
f_normal = 14 / 60; 
amplitude_normal = 1;
noise_normal = 0.1 * randn(size(t_normal));
signal_normal = amplitude_normal * sin(2*pi*f_normal*t_normal) + noise_normal;

% 20-40s: Fatigue/Drowsiness (8 BPM, shallow, low amplitude, some noise)
t_fatigue = t(t > 20 & t <= 40);
f_fatigue = 8 / 60;
amplitude_fatigue = 0.5; % Shallow depth
noise_fatigue = 0.2 * randn(size(t_fatigue));
signal_fatigue = amplitude_fatigue * sin(2*pi*f_fatigue*t_fatigue) + noise_fatigue;

% 40-60s: Stress/Alert (20 BPM, rapid, slightly erratic)
t_stress = t(t > 40);
f_stress = 20 / 60;
amplitude_stress = 1.2;
noise_stress = 0.05 * randn(size(t_stress));
% Add a second, slightly higher frequency component for "erratic" feeling
signal_stress = amplitude_stress * sin(2*pi*f_stress*t_stress) + ...
                0.3 * sin(2*pi*(f_stress*1.2)*t_stress) + noise_stress;


% --- COMBINE SIGNALS ---
simulated_breathing_signal = [signal_normal, signal_fatigue, signal_stress];
time_vector = t;

% --- VISUALIZATION ---
figure('Name', 'Simulated Breathing Data States');
plot(time_vector, simulated_breathing_signal);
title('Simulated Driver Breathing Signal (Amplitude Envelope)');
xlabel('Time (s)');
ylabel('Normalized Amplitude');
grid on;

% --- SAVE DATA (Optional) ---
% save('simulated_breathing_data.mat', 'simulated_breathing_signal', 'Fs', 'time_vector');
% disp('Data saved to simulated_breathing_data.mat');

disp('Simulation data generated. The graph shows the transition from Normal -> Fatigue -> Stress.');