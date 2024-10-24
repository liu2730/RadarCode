clc
close all
clear all

tic
% Set up the parameters for the audio recording
Fs = 48000; % Sampling rate (set according to your audio card's configuration)
nBits = 16; % Bit depth (optional)
nChannels = 1; % Number of channels (mono)
duration = 0.25; % Duration of each recording chunk in seconds
audioRecorder = audiorecorder(Fs, nBits, nChannels);

% Signal parameters
c = 3e8; % Speed of light
fc = 2.481e9; % Carrier frequency
Tp = duration-0.001; % Pulse duration
lambda = c / fc; % Wavelength

UPSAMPLE = 20;
N = Tp * Fs;
upsample = zeros(1, UPSAMPLE * N - N);
velocity_array = (0:N) / 2 * lambda * Fs / N;

% Initialize the figure for real-time plotting
figure;
h1 = imagesc(velocity_array, [], []);
colorbar;
xlim([0 10]);
%clim([-15 0]);
ylabel("Time [s]");
xlabel("Velocity [m/s]");

M = 10000; % Total time duration in seconds (or iterations)

recordblocking(audioRecorder, duration);
signal = getaudiodata(audioRecorder);
signal = transpose(signal(1:N, 1));
signal = signal - mean(signal, "all");

% FFT processing on the recorded audio chunk
fft_amplitudes = fft([signal, upsample]);
noise = fft_amplitudes;
old = fft_amplitudes;

data = 20 * log10(abs(fft_amplitudes));


% Normalize the FFT results
data = data - max(data);


time_array = [toc];
velocities = [0];
tic
for t = 1:M
    time_array(t+1) = time_array(t) + toc;
    tic

    % Record a chunk of audio
    recordblocking(audioRecorder, duration);
    signal = getaudiodata(audioRecorder);
    signal = transpose(signal(1:N, 1));

    % Clutter rejection (DC removal)
    signal = signal - mean(signal, "all");

    % FFT processing on the recorded audio chunk
    fft_amplitudes = fft([signal, upsample]);
    delta_fft = fft_amplitudes;
    
    %delta_fft = delta_fft - noise;
    %old = delta_fft;

    delta_fft_dB = 20 * log10(abs(delta_fft));

    % Normalize the FFT results
    %delta_fft_dB = delta_fft_dB - max(delta_fft_dB);

    % Calculate velocity for this chunk
    x = find(velocity_array > 10);
    [v, idx] = max(delta_fft_dB(1:(x(1) * UPSAMPLE)));

    if v > 30
        velocity = velocity_array(idx) / UPSAMPLE;
        velocities(t) = velocity;
    else
        velocities(t) = 0;
    end

    data = [data; delta_fft_dB];
    
    figure(1);
    lower = max(1, t-30);
    set(h1, 'YData', time_array(lower:t));
    set(h1, 'CData', data(lower:t, :));
    xlim([0 5])



    max_dB = max(data(lower:t, :), [], "all");
    clim([max_dB - 30, max_dB])

    figure(2)
    plot(time_array(1:t), movmean(velocities(1:t), 2), "LineWidth", 2);

    xlabel("Time [s]");
    ylabel("Velocity [m/s]");    
    refreshdata;
    % Pause for a short time to allow the plot to refresh
    drawnow;
end
