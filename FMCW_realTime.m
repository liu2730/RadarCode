clc
close all
clear all

tic
% Set up the parameters for the audio recording
Fs = 48000; % Sampling rate (set according to your audio card's configuration)
nBits = 16; % Bit depth (optional)
nChannels = 2; % Number of channels (stereo)
duration = 0.2; % Duration of each recording chunk in seconds
audioRecorder = audiorecorder(Fs, nBits, nChannels);

% Signal parameters
c = 3e8; % Speed of light
%fc = 2.481e9; % Carrier frequency
fc = 2.431e9; % Carrier frequency
Tp = 20e-3; % Pulse duration
lambda = c / fc; % Wavelength

UPSAMPLE = 10;
N = floor(Tp * Fs);
upsample = zeros(1, UPSAMPLE * N - N);

% Initialize the figure for real-time plotting
figure;

delta_f = 0.087e9;
range_array = (0:N) / 2 * c / delta_f;
h1 = imagesc(range_array, [], []);
colorbar;
xlim([0 10]);
%clim([-15 0]);
ylabel("Time [s]");
xlabel("Range [m]");

T = 10000; % Total time duration in seconds (or iterations)

time_array = [toc];
velocities = [0];
data = [ones(1, N*UPSAMPLE)] .* -100;
ranges = [0];
ranges_time = [0];

tic
for t = 1:T

    % Record a chunk of audio
    recordblocking(audioRecorder, duration);
    signal = getaudiodata(audioRecorder);

    sync = -sign(signal(:, 2));
    range_signal = signal(:, 1) .* sync;

    M = 0;
    idx = 1;
    while idx <= length(range_signal)
        while idx <= length(range_signal) && sync(idx) < 0
            idx = idx+1;
        end
        
        if idx >= length(range_signal)-N
            break
        end

        while idx <= length(range_signal) && sync(idx) >= 0
            idx = idx+1;
        end
        
        if idx >= length(range_signal)-N
            break
        end
    
        start_idx = idx;
        end_idx = idx+N;
        M = M+1;
        upchirps(M, :) = range_signal(start_idx:end_idx-1);
    
        while idx <= length(range_signal) && sync(idx) >= 0
            idx = idx+1;
        end
    end


    % MS clutter rejection
    for col=1:N
        upchirps(:, col) = upchirps(:, col) - mean(upchirps(:, col), "all");
    end


    % 3 Pulse MTI
    mti_matrix(1, :) = upchirps(1, :);
    mti_matrix(2, :) = upchirps(2, :);
    for idx=3:M
        mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :) - (upchirps(idx-1, :) - upchirps(idx-2, :));
    end

    mti_matrix = upchirps; 

    upsample = zeros(1,UPSAMPLE*N-N);

    %for idx=1:M
    %    upchirps(idx, :) = upchirps(idx, :) - max(upchirps(idx, :), [], "all");
    %end
    
    ifft_amplitudes = zeros(M, N*UPSAMPLE);
    for idx=1:M
        upsampled = [mti_matrix(idx, :), upsample];
        ifft_amplitudes(idx, :) = ifft(upsampled);
    end
    
    ifft_dB = 20*log10(abs(ifft_amplitudes));
    ifft_dB = ifft_dB - max(ifft_dB,[],"all");
    
    delta_f = 0.087e9;
    time_array = [time_array, toc + (1:M)*Tp];
    
    data = [data; ifft_dB];

    figure(1);
    lower = max(1, length(time_array)-50);
    set(h1, 'YData', time_array(lower:length(time_array)));
    set(h1, 'CData', data(lower:length(time_array), :));
    xlim([0 15])

    max_dB = max(data(lower:length(time_array), :), [], "all");
    clim([max_dB - 10, max_dB])

    
    sz = size(data);
    lastrow = data(sz(1), :);
    
    maxidx = find(range_array > 15);

    [x, idx] = max(lastrow(1:(maxidx * UPSAMPLE)));
    ranges(t+1) = range_array(min(idx, length(range_array))) / UPSAMPLE;
    ranges_time(t+1) = toc;

    figure(2)
    plot(ranges_time, movmean(ranges, 3), "LineWidth", 2);

    %xlabel("Time [s]");
    %ylabel("Velocity [m/s]");    
    refreshdata;
    % Pause for a short time to allow the plot to refresh
    drawnow;
end
