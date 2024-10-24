clc
close all
clear

[signal, Fs] = audioread("FMCW1_5.wav");
c = 3e8;
Tp = 20e-3;
N = floor(Tp*Fs);

sync = -sign(signal(:, 2));
signal(:, 1) = movmean(signal(:, 1), 5);
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

% 2 Pulse MTI
% mti_matrix(1, :) = upchirps(1, :);
% for idx=2:M
%    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :);
% end

% 3 Pulse MTI
%mti_matrix(1, :) = upchirps(1, :);
%mti_matrix(2, :) = upchirps(2, :);
%for idx=3:M
%    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :) - (upchirps(idx-1, :) - upchirps(idx-2, :));
%end

% No MTI
 mti_matrix = upchirps;

UPSAMPLE = 8;
upsample = zeros(1,UPSAMPLE*N-N);

ifft_amplitudes = zeros(M, N*UPSAMPLE);
fft_amplitudes = zeros(M, N*UPSAMPLE);
for idx=1:M
    upsampled = [mti_matrix(idx, :), upsample];
    ifft_amplitudes(idx, :) = ifft(upsampled);
    fft_amplitudes(idx, :) = fft(upsampled);
end

ifft_dB = 20*log10(abs(ifft_amplitudes));
ifft_dB = ifft_dB - max(ifft_dB,[],"all");

delta_f = 0.087e9;
range_array = (0:N) / 2 * c / delta_f;
time_array = (0:M-1)*Tp*2;

ranges = zeros(1, M);
for i=1:M
    row = ifft_dB(i, :);
    x = find(range_array > 25);
    [v, idx] = max(row(1:(x*UPSAMPLE)));
    ranges(i) = range_array(min(idx, length(range_array))) / UPSAMPLE; 
end


figure;
% Image plot
imagesc(range_array, time_array, ifft_dB)

axis xy;  % Correct axis orientation
title('3-Pulse MTI');
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;

xlim([0 40])
clim([-50 0])

% Range plot
figure;
%plot(time_array, smooth(ranges));
plot(time_array, ranges);
ranges = smooth(ranges);
xlabel("Time [s]");
ylabel("Range [m]")

dx = 5
velocities_fd = zeros(1, M);
velocities_fd((1+dx):(M-dx)) = (ranges((1+dx):(M-dx)) - ranges(1:(M-2*dx))) ./ (2*dx*Tp);
hold on;
plot(time_array, movmean(velocities_fd, 20));

%% Dangerous
M2 = 768
L = 12;
X = M2/L

fc = 2.431e9;
times_sw = zeros(1, X);
velocities_sw_up = zeros(1, X);
for idx=0:(X-1)
    mat = fft_amplitudes(1+idx*L:((idx+1)*L), :);
    for i=1:L
        fft_mat = fft(mat(:, i));
        [x, row] = max(fft_mat(1:ceil(L/2)));
    
        times_sw(idx+1) = time_array(1+idx*L) + L*Tp;
        velocities_sw_up(idx+1) = (row-1)*c/fc/2/Tp/L;
    end
end

plot(linspace(0, max(times_sw), 100), ones(1, 100)*c/fc/4/Tp, "k--");
plot(linspace(0, max(times_sw), 100), -ones(1, 100)*c/fc/4/Tp, "k--");

grid on;
%legend("Range", "Velocity FDM", "Vel slow")

% ---------------------------------------------------------- DOWNCHIRPS

sync = sign(signal(:, 2));
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

% 2 Pulse MTI
% mti_matrix(1, :) = upchirps(1, :);
% for idx=2:M
%    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :);
% end

% 3 Pulse MTI
%mti_matrix(1, :) = upchirps(1, :);
%mti_matrix(2, :) = upchirps(2, :);
%for idx=3:M
%    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :) - (upchirps(idx-1, :) - upchirps(idx-2, :));
%end

% No MTI
 mti_matrix = upchirps;

upsample = zeros(1,UPSAMPLE*N-N);

ifft_amplitudes = zeros(M, N*UPSAMPLE);
fft_amplitudes = zeros(M, N*UPSAMPLE);

for idx=1:M
    upsampled = [mti_matrix(idx, :), upsample];
    ifft_amplitudes(idx, :) = ifft(upsampled);
    fft_amplitudes(idx, :) = fft(upsampled);
end

ifft_dB = 20*log10(abs(ifft_amplitudes));
ifft_dB = ifft_dB - max(ifft_dB,[],"all");

delta_f = 0.087e9;
range_array = (0:N) / 2 * c / delta_f;
time_array = (0:M-1)*Tp*2;

ranges = zeros(1, M);
for i=1:M
    row = ifft_dB(i, :);
    x = find(range_array > 25);
    [v, idx] = max(row(1:(x*UPSAMPLE)));
    ranges(i) = range_array(min(idx, length(range_array))) / UPSAMPLE; 
end

fc = 2.431e9;
times_sw = zeros(1, X);
velocities_sw_down = zeros(1, X);
for idx=0:(X-1)
    mat = fft_amplitudes(1+idx*L:((idx+1)*L), :);
    for i=1:L
        fft_mat = fft(mat(:, i));
        [x, row] = max(fft_mat(1:ceil(L/2)));
    
        times_sw(idx+1) = time_array(1+idx*L) + L*Tp;
        velocities_sw_down(idx+1) = (row-1)*c/fc/2/Tp/L;
    end
end

plot(times_sw, velocities_sw_up-velocities_sw_down)

