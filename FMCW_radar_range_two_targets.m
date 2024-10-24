clc
close all
clear

[signal, Fs] = audioread("FMCW2_8.wav");
plot(signal(:, 1))
figure
plot((1:length(signal(:,2)))/Fs, signal(:, 2))
figure
c = 3e8;
Tp = 20e-3;
N = floor(Tp*Fs);

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
%mti_matrix(1, :) = upchirps(1, :);
%for idx=2:M
%    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :);
%end

% 3 Pulse MTI
mti_matrix(1, :) = upchirps(1, :);
mti_matrix(2, :) = upchirps(2, :);
for idx=3:M
    mti_matrix(idx, :) = upchirps(idx, :) - upchirps(idx-1, :) - (upchirps(idx-1, :) - upchirps(idx-2, :));
end

% No MTI
%mti_matrix = upchirps;

UPSAMPLE = 8;
upsample = zeros(1,UPSAMPLE*N-N);

ifft_amplitudes = zeros(M, N*UPSAMPLE);
for idx=1:M
    upsampled = [mti_matrix(idx, :), upsample];
    ifft_amplitudes(idx, :) = ifft(upsampled);
end

ifft_dB = 20*log10(abs(ifft_amplitudes));
ifft_dB = ifft_dB - max(ifft_dB,[],"all");

delta_f = 0.087e9;
range_array = (0:N) / 2 * c / delta_f;
time_array = (0:M-1)*Tp*2;

target1 = zeros(1, M);
target2 = zeros(1, M);

for i=1:M
    row = ifft_dB(i, :);
    range_out = find(range_array > 20);

    [v, idx] = max(row(1:(range_out(1) * UPSAMPLE)));
    t1 = range_array(idx) / UPSAMPLE;

    skip = 2*UPSAMPLE;
    row(max(1,idx-skip):min(length(row),idx+skip)) = -300;

    [v, idx] = max(row(1:(range_out(1) * UPSAMPLE)));
    idx = min(idx, length(range_array));
    t2 = range_array(idx) / UPSAMPLE;

    target1(i) = min(t1, t2);
    target2(i) = max(t1, t2);
end

% Image plot
imagesc(range_array, time_array, ifft_dB)

axis xy;  % Correct axis orientation
title('3-Pulse MTI');
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;

xlim([0 15])
clim([-40 0])

% Range plot
figure;
plot(time_array, smooth(movmean(target1, 10)));
hold on;
plot(time_array, smooth(movmean(target2, 10)));
legend("Target 1", "Target 2");
xlabel("Time [s]");
ylabel("Range [m]")
