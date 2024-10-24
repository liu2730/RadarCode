clc
close all
clear all

[signal, Fs] = audioread("CW2_5.wav");
signal = transpose(signal(1:1543500, 1));
%signal = transpose(signal);

c = 3e8;
fc = 2.481e9;
Tp = 0.2;
lambda = c / fc;

% clutter rejection
signal = signal - mean(signal, "all");

tot_samples = length(signal);
M = floor(tot_samples / Fs / Tp);
N = Tp*Fs;

UPSAMPLE = 10;
upsample = zeros(1,UPSAMPLE*N-N);
fft_amplitudes = zeros(M, N*UPSAMPLE);

for idx=0:(M-1)
    row = signal((1+idx*N):((idx+1)*N));
    upsampled = [row, upsample];
    fft_amplitudes(idx+1, :) = fft(upsampled);
end

delta_fft = zeros(M, N*UPSAMPLE);
delta_fft(1, :) = fft_amplitudes(1, :);
for idx=2:M
    delta_fft(idx, :) = fft_amplitudes(idx, :) - fft_amplitudes(idx-1, :);
end

delta_fft_dB = 20*log10(abs(delta_fft));

% Norm 1
%delta_fft_dB(:, :) = delta_fft_dB(:, :) - max(delta_fft_dB(:, :), [], "all");

% Norm 2
for idx=1:M
   delta_fft_dB(idx, :) = delta_fft_dB(idx, :) - max(delta_fft_dB(idx, :), [], "all");
end

velocity_array = (0:N) / 2 * lambda * Fs / N ;
time_array = (1:M)*Tp;


fastTarget = zeros(1, M);
slowTarget = zeros(1, M);

for i=1:M
    row = delta_fft_dB(i, :);
    [v, idx] = max(row);

    range_out = find(velocity_array > 15);

    [v1, idx] = max(row(1:(range_out(1) * UPSAMPLE)));
    t1 = velocity_array(idx) / UPSAMPLE;

    skip = 4*UPSAMPLE;
    row(max(1,idx-skip):min(length(row),idx+skip)) = -300;

    [v, idx] = max(row(1:(range_out(1) * UPSAMPLE)));
    %idx = min(idx, length(velocity_array));
    t2 = velocity_array(idx) / UPSAMPLE;
    
    if v < -5
        slowTarget(i) = min(t1, t2);
    else
        slowTarget(i) = min(t1, t2);
        fastTarget(i) = max(t1, t2);
    end

    %max_slow_idx = find(velocity_array > 0.75);

    %[v, idx] = max(row(1:(max_slow_idx(1) * UPSAMPLE)));
    %slowTarget(i) = velocity_array(min(idx, length(velocity_array))) / UPSAMPLE;
end

% Plot
imagesc(velocity_array, time_array, delta_fft_dB)
colorbar;
xlim([0 10])
clim([-15 0])
yticks(0:1:30)
ylabel("Time [s]")
xlabel("Velocity [m/s]")

figure;
plot(time_array, fastTarget);
hold on;
plot(time_array, smooth(slowTarget));
legend("Target 1", "Target 2")
title("1 Target");
xlabel("Time [s]");
ylabel("Velocity [m/s]");

