clc
close all
clear

%% Signal Input and Parsing Settings
[signal, Fs] = audioread("SAR_object.wav");
c = 3e8;
Tp = 20e-3;
N = floor(Tp*Fs);

sync = -signal(:, 2);
range_signal = -signal(:, 1);

Trp = 1; % s
Nrp = Trp * Fs;

fc = 2.481e9;
lambda = c / fc;

figure(1)
plot(range_signal)
xlabel('Data Sample Number')
ylabel('Amplitude')
title('Sample Down-converted Data')
%xlim([0 12E3])

figure(2)
plot(sync)
xlabel('Data Sample Number')
ylabel('Amplitude')
title('Sync Data')
%xlim([0 12E3])

%A = 1;
%B = 500000;
%plot(A:B, sync(A:B))

%% Position based data parsing
N_POSITIONS = 0;
%DISCARD_SAMPLES = [2, 7, 9, 17, 20, 21, 25, 28, 31, 32, 33, 35, ...
%    38, 39, 41, 43, 44, 46, 48, 49, 52, 61, 74, 79, 80, 83]; % starts from one for landscape
%DISCARD_SAMPLES = 0; % Use this for test file
DISCARD_SAMPLES = []; % For Object

idx = 1;
pos = 1;

% Graphing the samples over time
% time = (1:length(sync))/Fs;
% figure(3);
% plot(time, sync);
% figure(4);
% plot(time, range_signal);

while idx <= length(range_signal) && abs(sync(idx)) > 0.25
    idx = idx+1;
end
% idx is now at first sync=0

while idx <= length(range_signal) - 2*Nrp
    % skip sync = 0
    while idx <= length(range_signal) && abs(sync(idx)) < 0.25
        idx = idx+1;
    end

    if idx >= length(range_signal)-Nrp
        break
    end
    
    idx = idx+Nrp;
    
    start_idx = idx;
    end_idx = idx+Nrp;
        
    num_zeros = sum(abs(sync(start_idx:end_idx-1)) < 0.25);
    if num_zeros < Nrp/4
        if all(DISCARD_SAMPLES ~= pos)
            N_POSITIONS = N_POSITIONS + 1;
            data_matrix(N_POSITIONS, :) = range_signal(start_idx:end_idx-1);
            sync_matrix(N_POSITIONS, :) = sync(start_idx:end_idx-1);
        end

        pos = pos+1;
        
        idx = idx+Nrp;
    end
        
    while idx <= (length(range_signal) - 10) && (abs(sync(idx)) > 0.25 || abs(sync(idx+10)) > 0.25)
        idx = idx+1;
    end
end

% for i=1:(size(sync_matrix, 1))
%     for j=1:size(sync_matrix, 2) 
%         if (sync_matrix(i, j) > 0.05)
%             sync_matrix(i,j) = 1;
%         else
%             sync_matrix(i , j) = 0;
%         end
%     end
% end

% Shorter Rail Length
% length_remove = 3;
% start_bound = length_remove+1;
% end_bound = N_POSITIONS - length_remove + 1;
% data_matrix = data_matrix(start_bound:end_bound,:);
% sync_matrix = sync_matrix(start_bound:end_bound,:);
% N_POSITIONS = N_POSITIONS - (2 * length_remove) + 1;

figure(3)
plot(data_matrix(3, :))
hold on
plot(sync_matrix(3, :))
xlim([0 6E3])
ylim([-0.5 1.5])
xlabel('Data Sample Number')
ylabel('Amplitude')
title('Position parsed sync and down-converted data')

%% Integration of upchirps
for p=1:N_POSITIONS
    S = data_matrix(p, :);
    Sy = sync_matrix(p, :);
    output = zeros(1, N);

    M = 0;
    idx = 1;
    while idx <= length(S)
        while idx <= length(S) && Sy(idx) < 0
            idx = idx+1;
        end
        
        if idx >= length(S)-N
            break
        end

        start_idx = idx;
        end_idx = idx+N;
        M = M+1;
        output = output + S(start_idx:end_idx-1);

        while idx <= length(S) && Sy(idx) >= 0
            idx = idx+1;
        end
    end

    output = output / M;
    int_data_matrix(p, :) = output;
end

% figure;
% plot(int_data_matrix(4, :))

%% Hilbert transform
x = (-N/4:N/4-1);
hanns_window = (1 + cos(4*pi*x ./ N)) ./ 2;
% figure;
% plot(1:N/2, hanns_window);

for row=1:N_POSITIONS
    fft_row = fft(int_data_matrix(row, :));
    hil_data_matrix(row, :) = ifft(fft_row(1:length(fft_row)/2)) .* hanns_window;
    %hil_row = hilbert(int_data_matrix(row, :);
    
end

fmin = 2.4e9;
fmax = 2.5e9;
Kr = linspace(4*pi/c * fmin, 4*pi/c * fmax, N/2);
Xa = lambda/2 * ((0:N_POSITIONS-1) - N_POSITIONS/2);

figure(4)
imagesc(Kr, Xa, angle(hil_data_matrix));
xlabel('K_{r} (rad/m)', 'FontSize', 18)
ylabel('Synthetic Aperture Position, Xa (m)', 'FontSize', 18)
title('Phase Before Along Track FFT', 'FontSize', 24)
colorbar;

%% Cross-range fft
PAD = 1012;
pad_matrix = zeros(PAD, N/2);
zeropadded_data_matrix = [pad_matrix; hil_data_matrix; pad_matrix];

for col=1:N/2
    crfft_data_matrix(:, col) = fftshift(fft(zeropadded_data_matrix(:, col)));
end


S = size(crfft_data_matrix);
Zpad = S(1);
Kx = linspace(-2*pi/lambda, 2*pi/lambda, Zpad);

figure(5);
imagesc(Kr, Kx, 20*log10(abs(crfft_data_matrix)))
colorbar;
caxis([-27 12])
xlabel('K_{r} (rad/m)', 'FontSize', 18)
ylabel('K_{x} (rad/m)', 'FontSize', 18)
title('Magnitude After Along Track FFT', 'FontSize', 24)

figure(6);
imagesc(Kr, Kx, angle(crfft_data_matrix))
colorbar;
xlabel('K_{r} (rad/m)', 'FontSize', 18)
ylabel('K_{x} (rad/m)', 'FontSize', 18)
title('Phase After Along Track FFT', 'FontSize', 24)

for row=1:Zpad
    Ky (row, :) = sqrt(Kr.^2 - Kx(row)^2);
end

%% Interp data matrix
Kye = linspace(min(Ky, [], "all"), max(Ky, [], "all"), Zpad/2);
for row=1:S(1)
    interp_data_matrix(row, :) = interp1(Ky(row, :), crfft_data_matrix(row, :), Kye);
end

interp_data_matrix(isnan(interp_data_matrix))=1e-30;

figure(7);
imagesc(Kye, Kx, angle(interp_data_matrix))
colorbar;
xlabel('K_{y} (rad/m)', 'FontSize', 18)
ylabel('K_{x} (rad/m)', 'FontSize', 18)
title('Phase After Stolt Interpolation', 'FontSize', 24)
%caxis([-60 -50])

ifft2d_data_matrix = ifft2(interp_data_matrix, 4*Zpad, 2*Zpad);
ifft2d_data_matrix = rot90(ifft2d_data_matrix,3); % -90 deg

%% Final Image

delta_fy = c*(max(Ky, [], "all") - min(Ky, [], "all")) / 4 / pi;
delta_ky = 4*pi*delta_fy/c;

R_max = (Zpad*2)*c/2/delta_fy;
Rail_Rmax = Zpad * 6e-2;

d_range_1 = 1;
d_range_2 = 100;
c_range_1 = -25;
c_range_2 = 25;

d_index_1 = round(size(ifft2d_data_matrix, 1) / R_max * d_range_1);
d_index_2 = floor(round(size(ifft2d_data_matrix, 1) / R_max * d_range_2)*3.5);

c_index_1 = round((size(ifft2d_data_matrix, 2) / Rail_Rmax) * (c_range_1 + (Rail_Rmax/2)));
c_index_2 = round((size(ifft2d_data_matrix, 2) / Rail_Rmax) * (c_range_2 + (Rail_Rmax/2)));

cross_range = linspace(c_range_1, c_range_2, c_index_2-c_index_1);
down_range = linspace(-d_range_2, -d_range_1, d_index_2-d_index_1);

H = size(ifft2d_data_matrix, 1);
trunc_data_matrix = ifft2d_data_matrix((H-d_index_2):(H-d_index_1-1), c_index_1:c_index_2);

% for col=1:size(trunc_data_matrix, 2)
%     trunc_data_matrix(:, col) = trunc_data_matrix(:, col) .* (abs(down_range).^2)';
% end

trunc_data_matrix = 20*log10(abs(trunc_data_matrix));

figure(8);
imagesc(cross_range, down_range, trunc_data_matrix);
colorbar;
caxis([-140 -80])
xlim([-20 20])
xlabel('Cross range (m)', 'FontSize', 18)
ylabel('Down range (m)', 'FontSize', 18)
title('Final Image', 'FontSize', 24)
ylim([-40 -1])

