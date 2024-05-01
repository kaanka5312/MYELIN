function [bold, time] = hrfconv(v_E, time, dt)
%% Convolve excitatory firing rate with HRF function

c_time = 0:dt:10/dt; % Time of kernel (10 seconds)
tau_c = 1.25; % Timescale (in ms)
d = 2.25; % Delay (in ms)

H = (c_time - d) .* exp(-(c_time - d) / tau_c) / tau_c^2;

% Badass FFT convolution
Nconv = length(time) + length(H) - 1;
fH = fft(H, Nconv);
fH = fH ./ max(fH);
fV = fft(v_E', Nconv)';
for i = 1:360
    fBOLD(i,:) = fV(i,:) .* fH;
end
half_wav = floor( length(H)/2 )+1;
BOLD = ifft(fBOLD')';
bold = real(BOLD(:, half_wav-1:end-half_wav));
% Downsample BOLD to be equal to fMRI
% BOLD = downsample(BOLD', 2000)'; BOLD = BOLD(:, 1:300);
end