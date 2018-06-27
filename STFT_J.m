function [tfr, frequency] = STFT_J(x, Fs, hlength, hop, n, hf, lf)
% Computes the short-time Fourier transform of the signal x.
% INPUT
%    x      :  Signal (x should be a column vector).
%    Fs     :  Sampling rate of x.
%    hlength:  Window length (in samples).
%    hop    :  Calculate the fft every hop samples, starting at 1.
%    n      :  Number of pixels in the frequency axis.
%    lf     :  Crop output to display only frequencies larger than lf.
%    hf     :  Crop output to display only frequencies less than hf.
% OUTPUT
%    tfr    :  The STFT of the signal x.
% Written by John Malik on 2018.6.22, john.malik@duke.edu.

switch nargin
    case 6 
        lf = 0;
    case 5
        hf = inf;
        lf = 0;
    case 4
        n = pow2(nextpow2(length(x))) / 2 + 1;
        hf = inf;
        lf = 0;
    case 3
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
        hf = inf;
        lf = 0;
    case 2
        error('Select a window length.')
    case 1
        error('Select a sampling rate.')
    case 0
        Fs = 200;
        x = 2 * mod(1e-2:1e-2:1e2, 1) - 1;
        hlength = 1001;
        lf = 1;
        hop = 40;
        n = 8000;
        hf = 12;
        disp('Testing code on a 2 Hz sawtooth wave.')
end

% window bandwidth
sigma = 0.15;

% organize input
x = x(:);
if any(isnan(x))
    x = interp1(find(~isnan(x)), x(~isnan(x)), 1:length(x), 'pchip', 'extrap')';
end

% time (samples)
NN = length(x);
t = 1:hop:NN;
[~, tcol] = size(t);

% N-point fft
n = n + 1 - rem(n, 2);
N = 2 * (n - 1);

% make sure window length is odd
hlength = hlength + 1 - rem(hlength, 2);
Lh = (hlength - 1) / 2;

% gaussian window
ex = linspace(-0.5, 0.5, hlength)';
h = exp(-ex.^2 / 2 / sigma^2);

% non-negative frequency axis, cropped
frequency = Fs / 2 * linspace(0, 1, n)'; 
neta = sum(frequency <= hf);
frequency = frequency(frequency <= hf);

% perform convolution
tfr = zeros(N, tcol); 
for icol = 1:tcol
    ti = t(icol); 
    tau = -min([n - 1, Lh, ti - 1]):min([n - 1, Lh, NN - ti]);
    indices = rem(N + tau, N) + 1; 
    rSig = x(ti + tau, 1);
    tfr(indices, icol) = rSig .* h(Lh + 1 + tau);
end

% Fourier transform
tfr = fft(tfr);
tfr = tfr(1:neta, :);

% crop output
tfr(frequency < lf, :) = [];
frequency(frequency < lf) = [];

end