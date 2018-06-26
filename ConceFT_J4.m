function [cft, frequency] = ConceFT_J4(x, Fs, hlength, iter, hop, n, hf, lf)
% Computes the concentrated time-freq representation of the signal x.
% This serial implementation uses less memory but is slower than ConceFT_J3.
% Choose iter = 1 to compute the SST of the signal x.
% Choosing iter = 0 is not supported.
% INPUT
%    x      :  Signal (x should be a column vector).
%    Fs     :  Sampling rate of x.
%    hlength:  Window length.
%    iter   :  Number of SST computations.
%    hop    :  Calculate the fft every hop samples, starting at 1.
%    n      :  Number of pixels in the frequency axis.
%    lf     :  Crop output to display only frequencies larger than lf.
%    hf     :  Crop output to display only frequencies less than hf.
% OUTPUT
%    cft    :  The concentrated time-freq representation of the signal x.
% This implementation requires the function hermf.m.
% Written by John Malik on 2017.4.28, john.malik@duke.edu.

switch nargin
    case 7
        lf = 0;
    case 6
        hf = inf;
        lf = 0;
    case 5
        hf = inf;
        lf = 0;
        n = pow2(nextpow2(length(x))) / 2 + 1;
    case 4
        hf = inf;
        lf = 0;
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
    case 3
        hf = inf;
        lf = 0;
        iter = 1; 
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
    case 2
        error('Select a sampling rate.')
    case 1
        error('Select a window length.')
    case 0
        Fs = 200;
        x = 2 * mod(1e-2:1e-2:1e2, 1) - 1;
        x = x + random('Normal', zeros(size(x)), 0.33);
        hlength = 1001;
        iter = 10;
        lf = 1;
        hop = 40;
        n = 8000;
        hf = 5;
        disp('Testing code on a 2 Hz sawtooth wave.')
end

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

% orthogonal hermite windows
bandwidth = 6;
funx = 4;
[H, dH] = hermf(hlength, funx, bandwidth);

% non-negative frequency axis, cropped
frequency = Fs / 2 * linspace(0, 1, n)'; 
frequency(frequency > hf) = [];
neta = length(frequency);

% iterate
cft = zeros(neta, tcol);
for ii = 1:max(1, iter)

if iter == 0
    error('Choose a number of iterations greater than zero.');
else disp(['Iteration ' num2str(ii) ': Calculating SST.']);
end

if iter < 2
    h = H(:, 1); dh = dH(:, 1);
else
    z = randn(funx, 1) + 1i * randn(funx, 1);
    z = z ./ norm(z);
    h = conj(H * z);
    dh = conj(dH * z);
end

% STFT
tfr = zeros(neta, tcol); 
tfr2 = zeros(neta, tcol);
for icol = 1:tcol,
    
    tmp = zeros(N, 1);
    tmp2 = zeros(N, 1);
    
    ti = t(icol); 
    tau = -min([n - 1, Lh, ti - 1]):min([n - 1, Lh, NN - ti]);
    indices = rem(N + tau, N) + 1; 
    rSig = x(ti + tau, 1);
    
    tmp(indices) = rSig .* h(Lh + 1 + tau);
    tmp = fft(tmp);
    tfr(:, icol) = tmp(1:neta);
    
    tmp2(indices) = rSig .* dh(Lh + 1 + tau);
    tmp2 = fft(tmp2);
    tfr2(:, icol) = tmp2(1:neta);
    
end

% normalize
tfr = tfr / sqrt(N);
tfr2 = tfr2 / sqrt(N);

% omega operator
omega = -inf(neta, tcol);
ups = abs(tfr) > 0;
omega(ups) = round(N * imag(tfr2(ups) ./ tfr(ups) / (2 * pi)));

% mapped out of range
index = repmat((1:neta)', [1, tcol]);
omega = index - omega;
id = omega < 1 | omega > neta;
omega(id) = index(id);
sst = tfr; sst(id) = 0;

% reassignment
id = bsxfun(@plus, omega, neta * (0:tcol - 1));
cft = cft + reshape(accumarray(id(:), sst(:), [tcol * neta, 1]), [neta, tcol]);

end

% average
cft = cft / max(1, iter);

% crop output
cft(frequency < lf) = [];
frequency(frequency < lf) = [];

end