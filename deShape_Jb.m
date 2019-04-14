function [deshape, ceps, mask, tfr, frequency, quefrency] = deShape_Jb(x, Fs, hlength, hf, gamma, hop, n, lf, ths)
% Computes the de-shape synchrosqueezing transform of the signal x.
% INPUT
%    x      :  Signal (x should be a column vector).
%    Fs     :  Sampling rate of x.
%    hlength:  Window length (in samples).
%    hop    :  Calculate the fft every hop samples, starting at 1.
%    n      :  Number of pixels in the frequency axis.
%    lf     :  Crop output to display only frequencies larger than lf.
%    hf     :  Least upper bound on fundamental frequency of all components.
%    gamma  :  Cepstral power.
%    ths    :  Fraction of values to reassign. 
% OUTPUT
%    deshape:  The de-shape SST of the signal x. 
%    ceps   :  The short-time cepstral transform of the signal x.
%    mask   :  The inverted short-time cepstral transform of the signal x.
%    tfr    :  The STFT of the signal x.
% Written by John Malik on 2018.6.22, john.malik@duke.edu.

switch nargin
    case 8
        ths = 1;
    case 7
        ths = 1;
        lf = 0;
    case 6 
        n = pow2(nextpow2(length(x))) / 2 + 1;
        lf = 0;
        ths = 1;
    case 5
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
        lf = 0;
        ths = 1;
    case 4
        gamma = 0.2;
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
        lf = 0;
        ths = 1;
    case 3
        error('Select an upper bound for the fundamental frequency of all components.')
    case 2
        error('Select a window length.')
    case 1
        error('Select a sampling rate.')
    case 0
        Fs = 200;
        x = 2 * mod(1e-2:1e-2:1e2, 1) - 1;
        hlength = 1001;
        lf = 1;
        gamma = 0.2;
        hop = 40;
        n = 8000;
        hf = 5;
        ths = 0.1;
        disp('Testing code on a 2 Hz sawtooth wave.')
end

% cepstrum resampling factor
alpha = 20;

% window bandwidth
sigma = 0.15;

% do reassignment
squeeze_flag = 1;

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

% gaussian window and its derivative
ex = linspace(-0.5, 0.5, hlength)';
h = exp(-ex.^2 / 2 / sigma^2);
dh = -ex ./ sigma^2 .* h;

% perform convolution
tfr = zeros(N, tcol); 
tfr2 = zeros(N, tcol);
for icol = 1:tcol
    ti = t(icol); 
    tau = -min([n - 1, Lh, ti - 1]):min([n - 1, Lh, NN - ti]);
    indices = rem(N + tau, N) + 1; 
    rSig = x(ti + tau, 1);
    rSig = rSig - mean(rSig);
    tfr(indices, icol) = rSig .* h(Lh + 1 + tau);
    tfr2(indices, icol) = rSig .* dh(Lh + 1 + tau);
end

% Fourier transform
tfr = fft(tfr);
tfr = tfr(1:n, :);

% non-negative frequency axis
frequency = Fs / 2 * linspace(0, 1, n)'; 

if squeeze_flag
    tfr2 = fft(tfr2);
    tfr2 = tfr2(1:n, :);
end

% short-time cepstral transform and quefrency axis
ceps = ifft(abs(tfr).^gamma, N, 1); 
ceps = real(ceps(1:n, :));
quefrency = (0:n-1)' / Fs;

% remove envelope
ceps(quefrency < 1 / hf, :) = 0;

% step of unknown purpose
flo = tfr;
flo(abs(tfr) < real(fft(ceps)).^(1 / gamma)) = 0;

% upsample cepstrum by a factor of alpha
uceps = interp1(1:n, ceps, 1:1/alpha:n);
uquefrency = (0:1/alpha:n-1)' / Fs;

% inverted short-time cepstral transform
mask = zeros(size(tfr));
df = frequency(end) / (2 * n);
for i = 1:length(frequency)
    qj = find(uquefrency > 1 / (frequency(i) + df) & ...
        uquefrency <= 1 / (frequency(i) - df));
    if isempty(qj)
        continue
    end
    mask(i, :) = sum(diag(1 ./ qj) * uceps(qj, :), 1);
end

% force negative entries to be zero
ceps = max(ceps, 0);
mask = max(mask, 0);

% crop output
u = frequency <= hf;
tfr = tfr(u, :);
mask = mask(u, :);
frequency = frequency(u);
flo = flo(u, :);

% crop output (ceps)
ceps = ceps(quefrency <= 2 / lf, :);
quefrency = quefrency(quefrency <= 2 / lf);

% deshape short-time Fourier transform
deshape = flo .* mask;

if ~squeeze_flag
    return
end

% crop
tfr2 = tfr2(u, :);

% reassignment threshold
ths = quantile(abs(tfr), (1 - ths));

% omega operator
neta = length(frequency);
omega = -inf(neta, tcol);
ups = ~bsxfun(@le, abs(tfr), ths);
omega(ups) = round(N / hlength * imag(tfr2(ups) ./ tfr(ups) / (2 * pi)));

% mapped out of range
index = repmat((1:neta)', [1, tcol]);
omega = index - omega;
id = omega < 1 | omega > neta | ~ups;
omega(id) = index(id);
deshape(id) = 0;

% reassignment
id = bsxfun(@plus, omega, neta * (0:tcol - 1));
deshape = reshape(accumarray(id(:), deshape(:), [tcol * neta, 1]), [neta, tcol]);

% crop
tfr(frequency < lf, :) = [];
mask(frequency < lf, :) = [];
deshape(frequency < lf, :) = [];
frequency(frequency < lf) = [];

end