function [h, Dh, tt] = hermf(N, M, tm)
% P. Flandrin & J. Xiao, 2005
% Computes a set of orthonormal Hermite functions 
% INPUT
%   N  : Number of points (must be odd).
%   M  : Maximum order.
%   tm : Half-time support (>= 6 recommended).
% OUTPUT
%   h  : Hermite functions (N x M).
%   Dh : Derivatives (N x M).
%   tt : Time vector (N x 1).

dt = 2 * tm / (N - 1); 
tt = linspace(-tm, tm, N); 
g = exp(-tt.^2 / 2); 

P = zeros(M + 1, N);
P(1, :) = ones(1, N); 
P(2, :) = 2 * tt; 
for k = 3:M + 1 
    P(k, :) = 2 * tt .* P(k - 1, :) - 2 * (k - 2) * P(k - 2, :); 
end

H = zeros(M + 1, N);
for k = 1:M + 1    
    H(k, :) = P(k, :) .* g / sqrt(sqrt(pi) * 2^(k - 1) * gamma(k) / dt); 
end 
h = H(1:M, :); 

Dh = zeros(size(h));
for k = 1:M   
    Dh(k, :) = (tt .* H(k, :) - sqrt(2 * k) * H(k + 1, :)) * dt; 
end 
h = h'; Dh = Dh'; tt = tt';

end
