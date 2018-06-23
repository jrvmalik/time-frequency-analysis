function imageSQ(ytic, M, d)
% IMAGESQ Displays the time-frequency representation M as a black and white
% image. A clipping effect is used to normalize M.
% INPUT
%      M    : Time-frequency representation to be visualized.
%      ytic : Frequency axis.
%      d    : Range above which to clip, a number in (0, 1], e.g. 0.995
% Adapted by John Malik on 2018.6.22, john.malik@duke.edu.

q = quantile(M(:), d);
M(M > q) = q;

h = get(groot, 'defaultAxesTickLabelInterpreter');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

imagesc(1, ytic, M)
axis xy
colormap(1 - gray)

set(groot, 'defaultAxesTickLabelInterpreter', h);

end
