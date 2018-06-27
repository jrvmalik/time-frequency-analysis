function imageSQ(M, y, h, l)
% IMAGESQ Displays the time-frequency representation M as a black and white
% image. A clipping effect is used to normalize M.
% INPUT
%      M    : Time-frequency representation to be visualized.
%      y    : Frequency axis.
%      h    : Range above which to clip, a number in (0, 1], e.g. 0.995
%      l    : Discard numbers below this quantile (default = 0).
% Written by John Malik on 2018.6.22, john.malik@duke.edu.

switch nargin
    case 1
        y = 1;
        h = 1;
        l = 0;
    case 2
        h = 1;
        l = 0;
    case 3
        l = 0;
end

q = quantile(M(:), [h, l]);
M(M > q(1)) = q(1);
M(M < q(2)) = 0;

h = get(groot, 'defaultAxesTickLabelInterpreter');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

imagesc(1, y, M)
axis xy
colormap(1 - gray)

set(groot, 'defaultAxesTickLabelInterpreter', h);

end
