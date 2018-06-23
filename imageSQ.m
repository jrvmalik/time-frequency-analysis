function imageSQ(ytic, M, d)

q = quantile(M(:), d);
M(M > q) = q;

imagesc(1, ytic, M)
axis xy
colormap(1 - gray)

end
