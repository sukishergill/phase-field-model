function Z = even_basis(sz)
% Orthonormal basis of mean-zero functions even in x and y on the grid:
% one normalized indicator per reflection orbit of grid points, then the
% constant removed. About (N/2+1)^2 columns in 2D.
ny = sz(1);  nx = sz(2);
ix = [nx-1:-1:1, nx]';      % columns, so ix(C)/iy(R) keep the shape of C/R
iy = [ny-1:-1:1, ny]';      % (also when they are vectors, i.e. in 1D)
[C, R] = meshgrid(1:nx, 1:ny);
if nx == 1
    orbit = min(sub2ind(sz, R, C), sub2ind(sz, iy(R), C));
else
    orbit = min(cat(3, sub2ind(sz, R, C), sub2ind(sz, iy(R), C), ...
        sub2ind(sz, R, ix(C)), sub2ind(sz, iy(R), ix(C))), [], 3);
end
[~, ~, lab] = unique(orbit(:));
Z0 = sparse(1:numel(lab), lab, 1);
Z0 = Z0 ./ sqrt(sum(Z0, 1));            % orthonormal orbit indicators
c = full(Z0' * ones(size(Z0, 1), 1));   % constant in orbit coordinates
Z = full(Z0) * null(c');                % remove the mean
end
