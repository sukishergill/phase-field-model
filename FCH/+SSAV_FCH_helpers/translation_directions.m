function [T, gram_eigs] = translation_directions(u, Grid, rank_tol)
% Orthonormal basis (columns of T, numel(u) rows) for span{du/dx, du/dy},
% dropping directions in which u is (numerically) invariant: stripes have
% one, a constant state has none. Using the eigenvectors of the Gram matrix
% also handles tilted stripes, whose invariant direction is a combination
% of x and y.
if nargin < 3, rank_tol = 1e-8; end
U = fftn(u);
if size(u, 2) == 1
    D = real(ifft(1i*Grid.kxx .* U));
else
    ux = real(ifft2(1i*Grid.kxx .* U));
    uy = real(ifft2(1i*Grid.kyy .* U));
    D = [ux(:), uy(:)];
end
n = size(D, 1);
[V, L] = eig((D'*D) / n);
[gram_eigs, order] = sort(diag(L), 'descend');
V = V(:, order);
% relative threshold, with an absolute floor so a constant state
% (derivatives at roundoff level) gets no directions at all
keep = gram_eigs > max(rank_tol*gram_eigs(1), 1e-12);
T = D * V(:, keep);
T = T ./ vecnorm(T);
end
