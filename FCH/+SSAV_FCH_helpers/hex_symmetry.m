function [S, info] = hex_symmetry(Grid)
% Projection onto fields with the full symmetry of a hexagonal lattice, on
% one exact hex cell: width Lx = d, height Ly = d*sqrt(3), lattice points
% at the origin and at the cell corners (as in hex_lattice_state and
% SpotStripeTest). Returns a function handle S with S(v) the projection of
% v (an Ny x Nx array, or a stack Ny x Nx x K).
%
% The symmetry group is D6 (order 12): rotations by multiples of 60 degrees
% about a lattice point, and reflections. In Fourier space the cell's
% modes are k = 2*pi*(i/Lx, j/Ly), and a rotation by 60 degrees maps
% mode (i, j) to ((i - j)/2, (3i + j)/2), an exact integer map when i + j
% is even. Modes with i + j odd change sign under the half-cell shift
% (Lx/2, Ly/2), which is also a symmetry of the lattice, so they are zero.
% S therefore:
%   - zeroes the modes with i + j odd,
%   - replaces each remaining Fourier coefficient by the average over the
%     12 images of its mode under D6,
%   - zeroes modes whose orbit does not fit on the grid (it would need a
%     mode with |i| >= Nx/2 or |j| >= Ny/2); since D6 preserves |k|, these
%     are exactly the modes outside a circle of radius about
%     pi*min(Nx/Lx, Ny/Ly), i.e. only the highest wavenumbers.
% This is an exact, orthogonal projection (S(S(v)) = S(v), and S is
% symmetric in the grid inner product). Coefficients are taken relative to
% the lattice point at x = y = 0, not the first grid point.
%
% Patterns with this symmetry: round spots in a hex lattice (dark or
% light), honeycombs, and the uniform state. Ovals, dumbbells, stripes and
% wavy stripes are removed by S, so continuation restricted by S cannot
% switch onto them.
%
% The discretization (pointwise products on the rectangular grid) is not
% exactly D6-equivariant, so mu(S(v)) can have a small non-symmetric part
% (aliasing); restricted solvers should project the residual with S too.
%
% info: dim (dimension of the symmetric subspace, including the mean
% mode), kept (logical Ny x Nx mask of kept modes), M (the projection as a
% sparse matrix on Fourier coefficients relative to x = y = 0).

Lx = Grid.L(1);  Ly = Grid.L(2);
Nx = Grid.N(1);  Ny = Grid.N(2);
if abs(Ly/Lx - sqrt(3)) > 1e-10
    error('hex_symmetry: needs a hex cell with Ly/Lx = sqrt(3); got %.12g', Ly/Lx);
end
if mod(Nx, 2) || mod(Ny, 2)
    error('hex_symmetry: needs even grid sizes; got %d x %d', Nx, Ny);
end

% signed mode indices in fft2 order; columns are x (i), rows are y (j)
sx = [0:Nx/2-1, Nx/2, -Nx/2+1:-1];
sy = [0:Ny/2-1, Ny/2, -Ny/2+1:-1];
[I, J] = meshgrid(sx, sy);
n = Nx*Ny;

inside = @(i, j) abs(i) < Nx/2 & abs(j) < Ny/2;     % Nyquist modes excluded
pos = @(i, j) mod(i, Nx)*Ny + mod(j, Ny) + 1;        % linear index in the Ny x Nx array

% the 12 group elements as integer maps of (i, j): R^r and (reflection in x) * R^r
R = [1/2, -1/2; 3/2, 1/2];
maps = cell(12, 1);
Rr = eye(2);
for r = 0:5
    maps{r+1} = Rr;
    maps{r+7} = [-1, 0; 0, 1] * Rr;
    Rr = R * Rr;
end

ok = mod(I + J, 2) == 0 & inside(I, J);
img_i = zeros(n, 12);
img_j = zeros(n, 12);
for g = 1:12
    A = maps{g};
    ig = round(A(1,1)*I + A(1,2)*J);
    jg = round(A(2,1)*I + A(2,2)*J);
    img_i(:, g) = ig(:);
    img_j(:, g) = jg(:);
    ok = ok & inside(ig, jg);
end

rows = find(ok(:));
nr = numel(rows);
cols = pos(img_i(rows, :), img_j(rows, :));
M = sparse(repmat(rows, 12, 1), cols(:), 1/12, n, n);

% phase factors: coefficients relative to x = y = 0 instead of the first grid point
x1 = Grid.L(1)/Grid.N(1) - Grid.L(1)/2;
y1 = Grid.L(2)/Grid.N(2) - Grid.L(2)/2;
ph = exp(-1i*(2*pi*I/Lx*x1 + 2*pi*J/Ly*y1));

S = @(v) apply(v, M, ph, Ny, Nx);

info.dim = round(full(sum(diag(M))));
info.kept = ok;
info.n_modes = nr;
info.M = M;

end

function w = apply(v, M, ph, Ny, Nx)
sz = size(v);
C = fft2(reshape(v, Ny, Nx, [])) .* ph;
C = reshape(M * reshape(C, Ny*Nx, []), Ny, Nx, []);
w = reshape(real(ifft2(C ./ ph)), sz);
end
