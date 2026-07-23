function Grid = generate_Grid(L, N, dim)

Grid.L = L;
Grid.N = N;

Grid.d = Grid.L ./ Grid.N;

Grid.x = cell(dim, 1);
k = cell(dim, 1);
k_even = cell(dim, 1);

for i = 1:dim
    Grid.x{i} = Grid.L(i)*(1:Grid.N(i))' / Grid.N(i) - Grid.L(i)/2;
    k{i} = [ 0:Grid.N(i)/2-1, 0.0, -Grid.N(i)/2+1:-1]' / (Grid.L(i)/pi/2);
    % even-order (Laplacian/biharmonic) operators have no sign ambiguity at
    % the Nyquist mode, so unlike k{i} above it should NOT be zeroed out
    k_even{i} = [ 0:Grid.N(i)/2-1, Grid.N(i)/2, -Grid.N(i)/2+1:-1]' / (Grid.L(i)/pi/2);
end

if dim == 1
    Grid.x = Grid.x{1};
    Grid.kxx = k{1};

    Grid.k = sqrt(Grid.kxx.^2);
    k2_full = k_even{1}.^2;
end

if dim == 1
    Grid.x = Grid.x{1};
    Grid.k = Grid.k{1};
end

if dim == 2
    [Grid.xx, Grid.yy] = meshgrid(Grid.x{1}, Grid.x{2});
    [Grid.kxx, Grid.kyy] = meshgrid(k{1}, k{2});
    [kxx_even, kyy_even] = meshgrid(k_even{1}, k_even{2});

    Grid.k = sqrt(Grid.kxx.^2 + Grid.kyy.^2);
    k2_full = kxx_even.^2 + kyy_even.^2;

end

if dim == 3

    [Grid.xx, Grid.yy, Grid.zz] = meshgrid(Grid.x{1}, Grid.x{2}, Grid.x{3});
    [Grid.kxx, Grid.kyy, Grid.kzz] = meshgrid(k{1}, k{2}, k{3});
    [kxx_even, kyy_even, kzz_even] = meshgrid(k_even{1}, k_even{2}, k_even{3});

    Grid.k = sqrt(Grid.kxx.^2 + Grid.kyy.^2 + Grid.kzz.^2);
    k2_full = kxx_even.^2 + kyy_even.^2 + kzz_even.^2;
end

Grid.inv_k = 1./(Grid.k.^2);      Grid.inv_k(Grid.k == 0) = 1;

Grid.k2 = k2_full;
Grid.k4 = (Grid.k2).^2;
Grid.k6 = (Grid.k2).^3;

end