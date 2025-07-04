function Grid = generate_Grid(L, N, dim, model)

Grid.L = L;
Grid.N = N;

Grid.d = Grid.L ./ Grid.N;

x = cell(dim, 1);
k = cell(dim, 1);

for i = 1:dim
    x{i} = Grid.L(i)*(1:Grid.N(i))' / Grid.N(i) - Grid.L(i)/2;  
    k{i} = [ 0:Grid.N(i)/2-1, 0.0, -Grid.N(i)/2+1:-1]' / (Grid.L(i)/pi/2);
end

if dim == 2
    [Grid.xx, Grid.yy] = meshgrid(x{1}, x{2});
    [Grid.kxx, Grid.kyy] = meshgrid(k{1}, k{2});

    Grid.k = sqrt(Grid.kxx.^2 + Grid.kyy.^2);

end

if dim == 3

    [Grid.xx, Grid.yy, Grid.zz] = meshgrid(x{1}, x{2}, x{3});
    [Grid.kxx, Grid.kyy, Grid.kzz] = meshgrid(k{1}, k{2}, k{3});

    Grid.k = sqrt(Grid.kxx.^2 + Grid.kyy.^2 + Grid.kzz.^2);
end

Grid.inv_k = 1./(Grid.k.^2);      Grid.inv_k(Grid.k == 0) = 1;

Grid.k2 = (Grid.k).^2;
Grid.k4 = (Grid.k2).^2;

end