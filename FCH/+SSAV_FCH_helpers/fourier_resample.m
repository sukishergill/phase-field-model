function v = fourier_resample(u, Nnew)
% Spectral (trigonometric) interpolation of a periodic field from the
% generate_Grid grid x_j = L*j/N - L/2 (j = 1..N) onto the same-L grid with
% Nnew points. 1D column (Nnew scalar) or 2D Ny x Nx (Nnew = [Nx_new, Ny_new],
% generate_Grid's N ordering). Refinement only (Nnew >= N, Nnew a multiple
% of N), so the old grid points are a subset of the new ones.
%
% interpft samples from the first old point x_1 with the new spacing; on
% the new generate_Grid grid x_1 is point number r = Nnew/N, so the result
% is shifted by r - 1 to line up.

if size(u, 2) == 1
    N = size(u, 1);
    r = Nnew / N;
    assert(r >= 1 && r == round(r), 'fourier_resample: Nnew must be a multiple of N');
    v = circshift(interpft(u, Nnew), r - 1);
else
    [Ny, Nx] = size(u);
    rx = Nnew(1) / Nx;
    ry = Nnew(2) / Ny;
    assert(rx >= 1 && rx == round(rx) && ry >= 1 && ry == round(ry), ...
        'fourier_resample: Nnew must be multiples of the old size');
    v = interpft(u, Nnew(2), 1);
    v = interpft(v, Nnew(1), 2);
    v = circshift(v, [ry - 1, rx - 1]);
end
end
