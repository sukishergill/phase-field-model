function L = spectral_lap(v, Grid)
% Spectral Laplacian with Grid.k2 (Nyquist mode kept). Applies to a single
% field or to a stack of fields: columns of an N x B array in 1D
% (size(Grid.k2, 2) == 1), pages of an Ny x Nx x B array in 2D.
if size(Grid.k2, 2) == 1
    L = real(ifft(-Grid.k2 .* fft(v)));
else
    L = real(ifft2(-Grid.k2 .* fft2(v)));
end
end
