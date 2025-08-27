function G = compute_G(u, Para, Wq, Wq1, Grid)

v = fftn(u);

ux = real(ifft2(-1i*Grid.kxx.*v));         
uy = real(ifft2(-1i*Grid.kyy.*v));

G = sum(3*Para.epsilon^2 * u.^2 .* (ux.^2 + uy.^2) + ...
    0.5 * Wq1.^2 - Para.eta2 * Wq, 'all') * prod(Grid.d);

end