function E = compute_E(u, w, Para, Grid)

u = fft2(u);
ux = real(ifft2(-1i*Grid.kxx.*u));         
uy = real(ifft2(-1i*Grid.kyy.*u));

du = ux.^2 + uy.^2;

E = sum(sum(0.5*Para.epsilon^4 * (real(ifftn(-Grid.k.^2 .* u))).^2 - ...
    (Para.tau^2/3 + Para.eta1/2 + 1)*Para.epsilon^2*du))*prod(Grid.d)...
    + w^2 - Para.B;

end