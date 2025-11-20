function E = compute_E(u_fft, w, Para, Grid)

% u_fft = fft2(u_fft);
ux = real(ifft2(-1i*Grid.kxx.*u_fft));         
uy = real(ifft2(-1i*Grid.kyy.*u_fft));

du = ux.^2 + uy.^2;

E = sum(sum(0.5*Para.epsilon^4 * (real(ifftn(-Grid.k.^2 .* u_fft))).^2 - ...
    (Para.tau^2/3 + Para.eta1/2 + 1)*Para.epsilon^2*du))*prod(Grid.d)...
    + w^2 - Para.B;

end