function E = compute_E(u_fft, w, Para, Grid)

% u_fft = fft2(u_fft);
if size(u_fft, 2) == 1
    ux = real(ifftn(-1i*Grid.k.*u_fft));
    du = ux.^2;

else
    ux = real(ifft2(-1i*Grid.kxx.*u_fft));   
    uy = real(ifft2(-1i*Grid.kyy.*u_fft));

    du = ux.^2 + uy.^2;
end


E = sum(0.5*Para.epsilon^4 * (real(ifftn(-Grid.k.^2 .* u_fft))).^2 - ...
    (Para.tau^2/3 + Para.eta1/2 + 1)*Para.epsilon^2*du, 'all')*prod(Grid.d)...
    + w^2 - Para.B;

% E = E / prod(Grid.L);

end