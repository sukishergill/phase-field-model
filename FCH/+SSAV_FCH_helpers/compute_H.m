function H = compute_H(u, w, Para, F1, F2, Grid)

v = fft2(u);
v2 = fft2(u.^2);

ux = real(ifft2(-1i*Grid.kxx.*v));         
uy = real(ifft2(-1i*Grid.kyy.*v));

ux2 = real(ifft2(-1i*Grid.kxx.*v2));         
uy2 = real(ifft2(-1i*Grid.kyy.*v2));

u_grad1 = ux.^2 + uy.^2;

H = (6*Para.epsilon^2 * (u.* u_grad1 -  ux.*ux2 - uy.*uy2 - ...
    u.^2 .* real(ifft2(-Grid.k.^2 .* v))) + F1 .* F2 - ...
    Para.eta2 * F1) / w;

end