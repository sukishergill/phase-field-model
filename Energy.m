function E = Energy(eps,eta,u,kx,ky)

E = zeros(size(eps));

k2 = kx.^2 + ky.^2;

F = .25*(1-u.^2).^2;
Fp = u.^3-u;

v = fft2(u);
ux = real(ifft2(1i*kx.*v));
uy = real(ifft2(1i*ky.*v));
gradu2 = ux.^2+uy.^2;
lapu = real(ifft2(-k2.*v));

for i=1:length(eps)
     Et = .5*(-eps(i)^2*lapu+Fp).^2 - eta*(eps(i)^2*gradu2/2 + F);
    E(i) = mean(Et,'all');
end
