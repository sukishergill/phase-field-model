function EPS = EpsStar(u,eta,kx,ky)

f = u.^3-u;
v = fft2(u);
lapu = real(ifft2(-(kx.^2+ky.^2).*v));
ux = real(ifft2(1i*kx.*v));
uy = real(ifft2(1i*ky.*v));

numer = lapu.*f + eta*(ux.^2+uy.^2)/2;
denom = lapu.^2;

if mean(numer,'all') > 0
    EPS = sqrt(mean(numer,'all')/mean(denom,'all'));
else
    EPS = -1;
end