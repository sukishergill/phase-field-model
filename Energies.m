function [E1, Es, Ed] = Energies(uu, epsilon,eta,B,S,kx, ky, dA, A)

E1 = zeros(length(uu),1);
Ed = 0;
k2 = kx.^2 + ky.^2;
v = fft2(uu{1});

F = (uu{1}.^2-1).^2/4;
f = uu{1}.^3-uu{1};
ux = real(ifft2(1i*kx.*v));
uy = real(ifft2(1i*ky.*v));
lapu = real(ifft2(-k2.*v));
eps = epsilon(1);
E1(1) = compute_E;
Es = E1;

G = compute_G(uu{1});
w_old = sqrt(G+B);
    ux_old = ux;
    uy_old = uy;
    lapu_old = lapu;
for i=2:length(E1)

    u = uu{i};
    v = fft2(u);

    F = (u.^2-1).^2/4;
    f = u.^3-u;
    ux = real(ifft2(1i*kx.*v));
    uy = real(ifft2(1i*ky.*v));
    lapu = real(ifft2(-k2.*v));
    eps = epsilon(i);

    E1(i) = compute_E;
    G = compute_G(u);
    w = sqrt(G+B);
    Es(i) = compute_E2;
    Ed(i) = compute_Ediff;
    w_old = w;
    ux_old = ux;
    uy_old = uy;
    lapu_old = lapu;
end

    function E = compute_E

        Et = .5*(-eps^2*lapu+f).^2 - ... 
            eta*(eps^2*(ux.^2 + uy.^2)/2 + F);
        E = sum(Et*dA,'all')/A;
    end

    function G_loc = compute_G(U)

        G_loc = sum(3*eps^2 * U.^2 .* (ux.^2 + uy.^2) + ...
            0.5 * f.^2 - eta * F,'all') *  dA;

    end

    function E2 = compute_E2

        E2 = ...
            eps^4/2*(sum(lapu.^2,'all') + sum((2*lapu-lapu_old).^2,'all'))/2*dA + ...
            (S+eps^2*(2+eta))/2*sum((ux-ux_old).^2+(uy-uy_old).^2,'all')*dA - ...
            eps^2*(2+eta)/2*(sum(ux.^2+uy.^2,'all') + sum((2*ux-ux_old).^2+(2*uy-uy_old).^2,'all'))*dA/2 + ...
            (w^2+(2*w-w_old)^2)/2 - B;
        E2 = E2/A;

    end

    function Ed = compute_Ediff

        Ed = ...
            eps^4/2*(-sum(lapu.^2,'all') + sum((2*lapu-lapu_old).^2,'all'))/2*dA + ...
            (S+eps^2*(2+eta))/2*sum((ux-ux_old).^2+(uy-uy_old).^2,'all')*dA - ...
            eps^2*(2+eta)/2*(-sum(ux.^2+uy.^2,'all') + sum((2*ux-ux_old).^2+(2*uy-uy_old).^2,'all'))*dA/2 + ...
            (-w^2+(2*w-w_old)^2)/2;
        Ed = abs(Ed)/A;

    end

end