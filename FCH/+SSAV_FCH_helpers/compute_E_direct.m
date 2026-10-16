function E = compute_E_direct(u, uhat, eps, eta, Grid)
% Direct evaluation of the physical FCH energy functional straight from u
% (no SAV auxiliary variable). Used by arclength_FCH.m, and by the SAV
% time-stepper (alongside its own SAV-tracked Eu) so the two can be
% compared to see how much the SAV-tracked energy drifts from the true
% energy over the course of a run.

[F, dF, ~] = SSAV_FCH_helpers.compute_F(u, 0);

% |grad u|^2. First derivatives need the signed wavenumber Grid.kxx, not
% Grid.k = |k|: -1i*|k|.*uhat is anti-Hermitian, so its ifft is purely
% imaginary and real() of it is identically zero.
if size(uhat, 2) == 1
    ux = real(ifft(-1i*Grid.kxx.*uhat));
    du = ux.^2;
else
    ux = real(ifft2(-1i*Grid.kxx.*uhat));
    uy = real(ifft2(-1i*Grid.kyy.*uhat));
    du = ux.^2 + uy.^2;
end

% Laplacian (Grid.k.^2 = kx^2 + ky^2 in 2D)
lap_u = real(ifftn(-Grid.k.^2 .* uhat));

E = sum(0.5*(-eps^2*lap_u + dF).^2 - eta*(eps^2*0.5*du + F), 'all')*prod(Grid.d);

E = E / prod(Grid.L);

end
