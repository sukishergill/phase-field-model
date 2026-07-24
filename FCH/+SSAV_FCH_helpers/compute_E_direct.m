function E = compute_E_direct(u, uhat, eps, eta, Grid)
% Direct evaluation of the physical FCH energy functional straight from u
% (no SAV auxiliary variable). Used by arclength_FCH.m, and by the SAV
% time-stepper (alongside its own SAV-tracked Eu) so the two can be
% compared to see how much the SAV-tracked energy drifts from the true
% energy over the course of a run.

[F, dF, ~] = SSAV_FCH_helpers.compute_F(u, 0);

ux = real(ifftn(-1i*Grid.k.*uhat));
uxx = real(ifftn(-Grid.k.^2 .* uhat));

E = sum(0.5*(-eps^2*uxx + dF).^2 - eta*(eps^2*0.5*ux.^2 + F))*prod(Grid.d);

E = E / prod(Grid.L);

end
