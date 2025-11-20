function G = compute_G(u, Para, F, F1, Grid)

v = fftn(u);

ux = real(ifft2(-1i*Grid.kxx.*v));         
uy = real(ifft2(-1i*Grid.kyy.*v));

del_u = real(ifft2(-Grid.k2 .* u));

% if Para.W == 1
% 
%     G = sum(-Para.epsilon^2*del_u .* (Wq1.^2 + Para.zeta*u) + ...
%     0.5 * Wq1.^2 - Para.eta2 * Wq, 'all') * prod(Grid.d);
% 
% else
% 
%     G = sum(3*Para.epsilon^2 * u.^2 .* (ux.^2 + uy.^2) + ...
%     0.5 * Wq1.^2 - Para.eta2 * Wq, 'all') * prod(Grid.d);
% 
% end

G = sum(3*Para.epsilon^2 * u.^2 .* (ux.^2 + uy.^2) + ...
    0.5 * F1.^2 - Para.eta2 * F, 'all') * prod(Grid.d);
end