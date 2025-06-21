function innprod_Hu = compute_ip(H_snew, psi_r, psi_H, Grid)
% This fuction computes the inner product of H^{*,n+1} and u^{n+1}

innprod_Hu = sum(H_snew .* psi_r, 'all') * prod(Grid.d)...
    / (1 - 0.5*sum(H_snew .* psi_H, 'all') * prod(Grid.d));


end