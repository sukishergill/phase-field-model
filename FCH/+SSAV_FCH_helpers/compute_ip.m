function innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dxdy)

innprod_Hu = sum(sum(H_snew .* psi_r))*dxdy / ...
    (1 - 0.5*sum(sum(H_snew .* psi_H)) * dxdy);

end