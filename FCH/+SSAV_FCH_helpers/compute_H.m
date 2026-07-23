function H = compute_H(u, w, Para, F1, F2, Grid)

v = fftn(u);
v2 = fftn(u.^2);

if size(u, 2) == 1

    ux = real(ifft(-1i*Grid.k.*v));
    
    ux2 = real(ifft(-1i*Grid.k.*v2));   

    u_grad = ux.^2;

     H = (6*Para.epsilon^2 * (u.* u_grad -  ux.*ux2 - ...
        u.^2 .* real(ifftn(-Grid.k.^2 .* v))) + F1 .* F2 - ...
        Para.eta2 * F1) / w;
else
    
    ux = real(ifft2(-1i*Grid.kxx.*v));         
    uy = real(ifft2(-1i*Grid.kyy.*v));
    
    ux2 = real(ifft2(-1i*Grid.kxx.*v2));         
    uy2 = real(ifft2(-1i*Grid.kyy.*v2));
    
    u_grad1 = ux.^2 + uy.^2;
    
    H = (6*Para.epsilon^2 * (u.* u_grad1 -  ux.*ux2 - uy.*uy2 - ...
        u.^2 .* real(ifft2(-Grid.k.^2 .* v))) + F1 .* F2 - ...
        Para.eta2 * F1) / w;

end

end