function E_n = compute_En(u_fft, w, Para, Grid, D, model)

% num_fft_out = num_fft_in;
e = Para.epsilon^2/2;

% um = fft2(u - Para.m);      num_fft_out = num_fft_out + 1;
um_fft = u_fft;         um_fft(1,1) = 0;

v = -Grid.inv_k .* um_fft;
v(Grid.k == 0) = 0;

vx = real(ifft2(-1i*Grid.kxx.*v));         
vy = real(ifft2(-1i*Grid.kyy.*v));

dv = vx.^2 + vy.^2;

% u = fft2(u);        num_fft_out = num_fft_out + 1;

if model == 2

    du = real(ifft2(D.*u_fft));
else

    ux = real(ifft2(-1i*Grid.kxx.*u_fft));         
    uy = real(ifft2(-1i*Grid.kyy.*u_fft));

    du = ux.^2 + uy.^2;
end

E_n = sum(e *du + Para.alpha*e * dv, 'all')*prod(Grid.d) +  w^2 - Para.B;

end