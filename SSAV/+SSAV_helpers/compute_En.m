function E_n = compute_En(u_fft, w, Para, Grid, D, model, dim)

% num_fft_out = num_fft_in;
e = Para.epsilon^2/2;

% um = fftn(u - Para.m);      %num_fft_out = num_fft_out + 1;
% um_fft = u_fft;         um_fft(1,1) = 0;
um_fft = u_fft - Para.m;

v = -Grid.inv_k .* um_fft;
v(Grid.k == 0) = 0;

vx = real(ifftn(-1i*Grid.kxx.*v));         
vy = real(ifftn(-1i*Grid.kyy.*v));

dv = vx.^2 + vy.^2;

% u = fftn(u);        num_fft_out = num_fft_out + 1;

if model == 2

    du = real(ifftn(D.*u_fft));
else

    ux = real(ifftn(-1i*Grid.kxx.*u_fft));         
    uy = real(ifftn(-1i*Grid.kyy.*u_fft));

    du = ux.^2 + uy.^2;

    if dim == 3
        uz = real(ifftn(-1i*Grid.kzz.*u_fft));
        du = du + uz.^2;
    end

end

E_n = sum(e *du + Para.alpha*e * dv, 'all')*prod(Grid.d) +  w^2 - Para.B;

end