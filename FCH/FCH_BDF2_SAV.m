function Results = FCH_BDF2_SAV(Grid, Time, Para, u)

dt = Time.dt;

% define number of time steps and time steps for plotting purposes
nmax = round(Time.tf / dt);
nplt = floor( (Time.tf / 100) / dt);
tt = zeros(1, 101);
tt(1) = 0;

uu = cell(101, 1);
uu{1} = u;

mass = zeros(1, 101);
mass(1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));

Eu = zeros(1, 101);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Wq, Wq1, Wq2] = SSAV_FCH_helpers.compute_Wq(u, Para);

G = SSAV_FCH_helpers.compute_G(u, Para, Wq, Wq1, Grid);

w_old = SSAV_FCH_helpers.compute_w(G, Para.B);

Eu(1) = SSAV_FCH_helpers.compute_E(u, w_old, Para, Grid);

H = SSAV_FCH_helpers.compute_H(u, w_old, Para, Wq1, Wq2, Grid);
H_old = H;

r = w_old - 0.5 * sum(sum(H.*u))*prod(Grid.d);

H2 = fftn(H);

v = fftn(u);

r_hat = u/dt + (Para.epsilon^2*(2 + Para.eta1 + Para.tau^2/3) + Para.S) * ...
    real(ifftn(Grid.k4 .* v)) + r * real(ifftn(-Grid.k2 .* H));

% define linear operator P for BDF1
P = 1/dt + Para.epsilon^4 * Grid.k6 + Para.S * Grid.k4;

r_hat = fftn(r_hat);
psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

psi_H = P .\ (G .* H2);     psi_H = ifftn(psi_H, 'symmetric');

innprod_Hu = SSAV_FCH_helpers.compute_ip(H, psi_r, psi_H, prod(Grid.d));

w = innprod_Hu + r;

u_old = u;

u = 0.5 * innprod_Hu * psi_H + psi_r;

% define linear operator P for BDF2
P = P + 1/(2*dt);

if nplt == 1
    
    uu{2} = u;    
    
    tt(2) = dt;  

    w_vals(2) = w;
    
    Eu(2) = SSAV_FCH_helpers.compute_E(u, w, Para, Grid);
    
    mass(2) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2 : nmax

    t = j * dt;

    u_snew = SSAV_FCH_helpers.compute_u_snew(u, u_old);

    [Wq, Wq1, Wq2] = SSAV_FCH_helpers.compute_Wq(u_snew, Para);

    G = SSAV_FCH_helpers.compute_G(u, Para, Wq, Wq1, Grid);

    w_snew = SSAV_FCH_helpers.compute_w(G, Para.B);

    H_snew = SSAV_FCH_helpers.compute_H(u_snew, w_snew, Para, Wq1, Wq2,...
        Grid);

    H = fftn(H_snew);
    u_snew = fftn(u_snew);

    r = (4*w - w_old)/3 - 0.5 * sum(sum(H_snew .* (4*u - u_old)/3)) ...
        * prod(Grid.d);

    r_hat = (4*u - u_old)/(2*dt) + (Para.epsilon^2*(2 + Para.eta1 + ...
        2*Para.tau^2/3) + Para.S) * real(ifftn(Grid.k4 .* u_snew)) + ...
        r * real(ifftn(-Grid.k2 .* H));

    r_hat = fftn(r_hat);

    psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

    psi_H = P .\ (G .* H);     psi_H = ifftn(psi_H, 'symmetric');

    innprod_Hu = SSAV_FCH_helpers.compute_ip(H_snew, psi_r, psi_H, ...
        prod(Grid.d));

    w_new = 0.5*innprod_Hu + r;

    u_new = 0.5*innprod_Hu * psi_H + psi_r;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;

    if (mod(j, nplt) == 0)                  

        Eu(j/nplt + 1) = SSAV_FCH_helpers.compute_E(u, w, Para, Grid);
        
        mass(j/nplt + 1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));
        
        uu{j/nplt + 1} = u;
        
        tt(j/nplt + 1) = t;

        w_vals(j/nplt + 1) = w;
    end

end

Eu = Eu / prod(Grid.L);

Results.tt = tt;
Results.Eu = Eu;
Results.uu = uu;
Results.mass = mass;

end
