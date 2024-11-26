rerun = [1,5,6];

parfor i = rerun
    TEMP = ContSolsRun(Profiles(i),eta,tf);
    Profiles(i) = TEMP;
    x = TEMP.x;
    y = TEMP.y;
    u = TEMP.u;
    t = TEMP.t;
    E = TEMP.E - TEMP.Em;
    figure(1);
    subplot(2,3,i)
    contourf(x,y,u);
    axis equal;
    axis off;
    clim([-1 1]);
    colorbar;
    figure(2);
    subplot(2,3,i)
    plot(t(2:end),E(2:end));
    drawnow;
end

fname = ['StartSols_m_',num2str(m),'_eta',num2str(eta),'.mat'];
save(fname,'Profiles');

%-------------------------------------------------
function run = ContSolsRun(run, eta, tf)

u = run.u;
[nx, ny] = size(u);
Lx = run.Lx;
Ly = run.Ly;
eps0 = run.epsilon;
dt = 0.0005;
dk = 10;
B = (Lx*Ly)^2;
S = 2;

[ xx, yy, kxx, kyy, tt, uu, Eu, Em, epsilon] = ...
    FCH_BDF2_SAV_JF_kappa ( nx, ny, Lx, Ly, dt, tf, eps0, eta, B, S, u, dk);

run = struct;
run.u0 = uu{1};
run.u = uu{end};
run.epsilon = epsilon;
run.E = Eu;
run.t = tt;
run.Em = Em;
run.x = xx;
run.y = yy;
run.kx = kxx;
run.ky = kyy;
run.Lx = Lx;
run.Ly = Ly;
run.B = B;
run.S = S;
run.dt = dt;
run.eta = eta;
run.m = mean(run.u);

end