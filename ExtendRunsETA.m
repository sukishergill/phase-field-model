clearvars;

LV = 151;
k = ((2*LV)-1):-1:0;
MM = 0.1;

etavec = (1+cos(pi*(1+2*k)/4/LV));
etavec = etavec(1:LV);
etavec = etavec/etavec(end);
clear k;

ETA = etavec(1);
Em = .5*(MM-MM^3)^2-ETA*(1-MM^2)^2/4;
Em_vec = .5*(MM-MM.^3).^2-etavec*(1-MM.^2).^2/4;

run_Stripe(1) = StartSols(5,MM, etavec(1));
Evec_S(1) = run_Stripe(1).Run.Eu(end);
epsvec_S(1) = run_Stripe(1).Run.epsilon(end);
run_Cyl(1) = StartSols(6,MM, etavec(1));
Evec_C(1) = run_Cyl(1).Run.Eu(end);
epsvec_C(1) = run_Cyl(1).Run.epsilon(end);
run_Square(1) = StartSols(3,MM, etavec(1));
Evec_Sq(1) = run_Square(1).Run.Eu(end);
epsvec_Sq(1) = run_Square(1).Run.epsilon(end);

OldData(1) = run_Stripe;
OldData(2) = run_Cyl;
OldData(3) = run_Square;

fname = 'Run_m_0.mat';
tf = 200;

for i=2:LV
    for k=1:3
        [Domain, Run, Params] = TEMPFUNC(i,k,OldData, tf,MM,etavec);
        TDomain(k) = Domain;
        TRun(k) = Run;
        TParams(k) = Params;
    end
    ETA = etavec(i);
    Em = .5*(MM-MM^3)^2-ETA*(1-MM^2)^2/4;
    for l=1:3
        tt.Domain = TDomain(l);
        tt.Run = TRun(l);
        tt.Params = TParams(l);
        if l == 1
            run_Stripe(i) = tt;
            Evec_S(i) = run_Stripe(i).Run.Eu(end);
            epsvec_S(i) = run_Stripe(i).Run.epsilon(end);
            OldData(3) = tt;
        elseif l == 2
            run_Cyl(i) = tt;
            Evec_C(i) = run_Cyl(i).Run.Eu(end);
            epsvec_C(i) = run_Cyl(i).Run.epsilon(end);
            OldData(4) = tt;
        elseif l == 3
            run_Square(i) = tt;
            Evec_Sq(i) = run_Square(i).Run.Eu(end);
            epsvec_Sq(i) = run_Square(i).Run.epsilon(end);
            OldData(5) = tt;
        end
    end
    figure(1);
    plot(etavec(1:i),Evec_S(1:i),etavec(1:i),Evec_C(1:i), ...
        etavec(1:i),Evec_Sq(1:i));
    legend('Stripe','Cylinder', 'Square');
    ym = 1.2*min(0,min(Evec_S(1:i)));
    yM = 2*max(0,max(Evec_S(1:i)));
    ylim([ym yM])
    drawnow;
    figure(2);
    subplot(1,3,1)
    contourf(run_Cyl(i).Run.u{end});colorbar;axis equal; axis off;
        subplot(1,3,2)
    contourf(run_Square(i).Run.u{end});colorbar;axis equal; axis off;
        subplot(1,3,3)
    contourf(run_Stripe(i).Run.u{end});colorbar;axis equal; axis off;
    drawnow;
    figure(3);
    plot(etavec(1:i),epsvec_S(1:i),etavec(1:i),epsvec_C(1:i), ...
        etavec(1:i),epsvec_Sq(1:i));
    legend('Stripe','Cylinder', 'Square');
    drawnow;
    figure(1);
    save(fname);
end

function [Domain, Run, Params] = TEMPFUNC(I,K, OldData,tf,MM,etavec)

run = OldData(K);
run.Params.m = MM;
for eta_count = 1:10
    CurrentEta = etavec(I-1) + eta_count/10*(etavec(I) - etavec(I-1));
    run.Params.eta = CurrentEta;
    run = ContSolsRun(run, 5);
end
run.Params.eta = etavec(I);
OutRun = ContSolsRun(run, tf);

Domain = OutRun.Domain;
Run = OutRun.Run;
Params = OutRun.Params;

end