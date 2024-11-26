% clearvars;
% 
% LV = 151;
% k = ((2*LV)-1):-1:0;
% mvec = (1+cos(pi*(1+2*k)/4/LV));
% mvec = mvec(1:LV);
% mvec = mvec/mvec(end);
% mvec = 0*mvec;
% % etavec = linspace(.05,.05,LV);
% etavec = (1+cos(pi*(1+2*k)/4/LV));
% etavec = etavec(1:LV);
% etavec = etavec/etavec(end);
% clear k;
% tf = 200;
% 
% MM = mvec(1);
% ETA = etavec(1);
% Em = .5*(MM-MM^3)^2-ETA*(1-MM^2)^2/4;
% Em_vec = .5*(mvec-mvec.^3).^2-ETA*(1-mvec.^2).^2/4;

% run_Stripe(1) = StartSols(5,mvec(1), etavec(1));
% Evec_S(1) = run_Stripe(1).Run.Eu(end);
% epsvec_S(1) = run_Stripe(1).Run.epsilon(end);
% run_Cyl(1) = StartSols(6,mvec(1), etavec(1));
% Evec_C(1) = run_Cyl(1).Run.Eu(end);
% epsvec_C(1) = run_Cyl(1).Run.epsilon(end);
% run_Square(1) = StartSols(3,mvec(1), etavec(1));
% Evec_Sq(1) = run_Square(1).Run.Eu(end);
% epsvec_Sq(1) = run_Square(1).Run.epsilon(end);

% run_Rand2(1) = StartSols(1,mvec(1), etavec(1),.25);
% Evec_R2(1) = run_Rand2(1).Run.Eu(end);
% epsvec_R2(1) = run_Rand2(1).Run.epsilon(end);
% run_Rand1(1) = StartSols(2,mvec(1), etavec(1),.25);
% Evec_R1(1) = run_Rand1(1).Run.Eu(end);
% epsvec_R1(1) = run_Rand1(1).Run.epsilon(end);

% OldData(1) = run_Rand2;
% OldData(2) = run_Rand1;
% OldData(3) = run_Stripe;
% OldData(4) = run_Cyl;
% OldData(5) = run_Square;

% OldData(1) = run_Stripe;
% OldData(2) = run_Cyl;
% OldData(3) = run_Square;
% 
% %fname = ['Run_m_',num2str(ETA),'.mat'];
fname = 'Run_m_0.mat';

for i=2:LV
    parfor k=3:5
        [Domain, Run, Params] = TEMPFUNC(i,k-2,OldData, tf,mvec,etavec);
        TDomain(k-2) = Domain;
        TRun(k-2) = Run;
        TParams(k-2) = Params;
    end
    MM = mvec(i);
    ETA = etavec(i);
    Em = .5*(MM-MM^3)^2-ETA*(1-MM^2)^2/4;
    for l=3:5
        tt.Domain = TDomain(l-2);
        tt.Run = TRun(l-2);
        tt.Params = TParams(l-2);
        if l == 1
            run_Rand2(i) = tt;
            Evec_R2(i) = run_Rand2(i).Run.Eu(end);
            epsvec_R2(i) = run_Rand2(i).Run.epsilon(end);
            OldData(1) = tt;
        elseif l == 2
            run_Rand1(i) = tt;
            Evec_R1(i) = run_Rand1(i).Run.Eu(end);
            epsvec_R1(i) = run_Rand1(i).Run.epsilon(end);
            OldData(2) = tt;
        elseif l == 3
            run_Stripe(i) = tt;
            Evec_S(i) = run_Stripe(i).Run.Eu(end);
            epsvec_S(i) = run_Stripe(i).Run.epsilon(end);
            OldData(3) = tt;
        elseif l == 4
            run_Cyl(i) = tt;
            Evec_C(i) = run_Cyl(i).Run.Eu(end);
            epsvec_C(i) = run_Cyl(i).Run.epsilon(end);
            OldData(4) = tt;
        elseif l == 5
            run_Square(i) = tt;
            Evec_Sq(i) = run_Square(i).Run.Eu(end);
            epsvec_Sq(i) = run_Square(i).Run.epsilon(end);
            OldData(5) = tt;
        end
    end
    figure(1);
    % plot(mvec(1:i),Evec_S(1:i),mvec(1:i),Evec_C(1:i), ...
    %     mvec(1:i),Evec_R2(1:i),mvec(1:i),Evec_R1(1:i),mvec(1:i),Evec_Sq(1:i));
    % legend('Stripe','Cylinder','Rand 2d','Rand 1d', 'Square');
    plot(etavec(1:i),Evec_S(1:i),etavec(1:i),Evec_C(1:i), ...
        etavec(1:i),Evec_Sq(1:i));
    legend('Stripe','Cylinder', 'Square');
    ym = 1.2*min(0,min(Evec_S(1:i)));
    yM = 2*max(0,max(Evec_S(1:i)));
    ylim([ym yM])
    drawnow;
    figure(2);
    contourf(run_Cyl(i).Run.u{end});colorbar;axis equal;drawnow;
    figure(3);
    plot(etavec(1:i),epsvec_S(1:i),etavec(1:i),epsvec_C(1:i), ...
        etavec(1:i),epsvec_Sq(1:i));
    legend('Stripe','Cylinder', 'Square');
    drawnow;
    figure(1);
    save(fname);
end

function [Domain, Run, Params] = TEMPFUNC(I,K, OldData,tf,mvec,etavec)

run = OldData(K);
run.Params.m = mvec(I);
run.Params.eta = etavec(I);
% if K < 3
%     OutRun = StartSols(K,mvec(I), etavec(I), eps0);
% else
    OutRun = ContSolsRun(run, tf);
%end

Domain = OutRun.Domain;
Run = OutRun.Run;
Params = OutRun.Params;

end