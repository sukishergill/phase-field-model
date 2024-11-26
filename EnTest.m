%runS1a = StartSols(4,.1, .01);
%runSa = StartSols(5,.1, .01);
runCa = StartSols(6,.1, .01);

eta = .01;
[nx, ny] = size(runC.kx);
dx = runC.Lx/nx;
dy = runC.Ly/ny;
dA = dx*dy;
A = runC.Lx*runC.Ly;

[EC, ECs, Ediff] = Energies(runCa.uAll,runCa.epsilon,eta,runCa.B,runCa.S, ... 
    runCa.kx, runCa.ky, dA, A);