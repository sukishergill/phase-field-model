%clearvars;
%clf;
m = 0.05;
eta = .01;

parfor i = 1:6
    TEMP = StartSols(i,m,eta);
    Profiles(i) = TEMP;
    x = TEMP.Domain.x;
    y = TEMP.Domain.y;
    u = TEMP.Run.u{end};
    t = TEMP.Run.t;
    E = TEMP.Run.Eu - TEMP.Run.Em;
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