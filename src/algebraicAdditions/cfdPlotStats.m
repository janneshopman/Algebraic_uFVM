function cfdPlotStats

global Region;

if Region.foamDirectory.fvSolution.AlguFVM.plotSATStats

    dts = Region.stats.timestep;
    phis = Region.stats.angle;
    maxEV = Region.stats.maxEV;
        
    ite = 1:length(dts);
    
    figure( 2 );
    filename = strcat(pwd,'/time_evolution.pdf');
    
    plot(ite,dts,'-k')
    hold on
    grid on
    xlabel('Iteration','Interpreter','latex')
    ylabel('$\Delta t$','Interpreter','latex')
    ax = gca;
    exportgraphics(ax,filename);
    
    figure( 3 );
    filename = strcat(pwd,'/eigenvalue_evolution.pdf');
    yyaxis left
    plot(ite,maxEV,'-k')
    ylabel('$\lambda\tilde{\Delta t}$','Interpreter','latex')
    yyaxis right
    plot(ite,phis/(pi/2),'--k')
    hold on
    grid on
    xlabel('Iteration','Interpreter','latex')
    ylabel('$\frac{\varphi}{\pi/2}$','Interpreter','latex')
    ax = gca;
    exportgraphics(ax,filename);

end



