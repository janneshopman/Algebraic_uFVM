function cfdSetEigenvalues

global Region;

fprintf('\nSetting Stability Region\n');

nElem = Region.mesh.numberOfElements;
nBElem = Region.mesh.numberOfBElements;
nContr = nElem+nBElem;
innerElem = [1:nElem,nContr+1:nElem+nContr];

Region.mesh.innerElemVec = innerElem;

aTab = Region.foamDictionary.fvSchemes.ddtSchemes.RK.aTab;
nStages = Region.foamDictionary.fvSchemes.ddtSchemes.RK.nStages;
A = aTab(1:nStages,:);
b = aTab(end,:)';
Region.evals.R = @(z) 1 + z*b'*((eye(nStages)-z*A)\ones(nStages,1));
Region.evals.max_ev = 0;
Region.evals.maxhf = 0;

figure( 1 );
if Region.foamDictionary.fvSolution.AlguFVM.showRegion
    rr = linspace(0,3,300);
    phi = linspace(0,pi,300);
    for i=1:length(rr)
        for j=1:length(phi)
            zz = rr(i)*cos(phi(j))+1i*rr(i)*sin(phi(j));
            Rz = Region.evals.R(zz);
            if abs(abs(Rz)-1.0) < 1e-3
                scatter(real(zz),imag(zz),'k')
                hold on
                xlabel('Real','interpreter','Latex');
                ylabel('Imag','interpreter','Latex');
                axis equal
                grid on
            end
        end
    end
end


