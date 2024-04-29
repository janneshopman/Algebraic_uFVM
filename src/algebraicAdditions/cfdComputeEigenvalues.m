function cfdComputeEigenvalues

global Region

nite = Region.time.nTimeStep;
fs = Region.foamDictionary.fvSolution.AlguFVM.sampleFreq;
lim = Region.foamDictionary.fvSolution.AlguFVM.evalPercent;

if rem(nite,fs) == 0

    F = full(Region.evals.F);
    innerElem = Region.mesh.innerElemVec;
    F = F(innerElem,innerElem);
    deltaT = cfdGetDeltaT;
    [evecs,evals] = eig(F);
    Region.evals.ev = diag(evals);

    U = cfdGetField('U');
    Uproj = abs(evecs*U(Region.mesh.innerElemVec));
    maxProj = max(Uproj);
    Uproj = Uproj/maxProj;

    if Region.foamDictionary.fvSolution.AlguFVM.showRegion
        figure( 1 );
        rreal = real(Region.evals.ev(Uproj >= lim));
        iimag = abs(imag(Region.evals.ev(Uproj >= lim)));
        filename = strcat(pwd,'/evals.pdf');
        scatter(deltaT*rreal,deltaT*iimag,[],Uproj(Uproj >= lim),'filled');
        axis equal
        colorbar
        colormap jet
        ax = gca;
        exportgraphics(ax,filename);


    end
    Region.evals.eigenbounds.real = max(abs(real(Region.evals.ev)));
    Region.evals.eigenbounds.imag = max(imag(Region.evals.ev));

end
