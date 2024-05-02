function cfdSelfAdaptiveTimestep

global Region;

if Region.foamDictionary.fvSolution.AlguFVM.enableSAT

    nite = Region.time.nTimeStep;
    
    fs = Region.foamDictionary.fvSolution.AlguFVM.sampleFreq;
    
    if rem(nite,fs) == 0
        [max_ev,idx] = max(abs(Region.evals.ev));
        ev = Region.evals.ev(idx);
        Region.evals.max_ev = ev;
        phi = atan(imag(ev)/abs(real(ev)));
    
        rr = linspace(0.1,4,10000);
        for i=1:length(rr)
           zz = -rr(i)*cos(phi)+1i*rr(i)*sin(phi);
           Rz = Region.evals.R(zz);
           if abs(abs(Rz)-1.0) < 1e-3
                dT = rr(i)/max_ev;
                ii=i;
           end
        end
        Region.evals.maxhf = rr(ii);
        fdt = Region.foamDictionary.fvSolution.AlguFVM.fdT;
    
        Region.foamDictionary.controlDict.deltaT = fdt*dT;
    end
    
    Region.stats.maxEV(nite) = Region.evals.maxhf;
    Region.stats.timestep(nite) = Region.foamDictionary.controlDict.deltaT;
    Region.stats.angle(nite) = atan(imag(Region.evals.max_ev)...
        /abs(real(Region.evals.max_ev)));
    Region.stats.evmax(nite) = abs(Region.evals.max_ev);

end



