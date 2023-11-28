function cfdCreatepCorr

if ~cfdIsFieldAvailable('pCorr')
    global Region
    
    Region.fluid.pCorr = Region.fluid.p;
    Region.fluid.pCorr.name = 'pCorr';
    Region.fluid.pCorr.phi(:) = 0;
end

