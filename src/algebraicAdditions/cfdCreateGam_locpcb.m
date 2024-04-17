function cfdCreateGam_locpcb

if ~cfdIsFieldAvailable('Gam')
    global Region
    
    Region.fluid.Gam = Region.fluid.p;
    Region.fluid.Gam.name = 'Gam';    
    Region.fluid.Gam.dimensions = [0,0,0,0,0,0,0];
    Region.fluid.Gam.phi(:) = 0 * cfdGetField('p');
    cfdUpdateVectorFieldForAllBoundaryPatches('Gam');
end
