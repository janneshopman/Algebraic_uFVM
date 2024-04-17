function cfdCreateUMag_locpcb

if ~cfdIsFieldAvailable('UMag')
    global Region
    
    Region.fluid.UMag = Region.fluid.p;
    Region.fluid.UMag.name = 'UMag';    
    Region.fluid.UMag.dimensions = [0,2,-2,0,0,0,0];
    Region.fluid.UMag.phi(:) = 0 * cfdGetField('p');
    cfdUpdateVectorFieldForAllBoundaryPatches('UMag');
end
