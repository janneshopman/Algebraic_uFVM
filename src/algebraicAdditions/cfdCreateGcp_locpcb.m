function cfdCreateGcp_locpcb

if ~cfdIsFieldAvailable('Gcp')
    global Region
    
    Region.fluid.Gcp = Region.fluid.U;
    Region.fluid.Gcp.name = 'Gcp';    
    Region.fluid.Gcp.dimensions = [0,2,-2,0,0,0,0];
    Region.fluid.Gcp.phi(:) = Region.operators.Gc * cfdGetField('p');
    cfdUpdateVectorFieldForAllBoundaryPatches('Gcp');
end