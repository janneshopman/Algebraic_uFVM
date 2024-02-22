function cfdCreateUcp

if ~cfdIsFieldAvailable('Ucp')
    global Region
    
    Region.fluid.Ucp = Region.fluid.U;
    Region.fluid.Ucp.name = 'Ucp';    
    Region.fluid.Ucp.phi(:) = 0;
end

