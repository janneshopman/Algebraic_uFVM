function cfdCreateUf

if ~cfdIsFieldAvailable('Uf')
    global Region

    % Setup mesh field
    cfdSetupMeshField('Uf', 'surfaceScalarField');
    
    % Store initial value in mesh field
    theMeshField = cfdGetMeshField('Uf');
    
    % Read dimensions
    theMeshField.dimensions = [0, 1, -1, 0, 0, 0, 0];
    
    % Read and store internal field
    if cfdIsFieldAvailable('phi')
        phi = cfdGetField('phi');
        theMeshField.phi = Region.operators.Af\phi;
    else    
        U = cfdGetField('U');  
        theMeshField.phi = Region.operators.GamCSM*U;
    end

    % Store mesh field in data base
    cfdSetMeshField(theMeshField);
end

