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
    U = cfdGetField('U');  
    theMeshField.phi = Region.operators.GamCSM*U;
              
    % Store mesh field in data base
    cfdSetMeshField(theMeshField);
end

