function cfdSetPoisson

global Region;

M = Region.operators.M;
G = Region.operators.G;
Af = Region.operators.Af;
Dnf = Region.operators.Dnf;

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;

% - Laplacian coefficient matrix for pressure
% First set up the matrix only considering internal connections
Lap = M(1:theNumberOfElements, 1:theNumberOfInteriorFaces)*G(1:theNumberOfInteriorFaces, 1:theNumberOfElements);
addSource = sparse(theNumberOfElements, 1);
pRefRequired = true;

for iBPatch = 1:theNumberOfBoundaryPatches
    BType = Region.fluid.p.boundaryPatchRef{iBPatch}.type;
       
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);

    iBFaces = cfdGetBFaceIndicesForBoundaryPatch(iBPatch);
    owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);    

    if strcmp(BType, 'fixedValue') || strcmp(BType, 'noSlip')
        % Add diagonal coefficient to Lap
        % Add RHS modification
        % Set pRefRequired

        pRefRequired = false;

        % Initialize BC value as 0 for noSlip
        theBCValue = zeros(length(iBFaces), 1);
        
        % Change value if fixedValue
        if strcmp(BType, 'fixedValue')
            theBCValue = cfdValueForBoundaryPatch(Region.fluid.p.name, iBPatch);

            % Extend to vector if uniform value is provided
            if size(theBCValue, 1)==1
                theBCValue = theBCValue*ones(length(iBFaces),1);
            end
        end

        % Calculate coefficients
        faceCoefs = -2*(diag(Dnf).\diag(Af));

        % Modify diagonal and RHS
        Lap = Lap + sparse(owners_b, owners_b, faceCoefs(iBFaces), theNumberOfElements, theNumberOfElements);
        addSource(owners_b) = addSource(owners_b) + faceCoefs(iBFaces).*theBCValue;
    elseif strcmp(BType, 'cyclic')
        % Add internal connections
        % RHS b adds 0
        iNBPatch = theBCInfo.neighbourPatchId;
                
        owners_Nb = cfdGetOwnersSubArrayForBoundaryPatch(iNBPatch);  

        indicesFaceD = sub2ind(size(Af), iBFaces, iBFaces);
        indicesCellD = sub2ind(size(Lap), owners_b, owners_b);
        indicesCellOD = sub2ind(size(Lap), owners_b, owners_Nb);

        coeffs = Af(indicesFaceD)./Dnf(indicesFaceD);

        Lap(indicesCellOD) = Lap(indicesCellOD) + coeffs;
        Lap(indicesCellD) = Lap(indicesCellD) - coeffs;
    elseif strcmp(BType, 'zeroGradient') ...
            || strcmp(BType, 'slip') ...
            || strcmp(BType, 'outlet') ...
            || strcmp(BType, 'symmetry') ...
            || strcmp(BType, 'empty')
        % Do nothing, Lap(i,i) only contains internal coefficients
        % RHS b adds 0
    end
end

if pRefRequired
    iCell = Region.foamDictionary.fvSolution.SIMPLE.pRefCell;
    value = Region.foamDictionary.fvSolution.SIMPLE.pRefValue;

    addSource(iCell) = addSource(iCell) + Lap(iCell,iCell)*value; 
    Lap(iCell,iCell) = 2*Lap(iCell,iCell);
end

Region.operators.Pois.Lap = Lap;
Region.operators.Pois.addSource =  addSource;