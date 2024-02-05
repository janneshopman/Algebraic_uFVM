function cfdSetPoisson

tol = 1E-9;     % Cleaning tolerance

global Region;

stencil = Region.foamDictionary.fvSolution.AlguFVM.LapStencil;

M = Region.operators.M;
G = Region.operators.G;
Af = Region.operators.Af;
Dnf = Region.operators.Dnf;
GamCS = Region.operators.GamCS;
GamSC = Region.operators.GamSC;

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;
ncTot = theNumberOfElements + cfdGetNumberOfBFaces;

% - Laplacian coefficient matrix for pressure
% First set up the matrix only considering internal connections
MIn = M(1:theNumberOfElements, 1:theNumberOfInteriorFaces);
GIn = G(1:theNumberOfInteriorFaces, 1:theNumberOfElements);

internalElements = 1:theNumberOfElements;

GamCSIn = GamCS(1:theNumberOfInteriorFaces, [internalElements, double(ncTot)+internalElements, double(2*ncTot)+internalElements]);
GamSCIn = GamSC([internalElements, double(ncTot)+internalElements, double(2*ncTot)+internalElements], 1:theNumberOfInteriorFaces);

if strcmp(stencil, 'compact')
    Lap = MIn*GIn;
elseif strcmp(stencil, 'wide')
    Lap = MIn*GamCSIn*GamSCIn*GIn;
else
    error('set LapStencil in fvSolutions to  \"compact\" or \"wide\".')
end

% Trick to make symmetric
minDiag = abs(min(diag(Lap)));
Lap = Lap.*(abs(Lap)>tol*minDiag);
Lap = (Lap + Lap.')/2.0;

addSource = sparse(theNumberOfElements, 1);
pRefRequired = true;

for iBPatch = 1:theNumberOfBoundaryPatches
    BType = Region.fluid.p.boundaryPatchRef{iBPatch}.type;
       
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);

    iBFaces = cfdGetBFaceIndicesForBoundaryPatch(iBPatch);
    owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);    

    if strcmp(BType, 'fixedValue') || strcmp(BType, 'noSlip')
        if strcmp(stencil, 'wide')
            error('BC not implemented for wide stencil Lap')
        end

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
        if strcmp(stencil, 'wide')
            error('BC not implemented for wide stencil Lap')
        end

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
    iCell = Region.foamDictionary.fvSolution.AlguFVM.pRefCell;
    value = Region.foamDictionary.fvSolution.AlguFVM.pRefValue;

   addSource(iCell) = addSource(iCell) + Lap(iCell,iCell)*value; 
   Lap(iCell,iCell) = 2*Lap(iCell,iCell);
end

Region.operators.Pois.Lap = Lap;
Region.operators.Pois.addSource =  addSource;