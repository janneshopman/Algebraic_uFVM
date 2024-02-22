function cfdSetPoisson

tol = 1E-9;     % Cleaning tolerance

global Region;

LapStencil = Region.foamDictionary.fvSolution.AlguFVM.LapStencil;
interpolation = Region.foamDictionary.fvSolution.AlguFVM.interpolation;

Af = Region.operators.Af;
Dnf = Region.operators.Dnf;
Nf = Region.operators.Nf;
M = Region.operators.M;
OmegaS = Region.operators.OmegaS;
OmegaIn = Region.operators.OmegaIn;
GamCS = Region.operators.GamCS;

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBFaces = cfdGetNumberOfBFaces;
theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;

ncTot = theNumberOfElements + theNumberOfBFaces;
nfTot = theNumberOfInteriorFaces + theNumberOfBFaces;
internalElements = 1:theNumberOfElements;
internalFaces = 1:theNumberOfInteriorFaces;

% Laplacian coefficient matrix for pressure
% Matrices to set up and adjust
MAdj = M(internalElements, internalFaces);

OmegaSAdj = OmegaS(internalFaces, internalFaces);

GamCSAdj = GamCS(internalFaces, [internalElements, double(ncTot)+internalElements, double(2*ncTot)+internalElements]);

% Loop over patches and adjust matrices for cyclic boundaries
checkPatches = ones(theNumberOfBoundaryPatches, 1);

for iBPatch = 1:theNumberOfBoundaryPatches

    % Retrieve relevant owner patch information
    BType = Region.fluid.p.boundaryPatchRef{iBPatch}.type;
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);
    iBFaces = cfdGetBFaceIndicesForBoundaryPatch(iBPatch);
    owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);
    theNumberOfFacesInPatch = cfdGetNumberOfFacesForBoundaryPatch(iBPatch);

    if strcmp(BType, 'cyclic') && checkPatches(iBPatch) == 1

        % Retrieve relevant neighbour patch information
        iNBPatch = theBCInfo.neighbourPatchId;
        iNBFaces = cfdGetBFaceIndicesForBoundaryPatch(iNBPatch);             
        owners_Nb = cfdGetOwnersSubArrayForBoundaryPatch(iNBPatch);  

        % Don't check neighbour patch after this
        checkPatches(iNBPatch) = 0;

        % Calculate face volumes
        fVolOwns = diag(OmegaS(iBFaces, iBFaces));
        fVolNeis = diag(OmegaS(iNBFaces, iNBFaces));
        fVolSum = fVolOwns + fVolNeis;

        %%% Adjust GamCS   
        % Indices for x y z faces
        iBFacesxyz = reshape(double(nfTot) * [0, 1, 2] + iBFaces, [], 1);

        % Access face normals
        Nfs = Nf(iBFaces, iBFacesxyz).';
        
        % Set weights depending on interpolation method
        if strcmp(interpolation, 'linear')
            wOwnNeis = [fVolNeis./fVolSum, fVolOwns./fVolSum];
        elseif strcmp(interpolation, 'midpoint')
            wOwnNeis = ones(theNumberOfFacesInPatch, 2)/2;
        elseif strcmp(interpolation, 'volumetric')
            wOwnNeis = [fVolOwns./fVolSum, fVolNeis./fVolSum];
        end

        % Fill in row and column indices
        GamCSRows = repmat(1:theNumberOfFacesInPatch, 1, 6).';

        GamCSCols = ...
        [ ...
            reshape(double(theNumberOfElements) * [0, 1, 2] + owners_b,  [], 1); ...
            reshape(double(theNumberOfElements) * [0, 1, 2] + owners_Nb, [], 1)  ...
        ];

        % Calculate values using weights and face normals
        GamCSVals = reshape(Nfs*wOwnNeis, [], 1);

        % Create and add additional rows to GamCS
        GamCSNewRows = sparse(GamCSRows, GamCSCols, GamCSVals, theNumberOfFacesInPatch, 3*theNumberOfElements);

        GamCSAdj = [GamCSAdj; GamCSNewRows];

        %%% Adjust OmegaS
        newDim = size(OmegaSAdj,1) + theNumberOfFacesInPatch;

        newVols = [diag(OmegaSAdj); fVolSum/2];

        OmegaSAdj = sparse(1:newDim, 1:newDim, newVols, newDim, newDim);
        
        %%% Adjust M
        MRows = [owners_b; owners_Nb];

        MCols = repmat(1:theNumberOfFacesInPatch, 1, 2).';

        fAreas = (diag(Af(iBFaces, iBFaces)));
        MVals = [fAreas; -fAreas];

        MNewCols = sparse(MRows, MCols, MVals, theNumberOfElements, theNumberOfFacesInPatch);
        
        MAdj = [MAdj, MNewCols];
    end
end

% Create Lap for compact or wide stencil
if strcmp(LapStencil, 'compact')
    Lap = -MAdj*(OmegaSAdj\MAdj.');
    fprintf("Setting compact");
elseif strcmp(LapStencil, 'wide')
    McIn = MAdj*GamCSAdj;
    Lap = -McIn*(OmegaIn\McIn.');
    fprintf("Setting wide");
else
    error('set LapStencil in fvSolutions to  \"compact\" or \"wide\".')
end

% Trick to make symmetric
minDiag = abs(min(diag(Lap)));
Lap = Lap.*(abs(Lap)>tol*minDiag);
Lap = (Lap + Lap.')/2.0;

% Make remaining adjustments to Lap
addSource = sparse(theNumberOfElements, 1);
pRefRequired = true;

for iBPatch = 1:theNumberOfBoundaryPatches
    BType = Region.fluid.p.boundaryPatchRef{iBPatch}.type;

    iBFaces = cfdGetBFaceIndicesForBoundaryPatch(iBPatch);
    owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);    

    if strcmp(BType, 'fixedValue') || strcmp(BType, 'noSlip')
        if strcmp(LapStencil, 'wide')
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

    elseif strcmp(BType, 'zeroGradient') ...
            || strcmp(BType, 'cyclic') ...
            || strcmp(BType, 'slip') ...
            || strcmp(BType, 'outlet') ...
            || strcmp(BType, 'symmetry') ...
            || strcmp(BType, 'empty')

        % Boundary connections already excluded
        % Cyclic connections already set
        % Diagonal of Lap already adjusted accordingly
        % RHS b adds 0
    end
end


% Pressure reference cell and value
if pRefRequired
    iCell = Region.foamDictionary.fvSolution.AlguFVM.pRefCell + 1;
    value = Region.foamDictionary.fvSolution.AlguFVM.pRefValue;

   addSource(iCell) = addSource(iCell) + Lap(iCell,iCell)*value; 
   Lap(iCell,iCell) = 2*Lap(iCell,iCell);
end

% Store
Region.operators.Pois.Lap = Lap;
Region.operators.Pois.addSource =  addSource;