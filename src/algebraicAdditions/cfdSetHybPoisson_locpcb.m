function cfdSetHybPoisson_locpcb(pnPredCoefMat)

tol = 1E-9;     % Cleaning tolerance

global Region;

interpolation = Region.foamDictionary.fvSolution.AlguFVM.interpolation;

Af = Region.operators.Af;
Nf = Region.operators.Nf;
M = Region.operators.M;
OmegaS = Region.operators.OmegaS;
OmegaIn = Region.operators.OmegaIn;
GamCS = Region.operators.GamCS;
GamSC = Region.operators.GamSC;
PiCS = Region.operators.PiCS;

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBFaces = cfdGetNumberOfBFaces;
theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;

ncTot = theNumberOfElements + theNumberOfBFaces;
nfTot = theNumberOfInteriorFaces + theNumberOfBFaces;
internalElements = 1:theNumberOfElements;
internalFaces = 1:theNumberOfInteriorFaces;

% Get internal pnPredCoefMat
pnPredCoefMat3D = kron(eye(3), pnPredCoefMat);
pnPredCoefMatIn = pnPredCoefMat(internalElements, internalElements);
pnPredCoefMatIn3D = kron(eye(3), pnPredCoefMatIn);
pnPredCoefMatS = diag(PiCS*diag(pnPredCoefMat));
pnPredCoefMatSIn = pnPredCoefMatS(internalFaces, internalFaces);

pnPredCoefMatSAdj = pnPredCoefMatSIn;

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
        
        %%% Adjust pnPredCoefMatS
        newCoefs = [diag(pnPredCoefMatSAdj); diag(pnPredCoefMatS(iBFaces, iBFaces))];

        pnPredCoefMatSAdj = sparse(1:newDim, 1:newDim, newCoefs, newDim, newDim);

        %%% Adjust M
        MRows = [owners_b; owners_Nb];

        MCols = repmat(1:theNumberOfFacesInPatch, 1, 2).';

        fAreas = (diag(Af(iBFaces, iBFaces)));
        MVals = [fAreas; -fAreas];

        MNewCols = sparse(MRows, MCols, MVals, theNumberOfElements, theNumberOfFacesInPatch);
        
        MAdj = [MAdj, MNewCols];
    end
end

makeStencils = {'wide', 'compact'};

for iSt = 1:length(makeStencils)
    makeStencil = makeStencils{iSt};

    if strcmp(makeStencil, 'wide')
        McIn = MAdj*GamCSAdj;
        LapWtd = -McIn*pnPredCoefMatIn3D*(OmegaIn\McIn.');
    elseif strcmp(makeStencil, 'compact')
        LapWtd = -MAdj*pnPredCoefMatSAdj*(OmegaSAdj\MAdj.');
    end

    % Trick to make symmetric
    minDiag = abs(min(diag(LapWtd)));
    LapWtd = LapWtd.*(abs(LapWtd)>tol*minDiag);
    LapWtd = (LapWtd + LapWtd.')/2.0;
    
    for iBPatch = 1:theNumberOfBoundaryPatches
        BType = Region.fluid.p.boundaryPatchRef{iBPatch}.type;

        if strcmp(BType, 'fixedValue') || strcmp(BType, 'noSlip')
            error('BC not implemented for hybrid Lap')
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

    if strcmp(makeStencil, 'wide')
        LapWideWtd = LapWtd;
    elseif strcmp(makeStencil, 'compact')
        LapCompWtd = LapWtd;
    end
end

% Create interpolator for hybrid corrector
GamSCSWtd = GamCS*pnPredCoefMat3D*GamSC;
% Replace rows of boundary faces with 1's on diagonal
bFaceIndices = theNumberOfInteriorFaces + (1:theNumberOfBFaces);
GamSCSbFaces = speye(size(GamSCSWtd));
GamSCSWtd(bFaceIndices,:) = GamSCSbFaces(bFaceIndices,:);

LapHyb = Region.operators.Pois.LapComp + LapWideWtd - LapCompWtd; 

% Store
Region.operators.Pois.GamSCSWtd = GamSCSWtd;
Region.operators.Pois.LapWtd = LapWtd;
Region.operators.Pois.LapCompWtd = LapCompWtd;
Region.operators.Pois.LapWideWtd = LapWideWtd;
Region.operators.Pois.LapHyb = LapHyb;