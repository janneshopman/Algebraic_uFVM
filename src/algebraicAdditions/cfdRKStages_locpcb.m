function cfdRKStages_locpcb

global Region;

op = Region.operators;
fvSol = Region.foamDictionary.fvSolution;

RK = Region.foamDictionary.fvSchemes.ddtSchemes.RK;

if isfield(fvSol.solvers, 'pFinal')
    sol = fvSol.solvers.pFinal;
else
    sol = fvSol.solvers.p;
end

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBFaces = cfdGetNumberOfBFaces;
theNumberOfFaces = cfdGetNumberOfFaces;
ncTot = theNumberOfElements + theNumberOfBFaces;
owners = cfdGetOwners;

deltaT = cfdGetDeltaT;

nu = Region.foamDictionary.transportProperties.nu.propertyValue;

p = cfdGetField('p');
pCorr = cfdGetField('pCorr');
U = cfdGetField('U');
Uf = cfdGetField('Uf');
Gcp = cfdGetField('Gcp');

% Calculate global checkerboard coefficient
pLcp = -diag(op.OmegaIn)'*cfdGetInternalField(op.Gc*p, 'vvf').^2;
pLp = -diag(op.OmegaSIn)'*(op.G*p).^2;

if ~pLp == 0
    Ccb = 1 - pLcp/pLp;
else
    Ccb = 0;
end    

% Calculate local checkerboard coefficient
intMagGp = cfdGetInternalField(op.PiCS'*op.OmegaS*(op.G*p).^2, 'vsf');
magGcp = cfdCellDot(cfdGetInternalField(op.Gc*p, 'vvf'), op.OmegaIn*cfdGetInternalField(op.Gc*p, 'vvf'));

gammaCb = zeros(ncTot, 1);
iNonZero = abs(intMagGp) > 0;
gammaCb(iNonZero) = 1 - magGcp(iNonZero)./intMagGp(iNonZero);
% - Set value in ghost cells equal to boundary cells
%gammaCb(theNumberOfElements+(1:theNumberOfBFaces)) = gammaCb(owners(theNumberOfInteriorFaces+(1:theNumberOfBFaces)));
gammaCb = cfdUpdateScalarBocoLike(gammaCb, 'p');
GammaCb = spdiags(gammaCb, 0, ncTot, ncTot);

fprintf("avg(gammaCb): %f\n", mean(cfdGetInternalField(gammaCb, 'vsf')));

% Pressure predictor
pnPredCoef = fvSol.AlguFVM.pnPredCoef;

resetGcp = true;

if pnPredCoef > 1.0 
    nResets = 1;
    modWidth = nResets-1;
    
    if mod(Region.time.nTimeStep+modWidth, 100) > modWidth
        resetGcp = false;
    end

    if resetGcp
        pnPredCoefMat = sparse(ncTot,ncTot);
        fprintf("Reset pnPredCoefMat to scalar\n");
    else
        pnPredCoefMat = speye(ncTot) - GammaCb; 
    end    
else
    if pnPredCoef < 0
        pnPredCoefMatVal = 1 - Ccb;
    else
        pnPredCoefMatVal = pnPredCoef;
    end

    pnPredCoefMat = pnPredCoefMatVal*speye(ncTot);
end

% Lap width
LWidthCoef = fvSol.AlguFVM.LWidthCoef;

if LWidthCoef > 1.0    
    LWidthCoefMat = speye(ncTot) - GammaCb;
else
    if LWidthCoef < 0
        LWidthCoefMatVal = 1 - Ccb;
    else
        LWidthCoefMatVal = LWidthCoef;
    end

    LWidthCoefMat = LWidthCoefMatVal*speye(ncTot);
end

LWidthCoefMatS = diag(op.PiCS*diag(LWidthCoefMat));

cfdSetHybPoisson_locpcb(LWidthCoefMat);
op = Region.operators;
GamSCSWtd = op.Pois.GamSCSWtd;
Lap = op.Pois.LapHyb;

% Set pRef
addSource = op.Pois.addSource;

% Makes time-step dependent adjustments, so don't store changes
if Region.operators.Pois.pRefRequired
    iCell = Region.foamDictionary.fvSolution.AlguFVM.pRefCell;
    pCorrRefValue = Region.foamDictionary.fvSolution.AlguFVM.pRefValue - pnPredCoefMat(iCell,iCell)*p(iCell);

    addSource(iCell) = addSource(iCell) + Lap(iCell,iCell)*pCorrRefValue; 
    Lap(iCell,iCell) = 2*Lap(iCell,iCell);
end

% Reset values for RK stages
dtpCorr = deltaT * cfdGetInternalField(pCorr, 'vsf');
UOld = U;
dUs = zeros(size(U, 1), RK.nStages);

for iStage = 1:RK.nStages + 1
    % Calculate fractional time step
    cdt = sum(RK.aTab(iStage,:))*deltaT;

    % Stage increment of U
    % if second order pPred is used, than the value can be extrapolated for
    % stages in this line
    Ucp = UOld + dUs*RK.aTab(iStage, :)' - cdt * kron(eye(3), pnPredCoefMat) * Gcp;  
    Ucp = cfdBCUpdate(Ucp, 'Ucp');

    if ~((iStage == 1) && (RK.aTab(1, 1) == 0))     % Skip stage if (first & explicit)
        divUcp = cfdGetInternalField(op.Mc * Ucp, 'vsf');
        source = divUcp + cdt * addSource;

        % Include preconditioner (To do - 4)
        dtpCorr = cfdCheckpcg(-Lap, -source, sol.tolerance, sol.maxIter, speye(size(dtpCorr, 1)), speye(size(dtpCorr, 1)), dtpCorr);

        pCorr = cfdSetInternalField(pCorr, dtpCorr/cdt, 'vsf');
        pCorr = cfdBCUpdate(pCorr, 'pCorr');

        U = Ucp - cdt * op.Gc * pCorr;
        U = cfdBCUpdate(U, 'U');

        Uf = op.GamCS*Ucp - cdt * (speye(theNumberOfFaces) - LWidthCoefMatS + GamSCSWtd) * op.G * pCorr; 
    end

    if iStage <= RK.nStages
        % Set convective
        Con = kron(eye(3), op.M*spdiags(Uf, 0, theNumberOfFaces, theNumberOfFaces)*op.PiCSM);

        % Set delta U for this stage
        dUs(:, iStage) = -deltaT * (op.Omega\((Con - nu * kron(eye(3), op.L))*U));
    end
end

p = pnPredCoefMat*p + pCorr;
p = cfdBCUpdate(p, 'p');

if resetGcp
    Gcp = op.Gc*p;
else
    Gcp = kron(eye(3), pnPredCoefMat) * Gcp + op.Gc*pCorr;
end
Gcp = cfdBCUpdate(Gcp, 'Gcp');

cfdSetField(p, 'p');
cfdSetField(p, 'pCorr');
cfdSetField(U, 'U');
cfdSetField(Uf, 'Uf');
cfdSetField(Gcp, 'Gcp');