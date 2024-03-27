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

LapStencil = fvSol.AlguFVM.LapStencil;
PWIM = fvSol.AlguFVM.PWIM;

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
% Set value in ghost cells equal to boundary cells
gammaCb(theNumberOfElements+(1:theNumberOfBFaces)) = gammaCb(owners(theNumberOfInteriorFaces+(1:theNumberOfBFaces)));
GammaCb = spdiags(gammaCb, 0, ncTot, ncTot);

% Pressure predictor
pnPredCoef = fvSol.AlguFVM.pnPredCoef;
pnPredCoefMat = pnPredCoef*speye(ncTot);

% Set pnPredCoef and pnPredCoefMat
if pnPredCoef < 0.0    
    pnPredCoef = 1 - Ccb;

    pnPredCoefMat = speye(ncTot) - GammaCb;
end

fprintf("avg(gammaCb): %f\n", mean(cfdGetInternalField(gammaCb, 'vsf')));

altMeth = "LcL";

pPred   = 1*p;
GcpPred = pnPredCoef*Gcp;

Lap     = op.Pois.Lap;
LapWide = op.Pois.LapWide;
LapComp = op.Pois.LapComp;

if strcmp(altMeth, "pp")
    %pPred = p - (1 - pnPredCoef   ) * (speye(ncTot) - op.PiSC*op.PiCS)*p;
    %pPred = p - (speye(ncTot) - pnPredCoefMat) * (speye(ncTot) - op.PiSC*op.PiCS)*p;
    %pPred = pnPredCoefMat*p;
elseif strcmp(altMeth, "Gcp")
    pPred = 0*p;

    %GcpPred = kron(eye(3), pnPredCoefMat)*Gcp;
elseif strcmp(altMeth, "LcL")
    Lap = pnPredCoef*LapWide + (1-pnPredCoef)*LapComp;
elseif strcmp(altMeth, "LcLLoc")
    
end

% Set pRef
addSource = op.Pois.addSource;

% Makes time-step dependent adjustments, so don't store changes
if Region.operators.Pois.pRefRequired
    iCell = Region.foamDictionary.fvSolution.AlguFVM.pRefCell;
    pCorrRefValue = Region.foamDictionary.fvSolution.AlguFVM.pRefValue - pPred(iCell);

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
    if strcmp(altMeth, "Gcp")
        Ucp = UOld + dUs*RK.aTab(iStage, :)' - cdt * GcpPred;          
    else
        Ucp = UOld + dUs*RK.aTab(iStage, :)' - cdt * op.Gc * pPred;  
    end
    Ucp = cfdBCUpdate(Ucp, 'Ucp');

    if ~((iStage == 1) && (RK.aTab(1, 1) == 0))     % Skip stage if (first & explicit)
        divUcp = cfdGetInternalField(op.Mc * Ucp, 'vsf');
        source = divUcp + cdt * addSource;

        % Include preconditioner (To do - 4)
        dtpCorr = cfdCheckpcg(-Lap, -source, sol.tolerance, sol.maxIter, speye(size(dtpCorr, 1)), speye(size(dtpCorr, 1)), dtpCorr);

        pCorr = cfdSetInternalField(pCorr, dtpCorr/cdt, 'vsf');
        pCorr = cfdBCUpdate(pCorr, 'pCorr');

        % Velocity update
        if strcmp(LapStencil, 'wide')
            U = Ucp - cdt * op.Gc * pCorr;
            U = cfdBCUpdate(U, 'U');
    
            if PWIM
                Uf = op.GamCS*Ucp - cdt * op.G * pCorr;
            else
                Uf = op.GamCS*U;
            end

            if strcmp(altMeth, "LcL")
                Uf = op.GamCS(Ucp - pnPredCoef*cdt*op.GamSC*op.G*pCorr) - cdt*(1-pnPredCoef)*op.G*pCorr;
            end
        elseif strcmp(LapStencil, 'compact')
            Uf = op.GamCS*Ucp - cdt * op.G * pCorr;
    
            if PWIM
                U = Ucp - cdt * op.Gc * pCorr;
            else
                U = op.GamSC*Uf;
            end
            U = cfdBCUpdate(U, 'U');
        else
            error('LapStencil should be set to wide or compact')
        end
    end

    if iStage <= RK.nStages
        % Set convective
        Con = kron(eye(3), op.M*spdiags(Uf, 0, theNumberOfFaces, theNumberOfFaces)*op.PiCSM);

        % Set delta U for this stage
        dUs(:, iStage) = -deltaT * (op.Omega\((Con - nu * kron(eye(3), op.L))*U));
    end
end

p = pPred + pCorr;
p = cfdBCUpdate(p, 'p');

%Gcp = GcpPred + op.Gc*pCorr;
Gcp = op.Gc*p;
Gcp = cfdBCUpdate(Gcp, 'Gcp');

cfdSetField(p, 'p');
cfdSetField(p, 'pCorr');
cfdSetField(U, 'U');
cfdSetField(Uf, 'Uf');
cfdSetField(Gcp, 'Gcp');