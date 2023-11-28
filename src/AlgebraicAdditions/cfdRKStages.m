function cfdRKStages

global Region;

op = Region.operators;
fvSol = Region.foamDictionary.fvSolution;

RK = Region.foamDictionary.fvSchemes.ddtSchemes.RK;

if isfield(fvSol.solvers, 'pFinal')
    sol = fvSol.solvers.pFinal;
else
    sol = fvSol.solvers.p;
end

pnPredCoef = fvSol.AlguFVM.pnPredCoef;
LapStencil = fvSol.AlguFVM.LapStencil;
PWIM = fvSol.AlguFVM.PWIM;

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfFaces = cfdGetNumberOfFaces;

deltaT = cfdGetDeltaT;

nu = Region.foamDictionary.transportProperties.nu.propertyValue;

p = cfdGetField('p');
pCorr = cfdGetField('pCorr');
U = cfdGetField('U');
Uf = cfdGetField('Uf');

% Reset for RK scheme
pPred = pnPredCoef * p;
dtpCorr = deltaT * cfdGetInternalField(pCorr, 'vsf');
UOld = U;
dUs = zeros(size(U, 1), RK.nStages);

for iStage = 1:RK.nStages + 1
    % Stage increment of U
    Ucp = UOld + dUs*RK.aTab(iStage, :)' - sum(RK.aTab(iStage,:)) * deltaT * op.Gc * pPred;
    Ucp = cfdBCUpdate(Ucp, 'Ucp');

    if ~((iStage == 1) && (RK.aTab(1, 1) == 0))     % Skip stage if first & explicit
        divUcp = cfdGetInternalField(op.Mc * Ucp, 'vsf');
        source = divUcp + deltaT*op.Pois.addSource;

        % Include preconditioner (To do - 4)
        dtpCorr = cfdCheckpcg(-op.Pois.Lap, -source, sol.tolerance, sol.maxIter, speye(size(dtpCorr, 1)), speye(size(dtpCorr, 1)), dtpCorr);

        pCorr = cfdSetInternalField(pCorr, dtpCorr/deltaT, 'vsf');
        pCorr = cfdBCUpdate(pCorr, 'pCorr');

        % Velocity update
        U = Ucp - deltaT * op.Gc * pCorr;
        Uf = op.GamCS*Ucp - deltaT * op.G * pCorr;

        if ~PWIM
            if strcmpi(LapStencil, 'compact')
                U = U + (op.GamSC*op.GamCS - speye(theNumberOfElements)) * Ucp;
            elseif strcmpi(LapStencil, 'wide')
                Uf = Uf + deltaT * (speye(theNumberOfFaces) - op.GamCS*op.GamSC) * op.G * pCorr;
            end
        end

        U = cfdBCUpdate(U, 'U');
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

cfdSetField(p, 'p');
cfdSetField(p, 'pCorr');
cfdSetField(U, 'U');
cfdSetField(Uf, 'Uf');