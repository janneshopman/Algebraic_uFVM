function cfdRKStages

global Region;

op = Region.operators;

RK = Region.foamDictionary.fvSchemes.ddtSchemes.RK;

sol = Region.foamDictionary.fvSolution.solvers.p;

theNumberOfFaces = cfdGetNumberOfFaces;

deltaT = cfdGetDeltaT;

nu = Region.foamDictionary.transportProperties.nu.propertyValue;

p = cfdGetField('p');
U = cfdGetField('U');
Uf = cfdGetField('Uf');

%reset for RK scheme
UOld = U;
pOld = p;
dUs = zeros(size(U, 1), RK.nStages);

for iStage = 1:RK.nStages + 1
    %stage increment of U
    %Is there a better way to include a pressure prediction?
    %Perhaps using updated intermediate values of pressure?
    U = UOld + dUs*RK.aTab(iStage, :)'; % - sum(RK.aTab(iStage,:)) * deltaT * pnPredCoef * op.Gc * pOld (3)
    cfdBCUpdate(U, 'U');

    if ~((iStage == 1) && (RK.aTab(1, 1) == 0))     %skip stage if first & explicit

        divU = cfdGetInternalField(op.Mc * U, 'vsf');

        %Intermediate p used for predictor, is it optimal?
        dtp = deltaT * cfdGetInternalField(p, 'vsf');

        %Before Poisson equation, b has to be updated according to bc (1)
        %And reference pressure has to be set (2)
        %Include preconditioner (4)
        [dtp, flag, relres, iter, resvec] = pcg(-op.Lap, -divU, sol.tolerance, sol.maxIter, speye(size(dtp, 1)), speye(size(dtp, 1)), dtp);

        p = cfdSetInternalField(p, dtp/deltaT, 'vsf');
        p = cfdBCUpdate(p, 'p');

        %velocity update
        U = U - deltaT * op.Gc * p;
        U = cfdBCUpdate(U, 'U');

        if strcmp(Region.foamDictionary.fvSchemes.laplacianSchemes.default, 'compact')
            Uf = op.GamCS*U + deltaT * (op.GamCS*op.GamSC - speye(theNumberOfFaces)) * op.G * p;
        else
            Uf = op.GamCS*U;
        end

        %Experiment with mixed Laplacian
        %Us = op.GamCS*U+op.LMix*(op.GamCS*op.GamSC - speye(geo.nf))*op.G * dtp; (5)
    end

    if iStage <= RK.nStages
        %set convective
        Con = kron(eye(3), op.M*spdiags(Uf, 0, theNumberOfFaces, theNumberOfFaces)*op.PiCSM);

        %set delta U for this stage
        dUs(:, iStage) = -deltaT * (op.Omega\((Con - nu * kron(eye(3), op.L))*U));
    end
end

%p = pnPredCoef*pOld + p; (3)

cfdSetField(p, 'p');
cfdSetField(U, 'U');
cfdSetField(Uf, 'Uf');