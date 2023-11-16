function cfdStatistics

global Region;

op = Region.operators;

deltaT = cfdGetDeltaT;

p = cfdGetField('p');
U = cfdGetField('U');
Uf = cfdGetField('Uf');

%Calc checkerboarding
pLcp = -diag(op.OmegaIn)'*cfdGetInternalField(op.Gc*p, 'vvf').^2;
pLp = -diag(op.OmegaSIn)'*(op.G*p).^2;

if ~pLp == 0
    Ccb = 1 - pLcp/pLp;
else
    Ccb = 0;
end

% Calculate kinetic energy
EKin = 0.5*cfdGetInternalField(U, 'vvf')'*op.OmegaIn*cfdGetInternalField(U, 'vvf')/sum(diag(op.OmegaCIn));

%Calc Courant
CoField = 0.5*deltaT*cfdGetInternalField(op.OmegaC\(op.Tfc*abs(op.Af*Uf)), 'vsf');
maxCo = max(CoField);
meanCo = sum(op.OmegaCIn*CoField)/sum(diag(op.OmegaCIn));

fprintf('\nCcb: %4.3f, EKin: %10.9f, maxCo: %4.3f, meanCo: %4.3f\n', Ccb, EKin, maxCo, meanCo);
