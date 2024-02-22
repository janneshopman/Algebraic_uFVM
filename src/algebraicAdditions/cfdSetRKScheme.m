function cfdSetRKScheme

global Region;

fprintf('\nSetting RKScheme\n');

RKSchemeName = Region.foamDictionary.fvSolution.AlguFVM.scheme;

% Set up RK tableaus
if strcmp(RKSchemeName, 'ForwardEuler')
    aTab = zeros(2, 1);
    aTab(2, 1) = 1;
elseif strcmp(RKSchemeName, 'Kutta')
    aTab = zeros(4, 3);
    aTab(2, 1) = 1/2;
    aTab(3, 1:2) = [-1, 2];
    aTab(4, :) = [1/6, 2/3, 1/6];
else
    error('ddtScheme not implemented')
end

% Implicit schemes not implemented, so this holds
nStages = size(aTab, 1) - 1;

Region.foamDictionary.fvSchemes.ddtSchemes.RK.aTab = aTab;
Region.foamDictionary.fvSchemes.ddtSchemes.RK.nStages = nStages;
