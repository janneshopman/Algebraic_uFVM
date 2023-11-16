function cfdSetRKScheme

global Region;

fprintf('\nSetting RKScheme\n');

RKSchemeName = Region.foamDictionary.fvSchemes.ddtSchemes.default;

%prep RK tableau
if strcmp(RKSchemeName, 'Euler')
    aTab = zeros(2, 1);
    aTab(2, 1) = 1;
elseif strcmp(RKSchemeName, 'RK3')
    aTab = zeros(4, 3);
    aTab(2, 1) = 1/2;
    aTab(3, 1:2) = [-1, 2];
    aTab(4, :) = [1/6, 2/3, 1/6];
else
    error('ddtScheme not implemented')
end

nStages = size(aTab, 1) - 1;

Region.foamDictionary.fvSchemes.ddtSchemes.RK.aTab = aTab;
Region.foamDictionary.fvSchemes.ddtSchemes.RK.nStages = nStages;
