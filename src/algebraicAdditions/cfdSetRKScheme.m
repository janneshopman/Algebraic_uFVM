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
elseif strcmp(RKSchemeName, 'heunRK2')
    aTab = zeros(3, 2);
    aTab(2, 1) = 1;
    aTab(3, :) = [0.5, 0.5];
elseif strcmp(RKSchemeName, 'stdRK4')
    aTab = zeros(5, 4);
    aTab(2, 1) = 0.5;
    aTab(3, 1:2) = [0, 0.5];
    aTab(4, 1:3) = [0, 0, 1.0];
    aTab(5, :) = [1/6, 1/3, 1/3, 1/6];
elseif strcmp(RKSchemeName, 'nystromRK5')
    aTab = zeros(7, 6);
    aTab(2, 1) = 1/3;
    aTab(3, 1:2) = [4/25, 6/25];
    aTab(4, 1:3) = [1/4, -3, 15/4];
    aTab(5, 1:4) = [2/27, 10/9, -50/81, 8/81];
    aTab(6, 1:5) = [2/25, 12/25, 2/15, 8/75, 0];
    aTab(7, :) = [23/192, 0, 125/192, 0, -27/64, 125/192];
else
    error('ddtScheme not implemented')
end

% Implicit schemes not implemented, so this holds
nStages = size(aTab, 1) - 1;

Region.foamDictionary.fvSchemes.ddtSchemes.RK.aTab = aTab;
Region.foamDictionary.fvSchemes.ddtSchemes.RK.nStages = nStages;
