% Initialize
cfdStartSession;

cfdClearResults;

cfdReadOpenFoamFiles;

cfdSetOperators;

cfdSetPoisson;

cfdCreateFields;

cfdSetRKScheme;

cfdSetEigenvalues;

cfdInitTime;

cfdSetTGVFields;

cfdPerturbateField;

% Run Case
cfdRunCase;

cfdPlotStats;