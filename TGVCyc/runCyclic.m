% Initialize
cfdStartSession;

cfdClearResults;

cfdReadOpenFoamFiles;

cfdSetOperators;

cfdSetPoisson;

cfdCreateFields;

cfdSetRKScheme;

cfdInitTime;

cfdSetTGVFields;

% Run Case
cfdRunCase;
