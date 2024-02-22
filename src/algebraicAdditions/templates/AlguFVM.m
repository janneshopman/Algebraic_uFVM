% Initialize
cfdStartSession;

cfdClearResults;

cfdReadOpenFoamFiles;

cfdSetOperators;

cfdSetPoisson;

cfdCreateFields;

cfdSetRKScheme;

cfdInitTime;

% Run Case
cfdRunCase;
