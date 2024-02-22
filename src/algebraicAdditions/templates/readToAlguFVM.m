cfdStartSession;

cfdReadPolyMesh;

cfdReadControlDictFile;

cfdReadFvSchemesFile;

cfdReadFvSolutionFile;

readFields = {'p', 'U', 'phi'};

cfdReadTimeDirectory('0', readFields);

cfdReadTransportProperties;

cfdSetOperators;

cfdSetPoisson;

cfdCreateFields;

cfdStatistics;
