%AlgebraicAdjustment
function cfdReadOpenFoamFiles
cfdReadPolyMesh;

cfdReadControlDictFile;

cfdReadFvSchemesFile;

cfdReadFvSolutionFile;

cfdReadTimeDirectory;

cfdReadTransportProperties;

