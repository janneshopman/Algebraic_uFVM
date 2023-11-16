%AlgebraicAdjustment
function cfdUpdateRunTime

global Region;

Region.time.currentTime = Region.time.currentTime + Region.foamDictionary.controlDict.deltaT;
Region.time.nTimeStep = Region.time.nTimeStep + 1;
Region.time.cpuTime = toc;
