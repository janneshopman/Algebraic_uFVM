%AlgebraicAdjustment
function cfdInitTime

global Region;

Region.time.nTimeStep = 0;

if strcmp(Region.foamDictionary.controlDict.startFrom, 'startTime')
    Region.time.currentTime = Region.foamDictionary.controlDict.startTime;
elseif strcmp(Region.foamDictionary.controlDict.startFrom, 'latestTime')
    timeSteps = cfdGetTimeSteps;
    Region.time.currentTime = max(timeSteps);
elseif strcmp(Region.foamDictionary.controlDict.startFrom, 'firstTime')
    timeSteps = cfdGetTimeSteps;
    Region.time.currentTime = min(timeSteps);
end
