%AlgebraicAdjustment
function proceed = cfdDoTransientLoop

global Region;

currentTime = Region.time.currentTime;
endTime = Region.foamDictionary.controlDict.endTime;

deltaT = cfdGetDeltaT;

if currentTime<endTime-deltaT*1E-6
    proceed = true;
else
    proceed = false;
end