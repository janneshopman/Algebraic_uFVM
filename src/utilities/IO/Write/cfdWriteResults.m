%AlgebraicAdjustment
function cfdWriteResults

% Get time quantities
runTime = cfdGetCurrentTime;
cpuTime = cfdGetCPUTime;
nTimeStep = cfdGetnTimeStep;
deltaT = cfdGetDeltaT;


% Get control settings
foamDict = cfdGetFoamDict;
writeControl = foamDict.controlDict.writeControl;
writeInterval = foamDict.controlDict.writeInterval;
purgeWrite = foamDict.controlDict.purgeWrite;

% Get equation names
theFieldNames = cfdGetFields;

if cfdIsTransient
    if strcmp(writeControl, 'timeStep')
        if mod(nTimeStep, writeInterval)==0
            if ~cfdIsFolderExists(num2str(runTime))
                mkdir(num2str(runTime));
            end
            
            % Write all fields
            for iField=1:length(theFieldNames)
                if strcmp(theFieldNames{iField}, 'DU1') || ...
                        strcmp(theFieldNames{iField}, 'DU2') || ...
                        strcmp(theFieldNames{iField}, 'DU3') || ...
                        strcmp(theFieldNames{iField}, 'DUT1') || ...
                        strcmp(theFieldNames{iField}, 'DUT2') || ...
                        strcmp(theFieldNames{iField}, 'DUT3') || ...
                        strcmp(theFieldNames{iField}, 'pp')
                    continue;
                end
                cfdWriteOpenFoamField(theFieldNames{iField}, runTime);
            end
        end
    elseif strcmp(writeControl, 'runTime')
        if floor(runTime/writeInterval) == floor((runTime - deltaT)/writeInterval) + 1
            if ~cfdIsFolderExists(num2str(runTime))
                mkdir(num2str(runTime));
            end
            
            % Write all fields
            for iField=1:length(theFieldNames)
                if strcmp(theFieldNames{iField}, 'DU1') || ...
                        strcmp(theFieldNames{iField}, 'DU2') || ...
                        strcmp(theFieldNames{iField}, 'DU3') || ...
                        strcmp(theFieldNames{iField}, 'DUT1') || ...
                        strcmp(theFieldNames{iField}, 'DUT2') || ...
                        strcmp(theFieldNames{iField}, 'DUT3') || ...
                        strcmp(theFieldNames{iField}, 'pp')
                    continue;
                end
                cfdWriteOpenFoamField(theFieldNames{iField}, runTime);
            end
        end   
    elseif strcmp(writeControl, 'cpuTime')
        error('write on cpuTime not implemented')
    end
else
    error('not implemented for cfdTransient == false')
end

% Check the number of written time steps and delete first time steps
writtenTimeSteps = cfdGetTimeSteps;
if purgeWrite~=0 && length(writtenTimeSteps)>purgeWrite+1
    for iTime=2:length(writtenTimeSteps)-2
        rmdir(num2str(writtenTimeSteps(iTime)),'s');
    end
end

