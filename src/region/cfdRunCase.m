%AlgebraicAdjustment
function cfdRunCase

cfdStatistics;

global Region

% Start steady false transience loop
while (cfdDoTransientLoop)

    % Time settings
    cfdUpdateRunTime;
    
    % Copy current time fields to previous time fields
    cfdUpdatePrevTimeStep;
    
    % Print current simulation time
    cfdPrintCurrentTime;

    cfdComputeEigenvalues;

    cfdRKStages;

    cfdSelfAdaptiveTimestep;
    
    cfdStatistics;   

    cfdWriteResults;

    display(Region.time.nTimeStep)
end
