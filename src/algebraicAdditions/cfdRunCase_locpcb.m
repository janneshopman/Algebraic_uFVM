%AlgebraicAdjustment
function cfdRunCase_locpcb

cfdStatistics;

% Start steady false transience loop
while (cfdDoTransientLoop)

    % Time settings
    cfdUpdateRunTime;
    
    % Copy current time fields to previous time fields
    cfdUpdatePrevTimeStep;
    
    % Print current simulation time
    cfdPrintCurrentTime;

    cfdRKStages_locpcb;

    cfdStatistics;

    cfdWriteResults;
end
