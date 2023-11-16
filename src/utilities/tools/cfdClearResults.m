%AlgebraicAdjustment
function cfdClearResults

% Clear time steps
theTimeSteps = cfdGetTimeSteps;
for timeStep=theTimeSteps'
    if cfdIsFolderExists(num2str(timeStep))
        if timeStep~=0
            rmdir(num2str(timeStep), 's');
        end
    end
end

if cfdIsFolderExists('0.orig')
    if cfdIsFolderExists('0')
        rmdir('0', 's');
    end
    copyfile('0.orig', '0')
end

% Clear other additional folders
if cfdIsFolderExists('convergence')
    rmdir('convergence', 's');
end
