%AlgebraicAdjustment
function cfdReadFvSolutionFile

fprintf('\nReading fvSolution file ...\n');

global Region;

caseDirectoryPath = cfdGetCaseDirectoryPath;

fvSolutionFileDirectory = [caseDirectoryPath, filesep, 'system', filesep, 'fvSolution'];

% Check if "fvSolution" exists
if exist(fvSolutionFileDirectory, 'file')~=2
    error('\n%s\n','"fvSolution" file is not found in "~foamDirectory', filesep, 'system"');
end

% Read dictionary
fvSolutionDict = cfdReadFoamDictFile(fvSolutionFileDirectory);

% Store and manage dictionary in data base
% Solvers
fieldNamesToSolve = fieldnames(fvSolutionDict.solvers);
for iField=1:length(fieldNamesToSolve)
    fieldName = fieldNamesToSolve{iField};
    
    % Store in data base
    keys = fieldnames(fvSolutionDict.solvers.(fieldName));
    for iEntry=1:length(keys)
        key = keys{iEntry};
        value = fvSolutionDict.solvers.(fieldName).(key);
        if strcmp(key, 'solver') || strcmp(key, 'preconditioner') || strcmp(key, 'smoother')
            Region.foamDictionary.fvSolution.solvers.(fieldName).(key) = value;
        else
            Region.foamDictionary.fvSolution.solvers.(fieldName).(key) = eval(value);
        end
    end
    
    % Additional default settings
    if ~isfield(Region.foamDictionary.fvSolution.solvers.(fieldName), 'maxIter')
        Region.foamDictionary.fvSolution.solvers.(fieldName).maxIter = 1000;
    end
    
    if strcmp(Region.foamDictionary.fvSolution.solvers.(fieldName).solver, 'GAMG')
        if ~isfield(Region.foamDictionary.fvSolution.solvers.(fieldName), 'nPreSweeps')
            Region.foamDictionary.fvSolution.solvers.(fieldName).nPreSweeps = 0;
        end
        
        if ~isfield(Region.foamDictionary.fvSolution.solvers.(fieldName), 'nPostSweeps')
            Region.foamDictionary.fvSolution.solvers.(fieldName).nPostSweeps = 2;
        end
        
        if ~isfield(Region.foamDictionary.fvSolution.solvers.(fieldName), 'nFinestSweeps')
            Region.foamDictionary.fvSolution.solvers.(fieldName).nFinestSweeps = 2;
        end
    end
end

% AlguFVM control
if isfield(fvSolutionDict, 'AlguFVM')
    entryNames = fieldnames(fvSolutionDict.AlguFVM);
    for iEntry=1:length(entryNames)
        if ~isstruct(fvSolutionDict.AlguFVM.(entryNames{iEntry}))
            if strcmp(entryNames{iEntry}, 'pRefCell')
                % N.B. if pRefCell is given for an OpenFOAM case, add 1
                Region.foamDictionary.fvSolution.AlguFVM.(entryNames{iEntry}) = eval(fvSolutionDict.AlguFVM.(entryNames{iEntry})) + 1;
            elseif strcmp(entryNames{iEntry}, 'pRefPoint')
                pointString = fvSolutionDict.AlguFVM.(entryNames{iEntry});
                i = strfind(pointString, '(');
                pointArray = [str2double(strsplit(pointString(i+1:end-1), ' '))]';
                Region.foamDictionary.fvSolution.AlguFVM.(entryNames{iEntry}) = pointArray;
            else
                try
                    Region.foamDictionary.fvSolution.AlguFVM.(entryNames{iEntry}) = eval(fvSolutionDict.AlguFVM.(entryNames{iEntry}));
                catch
                    Region.foamDictionary.fvSolution.AlguFVM.(entryNames{iEntry}) = fvSolutionDict.AlguFVM.(entryNames{iEntry});
                end
            end
        else
            % Residual control
            if strcmp(entryNames{iEntry}, 'residualControl')
                fieldNames = fieldnames(fvSolutionDict.AlguFVM.residualControl);
                for iField=1:length(fieldNames)
                    fieldName = fieldNames{iField};
                    if ~isempty(find(strcmp(fieldNamesToSolve,fieldName)))
                        Region.foamDictionary.fvSolution.AlguFVM.residualControl.(fieldName) = eval(fvSolutionDict.AlguFVM.residualControl.(fieldName));
                    else
                        Region.foamDictionary.fvSolution.AlguFVM.residualControl.(fieldName) = 1e-6;
                    end
                end
            end
        end
    end
    
    % Default reference cell index and value
    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'pRefCell')
        if isfield(Region.foamDictionary.fvSolution.AlguFVM, 'pRefPoint')
            pRefPoint = Region.foamDictionary.fvSolution.AlguFVM.pRefPoint';
            pRefCell = cfdLineSampleIndices(Region.mesh, pRefPoint, pRefPoint, 1);
        else
            % If no pRefCell can be generated, directly set to 1 (equivalent to OpenFOAM 0)
            disp('No pRefCell nor -Point given, setting pRefCell = 1')
            pRefCell = 1;
        end
        Region.foamDictionary.fvSolution.AlguFVM.pRefCell = pRefCell;
    end
    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'pRefValue')
        Region.foamDictionary.fvSolution.AlguFVM.pRefValue = 0;
    end
   
    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'perturbate')
        disp('No perturbate, disabled by default')
        Region.foamDictionary.fvSolution.AlguFVM.perturbate = false;
    end

    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'computeEV')
        disp('No computeEV, disabled by default')
        Region.foamDictionary.fvSolution.AlguFVM.computeEV = false;
    end

    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'enableSAT')
        disp('No enableSAT, disabled by default')
        Region.foamDictionary.fvSolution.AlguFVM.enableSAT = false;
    end

    if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'plotSATStats')
        disp('No plotSATStats, disabled by default')
        Region.foamDictionary.fvSolution.AlguFVM.plotSATStats = false;
    end    

    if Region.foamDictionary.fvSolution.AlguFVM.computeEV
        if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'showRegion')
            disp('No showRegion, enabled by default')
            Region.foamDictionary.fvSolution.AlguFVM.showRegion = true;
        end
        if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'sampleFreq')
            disp(['No sampleFreq given, ' ...
                'setting sampleFreq=200 by default'])
            Region.foamDictionary.fvSolution.AlguFVM.sampleFreq = 200;
        end
        if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'evalPercent')
            disp(['No evalPercent given, ' ...
                'setting evalPercent=0.0 by default'])
            Region.foamDictionary.fvSolution.AlguFVM.evalPercent = 0.0;
        end        
    end

    if Region.foamDictionary.fvSolution.AlguFVM.enableSAT
        if ~isfield(Region.foamDictionary.fvSolution.AlguFVM, 'fdT')
            disp('No fDT given, setting fdT=1.0 by default')
            Region.foamDictionary.fvSolution.AlguFVM.fdT = 1.0;
        end
    end



    % If eigenvalues are not computed, SAT is disabled regardless of the 
    % user input
    if ~Region.foamDictionary.fvSolution.AlguFVM.computeEV && ...
            Region.foamDictionary.fvSolution.AlguFVM.enableSAT
        Region.foamDictionary.fvSolution.AlguFVM.enableSAT = false;
        disp(['computeEV is disabled by the user, enableSAT is ' ...
            'disabled automatically'])
    end
    
    % If SAT is disabled, plotStats is disabled regardless of the user 
    % input
    if ~Region.foamDictionary.fvSolution.AlguFVM.enableSAT && ...
            Region.foamDictionary.fvSolution.AlguFVM.plotSATStats
        Region.foamDictionary.fvSolution.AlguFVM.plotSATStats = false;
        disp(['enableSAT is disabled by the user, plotSATStats is ' ...
            'disabled automatically'])
    end



else
    disp('fvSolutionDict.AlguFVM dictionary not found, creating default')
    Region.foamDictionary.fvSolution.AlguFVM = struct(...
    'scheme', 'ForwardEuler', ...
    'LapStencil', 'compact', ...
    'PWIM', true, ...
    'interpolation', 'volumetric', ...
    'pnPredCoef', 0, ...
    'pRefCell', 1, ...
    'pRefValue', 0, ...
    'perturbate', false, ...
    'computeEV', false ...
    );
end




% Relaxation factors
if isfield(fvSolutionDict, 'relaxationFactors')
    % Equations
    if isfield(fvSolutionDict.relaxationFactors, 'equations')
        fieldNames = fieldnames(fvSolutionDict.relaxationFactors.equations);
        for iField=1:length(fieldNames)
            fieldName = fieldNames{iField};
            if ~isempty(find(strcmp(fieldNamesToSolve,fieldName)))
                Region.foamDictionary.fvSolution.relaxationFactors.equations.(fieldName) = eval(fvSolutionDict.relaxationFactors.equations.(fieldName));
            else
                Region.foamDictionary.fvSolution.relaxationFactors.equations.(fieldName) = 1.0;
            end
        end
    end
    
    % Fields
    if isfield(fvSolutionDict.relaxationFactors, 'fields')
        fieldNames = fieldnames(fvSolutionDict.relaxationFactors.fields);
        for iField=1:length(fieldNames)
            fieldName = fieldNames{iField};
            if ~isempty(find(strcmp(fieldNamesToSolve,fieldName)))
                Region.foamDictionary.fvSolution.relaxationFactors.fields.(fieldName) = eval(fvSolutionDict.relaxationFactors.fields.(fieldName));
            else
                Region.foamDictionary.fvSolution.relaxationFactors.fields.(fieldName) = 1.0;
            end
        end
    end
end
