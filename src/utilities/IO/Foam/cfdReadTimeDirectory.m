%AlgebraicAdjustment
function cfdReadTimeDirectory(varargin)

global Region;

timeDirectory = '';
readFields = [];

% Set arguments, both optional, but timeDirectory should come first
if nargin >= 1
    if ischar(varargin{1})
        timeDirectory = varargin{1};
    elseif iscell(varargin{1})
        readFields = string(varargin{1});
    end
end

if nargin == 2
    if iscell(varargin{2})
        readFields = string(varargin{2});
    end
end

% Set timeDirectory if not set previously
if strcmp(timeDirectory, '')
    controlDict = Region.foamDictionary.controlDict;
    
    if strcmp(controlDict.startFrom, 'startTime')
        timeDirectory = num2str(controlDict.startTime);
    elseif strcmp(controlDict.startFrom, 'firstTime')
        timeDirectory = '0';
    elseif strcmp(controlDict.startFrom, 'latestTime')
        timeDirectories = cfdGetTimeSteps;
        timeDirectory = num2str(max(timeDirectories));
    else
        error('\nError in controlDict. startFrom is not valid\n');
    end
end

% Case directory
caseDirectoryPath = cfdGetCaseDirectoryPath;

% Read field
files = dir([caseDirectoryPath, filesep, timeDirectory]);

if isempty(files)
    error('\n%s\n', ['Fields are not found in the ', timeDirectory,' directory']);
end

% Get mesh
theMesh = cfdGetMesh;
theNumberOfElements = cfdGetNumberOfElements;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;

for iFile=1:length(files)
    if (files(iFile).bytes)==0 || (files(iFile).isdir)
        continue;
    end

    if ~isempty(readFields)
        if ~ismember(files(iFile).name, readFields)
            continue;
        end
    end
    
    % get field name from file name
    fieldName = files(iFile).name;
    fieldFilePath = [caseDirectoryPath, filesep, timeDirectory, filesep, fieldName];
    fieldName = erase(fieldName, ["(", ")", ","]);
    % Open field in read mode
    fid = fopen(fieldFilePath, 'r');
    
    % Initialize header
    header = cell(0);
    
    % Scan/Read header
    while ~feof(fid)
        tline = fgetl(fid);
        
        % Skip empty lines
        tline = cfdSkipEmptyLines(fid, tline);
        
        % Skip macro-commented section
        tline = cfdSkipMacroComments(fid, tline);
        
        % Skip commented lines
        tline = cfdSkipCommentedLine(fid, tline);
        
        % read header block
        if cfdContains(tline, 'FoamFile')
            if isempty(header)
                header = cfdReadCfdDictionary(fid, 'FoamFile');
            else
                break;
            end
        else
            if ~isempty(header)
                break;
            end
        end
    end
    
    % Setup mesh field
    cfdSetupMeshField(fieldName, header.class);
    
    % Store initial value in mesh field
    theMeshField = cfdGetMeshField(fieldName);
    
    % Read dimensions
    theMeshField.dimensions = cfdGetKeyValue('dimensions', 'dimensions', fid);
    
    % Read and store internal field
    internalField = cfdGetKeyValue('internalField', 'string', fid);
    C = textscan(internalField, '%s', 1);
    valueType = C{1}{1};
    
    if strcmp(theMeshField.type, 'surfaceScalarField')
        
        if strcmp(valueType, 'uniform')
            value_str = textscan(internalField, 'uniform %f;', 1);
            theMeshField.phi(1:theNumberOfInteriorFaces) = value_str{1};
        elseif strcmp(valueType, 'nonuniform')
            theMeshField.phi(1:theNumberOfInteriorFaces) = cfdReadNonuniformList('internalField', fieldFilePath);
        end
        
        % Read and store bounday field
        boundaries = theMesh.cfdBoundaryPatchesArray;
        theNumberOfBPatches = length(boundaries);
        
        for iBPatch=1:theNumberOfBPatches
            
            % Get info
            theBCInfo = boundaries{iBPatch};
            iFaceStart = theBCInfo.startFaceIndex;
            iFaceEnd = iFaceStart+theBCInfo.numberOfBFaces-1;
            
            
            % Get field bounday condition
            type = cfdGetKeyValueFromBlock('type', ['boundaryField/', theBCInfo.name], fieldFilePath);
            
            % Store field bounday condition
            boundaryPatchRef.type = type;
            
            % Read bounday value
            value = cfdGetKeyValueFromBlock('value', ['boundaryField/', theBCInfo.name], fieldFilePath);
            
            if isempty(value) || strcmp(type, 'empty')
                boundaryPatchRef.value = 0;
                
                theMeshField.boundaryPatchRef{iBPatch} = boundaryPatchRef;
                
                % Clear boundaryPatchRef
                structureFields = fieldnames(boundaryPatchRef);
                for structureFieldName=structureFields
                    boundaryPatchRef = rmfield(boundaryPatchRef, structureFieldName);
                end
                
                continue;
            end
            
            C = textscan(value, '%s', 1);
            valueType = C{1}{1};
            
            boundaryPatchRef.valueType = valueType;
            
            if strcmp(valueType, 'uniform')
                value_str = textscan(value, 'uniform %f', 1);
                theMeshField.phi(iFaceStart:iFaceEnd) = value_str{1};
                boundaryPatchRef.value = value_str{1};
            elseif strcmp(valueType, 'nonuniform')
                valueList = cfdReadNonuniformList(theBCInfo.name, fieldFilePath);
                theMeshField.phi(iFaceStart:iFaceEnd) = valueList;
                boundaryPatchRef.value = valueList;
            end

            
            theMeshField.boundaryPatchRef{iBPatch} = boundaryPatchRef;

            % Clear boundaryPatchRef
            structureFields = fieldnames(boundaryPatchRef);
            for structureFieldName=structureFields
                boundaryPatchRef = rmfield(boundaryPatchRef, structureFieldName);
            end            
        end
        
    else
        
        if strcmp(valueType, 'uniform')
            if strcmp(theMeshField.type, 'volScalarField')
                value_str = textscan(internalField, 'uniform %f;', 1);
                theMeshField.phi(1:theNumberOfElements) = value_str{1};
            elseif strcmp(theMeshField.type, 'volVectorField')
                value_str = textscan(internalField, 'uniform (%f %f %f);', 1);
                for iElement=1:theNumberOfElements
                    theMeshField.phi(iElement,:) = [value_str{1}, value_str{2}, value_str{3}];
                end
            end
        elseif strcmp(valueType, 'nonuniform')
            theMeshField.phi(1:theNumberOfElements, :) = cfdReadNonuniformList('internalField', fieldFilePath);
        end
        
        % Read and store bounday field
        boundaries = theMesh.cfdBoundaryPatchesArray;
        theNumberOfBPatches = length(boundaries);
        
        for iBPatch=1:theNumberOfBPatches
            
            % Get info
            theBCInfo = boundaries{iBPatch};

            numberOfBFaces = theBCInfo.numberOfBFaces;
            iFaceStart = theBCInfo.startFaceIndex;
            
            iElementStart = theNumberOfElements+iFaceStart-theNumberOfInteriorFaces;
            iElementEnd = iElementStart+numberOfBFaces-1;
            
            
            % Get field bounday condition
            type = cfdGetKeyValueFromBlock('type', ['boundaryField/', theBCInfo.name], fieldFilePath);
            
            % Store field bounday condition
            boundaryPatchRef.type = type;
            
            % Read bounday value
            value = '';
            if ~strcmp(type, 'zeroGradient') ...
                && ~strcmp(type, 'slip') ...
                && ~strcmp(type, 'noSlip') ...
                && ~strcmp(type, 'outlet') ...
                && ~strcmp(type, 'symmetry') ...
                && ~strcmp(type, 'empty') ...
                && ~strcmp(type, 'cyclic')
                    value = cfdGetKeyValueFromBlock('value', ['boundaryField/', theBCInfo.name], fieldFilePath);
            end

            if isempty(value)
                if strcmp(theMeshField.type, 'volScalarField') || strcmp(theMeshField.type, 'surfaceScalarField')
                    boundaryPatchRef.value = 0;
                elseif strcmp(theMeshField.type, 'volVectorField')
                    boundaryPatchRef.value = [0,0,0];
                end
                
                theMeshField.boundaryPatchRef{iBPatch} = boundaryPatchRef;
                
                % Clear boundaryPatchRef
                structureFields = fieldnames(boundaryPatchRef);
                for structureFieldName=structureFields
                    boundaryPatchRef = rmfield(boundaryPatchRef, structureFieldName);
                end
                
                continue;
            end
            
            C = textscan(value, '%s', 1);
            valueType = C{1}{1};
            
            boundaryPatchRef.valueType = valueType;
            
            if strcmp(valueType, 'uniform')
                if strcmp(theMeshField.type, 'volScalarField') || strcmp(theMeshField.type, 'surfaceScalarField')
                    value_str = textscan(value, 'uniform %f', 1);
                    theMeshField.phi(iElementStart:iElementEnd) = value_str{1};
                    boundaryPatchRef.value = value_str{1};
                elseif strcmp(theMeshField.type, 'volVectorField')
                    value_str = textscan(value, 'uniform (%f %f %f)', 1);
                    for iBElement=iElementStart:iElementEnd
                        theMeshField.phi(iBElement,:) = [value_str{1}, value_str{2}, value_str{3}];
                    end
                    boundaryPatchRef.value = [value_str{1}, value_str{2}, value_str{3}];
                end
            elseif strcmp(valueType, 'nonuniform')
                valueList = cfdReadNonuniformList(theBCInfo.name, fieldFilePath);
                theMeshField.phi(iElementStart:iElementEnd, :) = valueList;
                boundaryPatchRef.value = valueList;
            end
            
            theMeshField.boundaryPatchRef{iBPatch} = boundaryPatchRef;
        
            % Clear boundaryPatchRef
            structureFields = fieldnames(boundaryPatchRef);
            for structureFieldName=structureFields
                boundaryPatchRef = rmfield(boundaryPatchRef, structureFieldName);
            end
        end
    end
    
    % Store mesh field in data base
    cfdSetMeshField(theMeshField);
    
    if ~contains(theMeshField.type, 'surface')
        if contains(theMeshField.type, 'Scalar')
            cfdUpdateScalarFieldForAllBoundaryPatches(theMeshField.name);
        else
            cfdUpdateVectorFieldForAllBoundaryPatches(theMeshField.name);
        end
    end    
end






