%AlgebraicAdjustment
function cfdReadTransportProperties

global Region;

fprintf('\nReading Transport Properties ...');

caseDirectoryPath = cfdGetCaseDirectoryPath;

transportPropertiesFile = [caseDirectoryPath, filesep, 'constant', filesep, 'transportProperties'];

% Check if "transportProperties" exists
if exist(transportPropertiesFile, 'file')~=2
    return;
end

% Open transportProperties file in read mode
fid = fopen(transportPropertiesFile, 'r');

% Collect constant cfdTransport properties from file
properties = cfdGetAllEntries(fid);

% Get mesh
theMesh = cfdGetMesh;

% Store attributes in the foam data base and create property fields
for iProperty=1:length(properties)
    
    propertyName = properties{iProperty}.key;
    
    % Skip if the entry is the transport model type
    if strcmp(propertyName, 'transportModel')
        continue;
    end
    
    fprintf(['\nReading ', propertyName,' ...']);
    
    duplicateIndex = strfind(properties{iProperty}.value, properties{iProperty}.key);
    if ~isempty(duplicateIndex)
        properties{iProperty}.value(duplicateIndex:duplicateIndex+length(properties{iProperty}.key)-1) = '';
    end
    
    if contains(properties{iProperty}.value, '[') && contains(line, ']')
        C = textscan(properties{iProperty}.value, '[%d %d %d %d %d %d %d] %f');
    
        % Value
        propertyValue = C{8};
       
        % Dimensions
        propertyDimensions = [C{1} C{2} C{3} C{4} C{5} C{6} C{7}];
        if isempty(propertyDimensions)
            error('%s does not have dimensions\n');
        end
    else
        C = textscan(properties{iProperty}.value, '%f');
        propertyValue = C{1};
        propertyDimensions = [];
    end

    if isempty(propertyValue)
        error('%s does not have a value\n');
    end        
    
    Region.foamDictionary.transportProperties.(propertyName).name = propertyName;
    Region.foamDictionary.transportProperties.(propertyName).dimensions = propertyDimensions;
    Region.foamDictionary.transportProperties.(propertyName).propertyValue = propertyValue;
end

propertyNames = cell(length(properties),1);
for iProp=1:length(properties)
    propertyNames{iProp} = properties{iProp}.key;
end

% Set mode
Region.compressible = false;
