%AlgebraicAdjustment
function list = cfdReadNonuniformList(key, fileDirectory, varargin)

% Read File
fileID = fopen(fileDirectory, 'r');

if nargin==2
    cfdSkipHeader = true;
else
    cfdSkipHeader = varargin{1};
end

% Initialize header
header = cell(0);

% Skip Header
if cfdSkipHeader
    % Scan/Read header
    while ~feof(fileID)         
        tline = fgetl(fileID);

        % Skip empty lines
        tline = cfdSkipEmptyLines(fileID, tline);

        % Skip macro-commented section
        tline = cfdSkipMacroComments(fileID, tline);

        % Skip commented lines
        tline = cfdSkipCommentedLine(fileID, tline);    

        % read header block
        if cfdContains(tline, 'FoamFile')   
            if isempty(header)
                header = cfdReadCfdDictionary(fileID, 'FoamFile');
            else
                break;
            end
        else
            if ~isempty(header)
                break;
            end        
        end
    end  
end

while(~feof(fileID))
    % Read each line
    tline = fgetl(fileID);
    
    % Skip empty lines
    if isempty(tline)
        continue;
    end
    
    % Skip commented lines
    if length(tline)>1
        if strcmp(tline(1:2), '//')
            continue;
        end
    end
    
    if strcmp(key, 'internalField')
        C = textscan(tline, '%s nonuniform List%s', 1);
        
        if strcmp(C{1}{1}, key)
            % Read class (scalar or cfdVector)
            fieldClass = C{2}{1};
            
            % Read list length
            C = textscan(fgetl(fileID), '%d', 1);
            listLength = C{1};
            
            fgetl(fileID);
            
            if strcmp(fieldClass, '<scalar>')
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '%f', 1);
                    list(i, 1) = C{1};
                end
            else
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '(%f %f %f)');
                    list(i, 1) = C{1};
                    list(i, 2) = C{2};
                    list(i, 3) = C{3};
                end
            end
        end
    else
        C = textscan(tline, '%s');
        if strcmp(C{1}{1}, key)
            
            while ~strcmp(C{1}{1}, 'value')
                line = fgetl(fileID);
                if contains(line, '(') && contains(line, ')')
                    C = textscan(line, '%s nonuniform List%s %d()', 1);
                    error('can''t read inline lists')
                else
                    C = textscan(line, '%s nonuniform List%s', 1);
                end
            end
            
            % Read class (scalar or cfdVector)
            fieldClass = C{2}{1};            
            
            % Read list length
            C = textscan(fgetl(fileID), '%d', 1);
            listLength = C{1};

            fgetl(fileID);
            
            if strcmp(fieldClass, '<scalar>')
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '%f', 1);
                    list(i, 1) = C{1};
                end
            else
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '(%f %f %f)');
                    list(i, 1) = C{1};
                    list(i, 2) = C{2};
                    list(i, 3) = C{3};
                end
            end
        end
    end
end

fclose(fileID);

