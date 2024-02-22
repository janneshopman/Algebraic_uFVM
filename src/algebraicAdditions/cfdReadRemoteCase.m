function caseStruct = cfdReadRemoteCase(caseDirectory, timeDirectory)
    if ~exist('caseDirectory', 'var')
        caseDirectory = pwd;
    end

    if ~exist('timeDirectory', 'var')
        files = dir(caseDirectory);
        i = 1;
        for iFile=1:length(files)
            fileName = files(iFile).name;
            if files(iFile).isdir && ~strcmp(fileName, '.') && ~strcmp(fileName, '..')
                if ~isnan(str2double(fileName))
                    timeSteps(i, 1) = eval(fileName);
                    i = i + 1;
                end
            end
        end
        
        timeSteps = sort(timeSteps, 1, 'ascend');
        timeDirectory = num2str(max(timeSteps));
    end

    global Region;
    
    if ~isempty(Region)
        global RegionStored;
        RegionStored = Region;
        clear global Region;
        global Region;
    end
       
    Region.caseDirectoryPath = caseDirectory;
    
    cfdReadPolyMesh;

    cfdReadTimeDirectory(timeDirectory);
    
    caseStruct = Region;
    
    clear global Region;
    
    global RegionStored;
    
    if ~isempty(RegionStored)
        global Region 
        Region = RegionStored;
    end

    clear global RegionStored;    
end