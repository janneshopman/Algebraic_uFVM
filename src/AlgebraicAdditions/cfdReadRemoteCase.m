function caseStruct = cfdReadRemoteCase(caseDirectory, timeDirectory)
    if ~exist('caseDirectory', 'var')
        caseDirectory = pwd;
    end

    if ~exist('timeDirectory', 'var')
        timeDirectories = cfdGetTimeSteps;
        timeDirectory = num2str(max(timeDirectories));
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