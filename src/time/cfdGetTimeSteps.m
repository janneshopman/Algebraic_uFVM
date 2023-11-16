%AlgebraicAdjustment
function timeSteps = cfdGetTimeSteps

timeSteps = [];

files = dir;
i = 1;
for iFile=1:length(files)
    fileName = files(iFile).name;
    if files(iFile).isdir && ~strcmp(fileName, '.') && ~strcmp(fileName, '..')
        if ~isempty(str2num(fileName))
            timeSteps(i, 1) = eval(fileName);
            i = i + 1;
        end
    end
end

timeSteps = sort(timeSteps, 1, 'ascend');
