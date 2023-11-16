function cfdUpdateAllFieldBoundaries

global Region

fields = fieldnames(Region.fluid);

for iField = 1:numel(fields)
    fieldName = fields{iField};

    fieldDims = size(Region.fluid.(fieldName).phi, 2);

    %surface fields never have to be updated because of the use of ghost
    %cells
    if ~contains(Region.fluid.(fieldName).type, 'surface')
        if fieldDims == 1
            cfdUpdateScalarFieldForAllBoundaryPatches(fieldName);
        else
            cfdUpdateVectorFieldForAllBoundaryPatches(fieldName);
        end
    end
end