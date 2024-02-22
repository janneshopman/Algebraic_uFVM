function cfdUpdateAllFieldBoundaries

global Region

fields = fieldnames(Region.fluid);

for iField = 1:numel(fields)
    fieldName = fields{iField};

    fieldDims = size(Region.fluid.(fieldName).phi, 2);

    % Surface fields never updated because of ghost cells
    if ~contains(Region.fluid.(fieldName).type, 'surface')
        if fieldDims == 1
            cfdUpdateScalarFieldForAllBoundaryPatches(fieldName);
        else
            cfdUpdateVectorFieldForAllBoundaryPatches(fieldName);
        end
    end
end
