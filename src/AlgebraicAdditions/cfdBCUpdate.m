function field = BCUpdate(field, fieldName)

global Region;

fieldDims = size(Region.fluid.(fieldName).phi, 2);

meshField = cfdGetMeshField(fieldName);

meshField.phi(:) = field;

cfdSetMeshField(meshField);

if fieldDims == 1
    cfdUpdateScalarFieldForAllBoundaryPatches(fieldName);
else
    cfdUpdateVectorFieldForAllBoundaryPatches(fieldName);
end

field = cfdGetField(fieldName);