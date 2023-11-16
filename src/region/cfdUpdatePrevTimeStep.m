%AlgebraicAdjustment
function cfdUpdatePrevTimeStep

fieldNames = cfdGetFields;
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    field = cfdGetMeshField(fieldName);
    field.prevTimeStep.phi = field.phi;
    cfdSetMeshField(field);
end