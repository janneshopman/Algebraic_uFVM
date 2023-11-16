function cfdSetField(field, fieldName)

global Region;

Region.fluid.(fieldName).phi(:) = field;