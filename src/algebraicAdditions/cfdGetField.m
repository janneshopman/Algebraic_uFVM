function theField = cfdGetField(theFieldName)

global Region;

if isfield(Region.fluid, theFieldName)
    theField = Region.fluid.(theFieldName);
else
    theField = -1;
end

theField = theField.phi(:);