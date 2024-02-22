function field = cfdSetInternalField(field, inField, type)
nb = cfdGetNumberOfBFaces;

nc = cfdGetNumberOfElements;
ncTot = nc + nb;

nf = cfdGetNumberOfInteriorFaces;
nfTot = nf + nb;

if strcmp(type, 'vsf')
    indices = 1:nc;
elseif strcmp(type, 'ssf')
    indices = 1:nf;
elseif strcmp(type, 'vvf')
    indices = [int32(1:nc), ncTot+int32(1:nc), 2*ncTot+int32(1:nc)];
elseif strcmp(type, 'svf')
    indices = [int32(1:nf), nfTot+int32(1:nf), 2*nfTot+int32(1:nf)];
else
    error('Unable to retrieve internal field')
end

field(indices) = inField;