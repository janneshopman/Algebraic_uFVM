function dotField = cfdCellDot(vField1, vField2)
% Check if the length of the input vectors is divisible by 3
if mod(numel(vField1), 3) ~= 0 || mod(numel(vField2), 3) ~= 0
    error('Input vectors must have a length divisible by 3.');
end

% Calculate the length of one dimension
len1D = numel(vField1) / 3;

% Extract sub-vectors
vField1x = vField1(1:len1D);
vField1y = vField1(len1D+1:2*len1D);
vField1z = vField1(2*len1D+1:end);

vField2x = vField2(1:len1D);
vField2y = vField2(len1D+1:2*len1D);
vField2z = vField2(2*len1D+1:end);

% Perform dot product for each component
dotx = vField1x .* vField2x;
doty = vField1y .* vField2y;
dotz = vField1z .* vField2z;

% Sum the dot products
dotField = dotx + doty + dotz;
