function localpcb

global Region;

op = Region.operators;

nc = cfdGetNumberOfElements;
nx = sqrt(nc);
ny = nx;

nf = cfdGetNumberOfInteriorFaces;
nBf = cfdGetNumberOfBFaces;
nfTot = nf + nBf;
ncTot = nc + nBf;

deltaT = cfdGetDeltaT;

p = cfdGetField('p');
pSmth = op.PiSC*op.PiCS*p;
pSmth(nc + (1:nf)) = p(nc + (1:nf));
U = cfdGetField('U');
Uf = cfdGetField('Uf');

% p  = [1:length(p )].';
% U  = [1:length(U )].';
% Uf = [1:length(Uf)].';

UOmGcp = cfdCellDot(cfdGetInternalField(U, 'vvf'), op.OmegaIn*cfdGetInternalField(op.Gc*p, 'vvf'));

pMcU = cfdGetInternalField(p, 'vsf').*cfdGetInternalField(op.Mc*U, 'vsf');

pLp = cfdGetInternalField(p, 'vsf').*cfdGetInternalField(op.L*p, 'vsf');
pLcp = cfdGetInternalField(p, 'vsf').*cfdGetInternalField(op.M*op.GamSCS*op.G*p, 'vsf');

TOL = 1E-12;
for i = 1:length(pLp)
    if pLp(i) < TOL
        ipcbLcL(i) = 0;
    else
        ipcbLcL(i) = 1 - pLcp(i)/pLp(i);
    end
end

intMagGp = cfdGetInternalField(op.PiCS'*op.OmegaS*(op.G*p).^2, 'vsf');
magGcp = cfdCellDot(cfdGetInternalField(op.Gc*p, 'vvf'), op.OmegaIn*cfdGetInternalField(op.Gc*p, 'vvf'));

TOL = 1E-12;
for i = 1:length(intMagGp)
    if intMagGp(i) < TOL
        ipcbMag(i) = 0;
    else
        ipcbMag(i) = 1 - magGcp(i)/intMagGp(i);
    end
end

intMagGpSmth = cfdGetInternalField(op.PiCS'*op.OmegaS*(op.G*pSmth).^2, 'vsf');
magGcpSmth = cfdCellDot(cfdGetInternalField(op.Gc*pSmth, 'vvf'), op.OmegaIn*cfdGetInternalField(op.Gc*pSmth, 'vvf'));

TOL = 1E-12;
for i = 1:length(intMagGpSmth)
    if intMagGpSmth(i) < TOL
        ipcbMagSmth(i) = 0;
    else
        ipcbMagSmth(i) = 1 - magGcpSmth(i)/intMagGpSmth(i);
    end
end

global Region
if mod(Region.time.nTimeStep, 100) == 0
    UMag = cfdCellDot(cfdGetInternalField(U, 'vvf'), cfdGetInternalField(U, 'vvf'));
    UMagMat    = reshape(UMag,    nx, ny);
    
    [maxvalue, idx] = max(ipcbMag);
    %p([idx-66, idx-62, idx+62, idx+66]) = [0.5, 0.5, 0.5, 0.5];
    pMat           = reshape(cfdGetInternalField(p, 'vsf'), nx, ny);
    pSmthMat       = reshape(cfdGetInternalField(pSmth, 'vsf'), nx, ny);

    intMagGpMat    = reshape(intMagGp,    nx, ny);
    magGcpMat      = reshape(magGcp,      nx, ny);

    UOmGcpMat      = reshape(UOmGcp,      nx, ny);
    pMcUMat        = reshape(pMcU,        nx, ny);
    ipcbLcLMat     = reshape(ipcbLcL,     nx, ny);
    ipcbMagMat     = reshape(ipcbMag,     nx, ny);
    ipcbMagSmthMat = reshape(ipcbMagSmth, nx, ny);

    % Line plots
    nodeCs = Region.mesh.nodeCentroids;

    Lx = max(nodeCs(:, 1)) - min(nodeCs(:, 1));
    Ly = max(nodeCs(:, 2)) - min(nodeCs(:, 2));
    Lz = max(nodeCs(:, 3)) - min(nodeCs(:, 3));

    dx = Lx/nx;

    horBegCoor = [   dx/2, Ly/2, 0];
    horEndCoor = [Lx-dx/2, Ly/2, 0];    
    horInd = cfdLineSampleIndices(Region.mesh, horBegCoor, horEndCoor, nx);

    % Create a figure
    figure(1);
    
    figx = 4;
    figy = 3;
    figi = 1;

    subplot(figx,figy,figi);
    imagesc(UOmGcpMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('UOmGcp');
    figi = figi + 1;
    
    subplot(figx,figy,figi);
    imagesc(pMcUMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('pMcU');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(ipcbLcLMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('ipcbLcL');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(intMagGpMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('intMagGp');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(magGcpMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('magGcp');
    figi = figi + 1;    

    subplot(figx,figy,figi);
    imagesc(ipcbMagMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('ipcbMag');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(ipcbMagSmthMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('ipcbMagSmth');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(UMagMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('UMag');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(pMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('p');
    figi = figi + 1;

    subplot(figx,figy,figi);
    imagesc(pSmthMat);
    colormap('jet'); 
    colorbar;
    axis equal;
    title('pSmth');
    figi = figi + 1;

%     p5Lp = ...
%     [ ...
%         -0.0105   -0.0123   -0.0133   -0.0141   -0.0148   -0.0155   -0.0162 ...
%         -0.0171   -0.0180   -0.0191   -0.0202   -0.0214   -0.0226   -0.0237 ...
%         -0.0246   -0.0252   -0.0254   -0.0250   -0.0240   -0.0223   -0.0198 ...
%         -0.0166   -0.0128   -0.0088   -0.0050   -0.0016    0.0008    0.0018 ...
%          0.0014   -0.0004   -0.0032   -0.0081 ...
%     ];
% 
%     subplot(figx,figy,figi);
%     plot(Region.fluid.p.phi(horInd));
%     hold on;
%     plot(p5Lp);
%     hold off;
%     title('p - horizontal midline');
%     figi = figi + 1;



    fprintf('pause');
end