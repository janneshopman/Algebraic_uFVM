function cfdSetOperators

tol = 1E-9;     % Cleaning tolerance

global Region;

fprintf('\n\nCreating operators\n');

theNumberOfElements = cfdGetNumberOfElements;
theNumberOfFaces = cfdGetNumberOfFaces;
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
theNumberOfBFaces = cfdGetNumberOfBFaces;
theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;

ncTot = theNumberOfElements + theNumberOfBFaces;

owners = cfdGetOwners;
allNeighbours = cfdGetAllNeighbours;

Sfs = cfdGetFaceSf;
Sfs = Sfs.*(abs(Sfs)>tol);      % Clean small values

CcfsIn = Region.mesh.faceCF;
CcfsIn = CcfsIn.*(abs(CcfsIn)>tol);   % Clean small values

Cofs = Region.mesh.faceCf;
Cofs = Cofs.*(abs(Cofs)>tol);   % Clean small values

CVolsIn = cfdGetVolumesForElements;

% Erase connection to empty type ghost cells
% Set boundary-ghost distance and cyclic boundary pair distance
% Set ghost cell volumes
Tfo = sparse(owners, 1:length(owners), 1, ncTot, theNumberOfFaces);
Tfn = sparse(allNeighbours, 1:length(allNeighbours), 1, ncTot, theNumberOfFaces);

Tff = speye(theNumberOfFaces);

CVols = zeros(ncTot, 1);
CVols(1:theNumberOfElements) = CVolsIn;

for iBPatch=1:theNumberOfBoundaryPatches
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);

    iBFaces = cfdGetBFaceIndicesForBoundaryPatch(iBPatch);
    iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
    owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);    

    if strcmp(theBCInfo.type, 'cyclic')
        for iNBPatch=1:theNumberOfBoundaryPatches
            theNBCInfo = cfdGetBoundaryPatchRef(iNBPatch);
            if strcmp(theNBCInfo.name, theBCInfo.neighbourPatch)                  
                iNBFaces = cfdGetBFaceIndicesForBoundaryPatch(iNBPatch);
                owners_Nb = cfdGetOwnersSubArrayForBoundaryPatch(iNBPatch);  

                Region.mesh.cfdBoundaryPatchesArray{iBPatch}.neighbourPatchId = iNBPatch;

                break
            end
        end   

        Tff = Tff + sparse(iBFaces, iNBFaces, -1, theNumberOfFaces, theNumberOfFaces);        

        CVols(iBElements) = CVolsIn(owners_Nb);

    else
        if strcmp(theBCInfo.type, 'empty')
            Tfo(sub2ind(size(Tfo), owners_b, iBFaces)) = 0;
            Tfn(sub2ind(size(Tfn), iBElements, iBFaces)) = 0;        
        end
        
        Tff = Tff + sparse(iBFaces, iBFaces, 1, theNumberOfFaces, theNumberOfFaces);        

        CVols(iBElements) = CVolsIn(owners_b);
    end
end

Tfc = Tfo + Tfn;

Ccfs = Tff*CcfsIn;

% Construct geometric matrices
% N.B empty cells do have geometric values, just no connection
Sf = spdiags(Sfs, double(theNumberOfFaces)*[0, 1, 2], theNumberOfFaces, 3*theNumberOfFaces);
Af = spdiags(sqrt(sum(Sfs.^2,2)), 0, theNumberOfFaces, theNumberOfFaces);
Nf = Af\Sf;
Ccf = spdiags(Ccfs, double(theNumberOfFaces)*[0, 1, 2], theNumberOfFaces, 3*theNumberOfFaces);
Cof = spdiags(Cofs, double(theNumberOfFaces)*[0, 1, 2], theNumberOfFaces, 3*theNumberOfFaces);
Cnf = Ccf - Cof;
Dnf = Nf*Ccf';
Dnof = Nf*Cof';
Dnnf = Dnf - Dnof;
Wfs = diag(Dnof)./diag(Dnf);
Wof = spdiags(Wfs, 0, theNumberOfFaces, theNumberOfFaces);
Wnf = spdiags(1-Wfs, 0, theNumberOfFaces, theNumberOfFaces);
CcfIn = spdiags(CcfsIn, double(theNumberOfFaces)*[0, 1, 2], theNumberOfFaces, 3*theNumberOfFaces);
DnfIn = Nf*CcfIn';

OmegaC = sparse(1:ncTot, 1:ncTot, CVols, ncTot, ncTot);
Omega  = kron(eye(3), OmegaC);
OmegaCIn = sparse(1:theNumberOfElements, 1:theNumberOfElements, CVolsIn, theNumberOfElements, theNumberOfElements);
OmegaIn = kron(eye(3), OmegaCIn);

OmegaS = Af*Dnf;
OmegaSIn = Af*DnfIn;


% Construct derived matrices
% - Compact differential operators
M = (Tfo-Tfn)*Af;
G = -OmegaS\M.';
Gg = Omega\(kron(eye(3), M)*Nf.');
L = M*G;



% - Interpolators
% N.B. PiSC is not straightforward, especially in 1 or 2D
% 1/NDim is used to correct for this
% Not sure if PiSC should ever be used from differential geometry viewpoint

% Find if mesh is 2D or 3D
% Assume 2D if empty patch is found, 1D not implemented
NDim = 3;
for iBPatch=1:theNumberOfBoundaryPatches
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);
    thePhysicalPatchType = theBCInfo.type;

    if strcmp(thePhysicalPatchType,'empty')
        NDim = 2;
        break;
    end
end

PiCSM = Tfc.'/2;
PiSCM = 1/NDim * OmegaC\(PiCSM.'*OmegaS);
GamCSM = Nf*kron(eye(3), PiCSM);
GamSCM = Omega\(GamCSM.'*OmegaS);

PiCSL = Wnf*Tfo.' + Wof*Tfn.';
PiSCL = 1/NDim * OmegaC\(PiCSL.'*OmegaS);
GamCSL = Nf*kron(eye(3), PiCSL);
GamSCL = Omega\(GamCSL.'*OmegaS);

PiCSV = Wof*Tfo.' + Wnf*Tfn.';
PiSCV = 1/NDim * OmegaC\(PiCSV.'*OmegaS);
GamCSV = Nf*kron(eye(3), PiCSV);
GamSCV = Omega\(GamCSV.'*OmegaS);

% - Collocated differential operators
if strcmp(Region.foamDictionary.fvSolution.AlguFVM.interpolation, 'linear')
    PiCS = PiCSL;
elseif strcmp(Region.foamDictionary.fvSolution.AlguFVM.interpolation, 'midpoint')
    PiCS = PiCSM;
elseif strcmp(Region.foamDictionary.fvSolution.AlguFVM.interpolation, 'volumetric')
    PiCS = PiCSV;
end

PiSC = 1/NDim * PiCS.';
GamCS = Nf*kron(eye(3), PiCS);
GamSC = Omega\(GamCS.'*OmegaS);

% - Back and forth interpolators 
% N.B. not symmetric because of boundary corrections

GamCSC = GamSC*GamCS;
% Replace rows of ghost cells with 1's on diagonal
bElementIndices = theNumberOfElements + (1:theNumberOfBFaces);
bElementIndicesXYZ = reshape(double(ncTot) * [0, 1, 2] + double(bElementIndices'), 1, []);
GamCSCbElementsXYZ = speye(size(GamCSC));
GamCSC(bElementIndicesXYZ,:) = GamCSCbElementsXYZ(bElementIndicesXYZ,:);

GamSCS = GamCS*GamSC;
% Replace rows of boundary faces with 1's on diagonal
bFaceIndices = theNumberOfInteriorFaces + (1:theNumberOfBFaces);
GamSCSbFaces = speye(size(GamSCS));
GamSCS(bFaceIndices,:) = GamSCSbFaces(bFaceIndices,:);

% Collocated operators
Mc = M*GamCS;
Gc = -Omega\(GamCS.'*M.');

% Set Region variables
Region.operators.Tfo = Tfo;
Region.operators.Tfn = Tfn;
Region.operators.Tfc = Tfc;
Region.operators.Tff = Tff;

Region.operators.Sf = Sf;
Region.operators.Nf = Nf;
Region.operators.Af = Af;
Region.operators.Ccf = Ccf;
Region.operators.Cof = Cof;
Region.operators.Cnf = Cnf;
Region.operators.Dnf = Dnf;
Region.operators.Dnof = Dnof;
Region.operators.Dnnf = Dnnf;
Region.operators.Wof = Wof;
Region.operators.Wnf = Wnf;
Region.operators.CcfIn = CcfIn;
Region.operators.DnfIn = DnfIn;

Region.operators.OmegaC = OmegaC;
Region.operators.Omega  = Omega;
Region.operators.OmegaCIn = OmegaCIn;
Region.operators.OmegaIn  = OmegaIn;
Region.operators.OmegaS = OmegaS;
Region.operators.OmegaSIn = OmegaSIn;

Region.operators.M = M;
Region.operators.G = G;
Region.operators.Gg = Gg;
Region.operators.L = L;

Region.operators.PiCSM = PiCSM;
Region.operators.PiSCM = PiSCM;
Region.operators.GamCSM = GamCSM;
Region.operators.GamSCM = GamSCM;
Region.operators.PiCSL = PiCSL;
Region.operators.PiSCL = PiSCL;
Region.operators.GamCSL = GamCSL;
Region.operators.GamSCL = GamSCL;
Region.operators.PiCSV = PiCSV;
Region.operators.PiSCV = PiSCV;
Region.operators.GamCSV = GamCSV;
Region.operators.GamSCV = GamSCV;
Region.operators.PiCS = PiCS;
Region.operators.PiSC = PiSC;
Region.operators.GamCS = GamCS;
Region.operators.GamSC = GamSC;
Region.operators.GamCSC = GamCSC;
Region.operators.GamSCS = GamSCS;

Region.operators.Mc= Mc;
Region.operators.Gc = Gc;
