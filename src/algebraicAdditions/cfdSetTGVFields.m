function cfdSetTGV

global Region;

runTime = cfdGetCurrentTime;

U = cfdGetMeshField('U');
p = cfdGetMeshField('p');
Uf = cfdGetMeshField('Uf');

theNumberOfElements = cfdGetNumberOfElements;
elementCentroids = Region.mesh.elementCentroids;

for iElement=1:theNumberOfElements
    x = elementCentroids(iElement, 1);
    y = elementCentroids(iElement, 2);

    U.phi(iElement, 1) =  sin(x)*cos(y);
    U.phi(iElement, 2) = -cos(x)*sin(y);

    p.phi(iElement)    =  0.25*(cos(2*x) + cos(2*y));
end

theNumberOfFaces = cfdGetNumberOfFaces;
faceCentroids = cfdGetFaceCentroidsSubArrayForFaces;

faceSfs = cfdGetFaceSf;

for iFace=1:theNumberOfFaces
    x = faceCentroids(iFace, 1);
    y = faceCentroids(iFace, 2);

    Us =  [sin(x)*cos(y), -cos(x)*sin(y), 0];

    nf = faceSfs(iFace,:)/norm(faceSfs(iFace,:));
    
    Uf.phi(iFace) = dot(nf, Us);
end

cfdSetMeshField(U);
cfdSetMeshField(p);
cfdSetMeshField(Uf);

cfdUpdateVectorFieldForAllBoundaryPatches('U');
cfdUpdateScalarFieldForAllBoundaryPatches('p');

cfdWriteOpenFoamField('U', runTime);
cfdWriteOpenFoamField('p', runTime);
cfdWriteOpenFoamField('Uf', runTime);