function allNeighbours = cfdGetAllNeighbours

global Region;

neighbours = Region.mesh.neighbours;

theNumberOfElements = cfdGetNumberOfElements;

theNumberOfBFaces = cfdGetNumberOfBFaces;

bNeighbours = theNumberOfElements+1:theNumberOfElements+theNumberOfBFaces;

allNeighbours = double([neighbours; bNeighbours']);
