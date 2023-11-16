function allNeighbours = cfdGetAllNeighbours
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

global Region;

neighbours = Region.mesh.neighbours;

theNumberOfElements = cfdGetNumberOfElements;

theNumberOfBFaces = cfdGetNumberOfBFaces;

bNeighbours = theNumberOfElements+1:theNumberOfElements+theNumberOfBFaces;

allNeighbours = double([neighbours; bNeighbours']);