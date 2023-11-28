function sampleIndices = cfdLineSampleIndices(mesh, beginCoordinate, endCoordinate, numberOfPoints)
    linePoints = zeros(numberOfPoints, 3);
    sampleIndices = zeros(numberOfPoints, 1);
    
    for i=1:3
        linePoints(:, i) = linspace(beginCoordinate(i), endCoordinate(i), numberOfPoints);
    end
    
    theNumberOfElements = mesh.numberOfElements;
    theNumberOfFaces = mesh.numberOfFaces;
    theNumberOfInteriorFaces = mesh.numberOfInteriorFaces;
    
    owners = mesh.owners;
    neighbours = mesh.neighbours;
    
    faceCentroids = mesh.faceCentroids;
    faceSfs =  mesh.faceSf;
    
    FO = sparse(owners, 1:double(theNumberOfFaces), 1, theNumberOfElements, theNumberOfFaces);
    FN = sparse(neighbours, 1:double(theNumberOfInteriorFaces), 1, theNumberOfElements, theNumberOfFaces);
    FC = FO-FN;
    
    for iLP = 1:size(linePoints,1)
        LP = linePoints(iLP,:);
    
        for iEl = 1:theNumberOfElements
            elementFaces = find(FC(iEl, :));

            allInside = true;
    
            for iElFa = 1:length(elementFaces)
                iFa = elementFaces(iElFa);
                
                nfOut = FC(iEl,iFa)*faceSfs(iFa,:)/vecnorm(faceSfs(iFa,:));

                LPFaceVec = faceCentroids(iFa,:) - LP;

                distanceToOutside = LPFaceVec*nfOut';

                if distanceToOutside < -1E-12
                    allInside = false;
                    break;
                end
            end

            if allInside
                sampleIndices(iLP) = iEl;
                break;
            end
        end
    end
end
