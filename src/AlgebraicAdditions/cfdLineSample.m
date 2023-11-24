Region = OF;
field = OF.fluid.U.phi;
% beginCoordinate = [0.475,0.025,0.5];
% endCoordinate = [0.475,0.975,0.5];
beginCoordinate = [0.5,0.0,0.0];
endCoordinate = [0.5,1.0,0.0];
numberOfPoints = 21;

sampleIndices = cfdLineSampleIndices(Region.mesh, beginCoordinate, endCoordinate, numberOfPoints);

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
    
    elementFaces = mesh.elementFaces;
    faceCentroids = mesh.faceCentroids;
    faceSfs =  mesh.faceSf;
    elementCentroids = mesh.elementCentroids;
    
    FO = sparse(owners, 1:double(theNumberOfFaces), 1, theNumberOfElements, theNumberOfFaces);
    FN = sparse(neighbours, 1:double(theNumberOfInteriorFaces), 1, theNumberOfElements, theNumberOfFaces);
    FC = FO-FN;
    
    for iLP = 1:length(linePoints)
        LP = linePoints(iLP,:);
    
        for iEl = 1:theNumberOfElements
            elementFaces = find(FC(iEl, :));

            allPositiveDots = true;
    
            for iElFa = 1:length(elementFaces)
                iFa = elementFaces(iElFa);
                
                nf = faceSfs(iFa,:)/vecnorm(faceSfs(iFa,:));

                LPFaceVec = faceCentroids(iFa,:) - LP;
    
                LPCentroidVec = elementCentroids(iEl, :) - LP;

                distance = LPFaceVec*FC(iEl,iFa)*nf';

                if distance < 0
                    allPositiveDots = false;
                    break;
                elseif distance/vecnorm(LPCentroidVec) < 1E-9
                    %found point on face, nudge the line a bit to center
                    fprintf('NUDGED');
                    linePoints(1:end-1,:) = linePoints(1:end-1,:) + LPCentroidVec*1E-2;
                    LP = linePoints(iLP,:);
                end
            end

            if allPositiveDots
                sampleIndices(iLP) = iEl;
                break;
            end
        end
    end
end
