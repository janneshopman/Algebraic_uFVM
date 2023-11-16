%AlgebraicAdjustment
function cfdUpdateVectorFieldForAllBoundaryPatches(theFieldName)

theNumberOfBPatches = cfdGetNumberOfBPatches;

for iBPatch=1:theNumberOfBPatches
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);
    thePhysicalPatchType = theBCInfo.type;
    theBCType = cfdBcForBoundaryPatch(theFieldName, iBPatch);
    %
    % WALL
    %
    if strcmp(thePhysicalPatchType,'wall')
        if strcmp(theBCType,'fixedValue')
            updateFixedValue(iBPatch,theFieldName);
        elseif strcmp(theBCType,'zeroGradient')
            updateZeroGradient(iBPatch,theFieldName);
        elseif strcmp(theBCType,'slip')
            updateSymmetry(iBPatch,theFieldName);
        elseif strcmp(theBCType,'noSlip')
            updateNoSlip(iBPatch,theFieldName);
        else
            error([theBCType ' bc not defined']);
        end
        %
        % INLET
        %
    elseif (strcmp(thePhysicalPatchType,'inlet'))
        if (strcmp(theBCType,'fixedValue'))
            updateFixedValue(iBPatch,theFieldName);
        elseif strcmp(theBCType,'zeroGradient')
            updateZeroGradient(iBPatch,theFieldName);
        else
            error('Inlet bc not defined');
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalPatchType,'outlet')
        if strcmp(theBCType,'fixedValue')
            updateFixedValue(iBPatch,theFieldName);            
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'outlet')
            updateZeroGradient(iBPatch,theFieldName);            
        else
            error([theBCType 'Outlet bc not defined']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalPatchType,'symmetry')
        updateSymmetry(iBPatch,theFieldName);
        %
        % EMPTY
        %        
    elseif strcmp(thePhysicalPatchType,'empty')
        % Do nothing
        %
        % CYCLIC
        %        
    elseif strcmp(thePhysicalPatchType,'cyclic')
        updateCyclic(iBPatch, theFieldName);     
    else
        error([thePhysicalPatchType '<<<< Physical Condition bc not defined']);
        
    end
                
end

end



%===================================================
% Fixed Value
%===================================================
function updateFixedValue(iBPatch, theFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Get field
theVectorField = cfdGetMeshField(theFieldName);

% Apply cfdBoundary condition
theBCValue = cfdValueForBoundaryPatch(theFieldName, iBPatch);
if size(theBCValue,1)>1
    theVectorField.phi(iBElements,:) = 2.0*theBCValue - theVectorField.phi(owners_b,:);
else
    theVectorField.phi(iBElements,:) = 2.0*repmat(theBCValue,length(iBElements),1) - theVectorField.phi(owners_b,:);
end

% Store
cfdSetMeshField(theVectorField);

end


%===================================================
% Zero Gradient
%===================================================
function updateZeroGradient(iBPatch, theFieldName)

% Get info
theVectorField = cfdGetMeshField(theFieldName);

iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Apply cfdBoundary condition
theVectorField.phi(iBElements,:) = theVectorField.phi(owners_b,:);

% Store
cfdSetMeshField(theVectorField);

end


%===================================================
% No Slip
%===================================================
function updateNoSlip(iBPatch, theFieldName)

% Get info
theVectorField = cfdGetMeshField(theFieldName);

iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Apply cfdBoundary condition
theVectorField.phi(iBElements,:) = -theVectorField.phi(owners_b,:);

% Store
cfdSetMeshField(theVectorField);

end


%===================================================
% Symmetry
%===================================================
function updateSymmetry(iBPatch, theFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Get field
theVectorField = cfdGetMeshField(theFieldName);

% Update bounadry condition
Sb = cfdGetFaceSfSubArrayForBoundaryPatch(iBPatch);
normSb = cfdMag(Sb);
n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

U_normal_cfdMag = dot(theVectorField.phi(owners_b,:)',n')';
U_normal = [U_normal_cfdMag .* n(:,1),U_normal_cfdMag .* n(:,2),U_normal_cfdMag .* n(:,3)];

theVectorField.phi(iBElements,:) = theVectorField.phi(owners_b,:) - U_normal;

% Store
cfdSetMeshField(theVectorField);

end

%===================================================
% Cyclic
%===================================================
function updateCyclic(iBPatch, theFieldName)

theBCInfo = cfdGetBoundaryPatchRef(iBPatch);

iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);  

theNumberOfBoundaryPatches = cfdGetNumberOfBPatches;

for iNBPatch=1:theNumberOfBoundaryPatches
    theNBCInfo = cfdGetBoundaryPatchRef(iNBPatch);
    if strcmp(theNBCInfo.name, theBCInfo.neighbourPatch)                  
        owners_Nb = cfdGetOwnersSubArrayForBoundaryPatch(iNBPatch);  

        break
    end
end   

theVectorField = cfdGetMeshField(theFieldName);

theVectorField.phi(iBElements,:) = theVectorField.phi(owners_Nb,:);

cfdSetMeshField(theVectorField);

end
