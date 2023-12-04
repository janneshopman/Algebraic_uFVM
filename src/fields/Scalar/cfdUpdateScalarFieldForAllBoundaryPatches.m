%AlgebraicAdjustment
function cfdUpdateScalarFieldForAllBoundaryPatches(theFieldName)

theNumberOfBPatches = cfdGetNumberOfBPatches;

for iBPatch=1:theNumberOfBPatches
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);
    thePhysicalPatchType = theBCInfo.type;
    theBCType = cfdBcForBoundaryPatch(theFieldName, iBPatch);
    %
    % WALL
    %
    if strcmp(thePhysicalPatchType,'wall')
        if strcmp(theBCType,'fixedValue') || strcmp(theBCType,'calculated')
            cfdUpdateFixedValue(iBPatch,theFieldName);
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'slip')
            cfdUpdateZeroGradient(iBPatch,theFieldName);
        elseif strcmp(theBCType,'noSlip')
            cfdUpdateNoSlip(iBPatch,theFieldName);
        else
            error([theBCType ' bc not defined']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalPatchType,'inlet')
        if strcmp(theBCType,'fixedValue')
            cfdUpdateFixedValue(iBPatch,theFieldName);
        elseif strcmp(theBCType,'zeroGradient')
            cfdUpdateZeroGradient(iBPatch,theFieldName);
        else
            error('Inlet bc not defined');
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalPatchType,'outlet')
        if strcmp(theBCType,'fixedValue')
            cfdUpdateFixedValue(iBPatch,theFieldName);            
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'outlet')
            cfdUpdateZeroGradient(iBPatch,theFieldName);            
        else
            error([theBCType 'Outlet bc not defined']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalPatchType,'symmetry')
        cfdUpdateZeroGradient(iBPatch,theFieldName);
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
function cfdUpdateFixedValue(iBPatch, theFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Get field
theScalarField = cfdGetMeshField(theFieldName);

% Apply cfdBoundary condition
theBCValue = cfdValueForBoundaryPatch(theFieldName, iBPatch);
if size(theBCValue, 1)>1
    theScalarField.phi(iBElements) = 2.0*theBCValue - theScalarField.phi(owners_b);
else
    theScalarField.phi(iBElements) = 2.0*theBCValue*ones(length(iBElements),1) - theScalarField.phi(owners_b);
end

% Store
cfdSetMeshField(theScalarField);

end


%===================================================
% No Slip
%===================================================
function cfdUpdateNoSlip(iBPatch, theFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Get field
theScalarField = cfdGetMeshField(theFieldName);

% Apply cfdBoundary condition
theBCValue = cfdValueForBoundaryPatch(theFieldName, iBPatch);
theScalarField.phi(iBElements) = -theScalarField.phi(owners_b);

% Store
cfdSetMeshField(theScalarField);

end


%===================================================
% Zero Gradient
%===================================================
function cfdUpdateZeroGradient(iBPatch, theFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Get field
theScalarField = cfdGetMeshField(theFieldName);

% Copy owner values
theScalarField.phi(iBElements) = theScalarField.phi(owners_b);

% Store
cfdSetMeshField(theScalarField);

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

theScalarField = cfdGetMeshField(theFieldName);

theScalarField.phi(iBElements) = theScalarField.phi(owners_Nb);

cfdSetMeshField(theScalarField);

end
