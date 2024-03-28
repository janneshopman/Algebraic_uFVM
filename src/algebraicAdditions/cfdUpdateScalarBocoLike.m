%AlgebraicAdjustment
function theField = cfdUpdateScalarBocoLike(theField, theLikeFieldName)

theNumberOfBPatches = cfdGetNumberOfBPatches;

for iBPatch=1:theNumberOfBPatches
    theBCInfo = cfdGetBoundaryPatchRef(iBPatch);
    thePhysicalPatchType = theBCInfo.type;
    theBCType = cfdBcForBoundaryPatch(theLikeFieldName, iBPatch);
    %
    % WALL
    %
    if strcmp(thePhysicalPatchType,'wall')
        if strcmp(theBCType,'fixedValue') || strcmp(theBCType,'calculated')
            theField = cfdUpdateFixedValue(iBPatch, theField, theLikeFieldName);
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'slip')
            theField = cfdUpdateZeroGradient(iBPatch, theField);
        elseif strcmp(theBCType,'noSlip')
            theField = cfdUpdateNoSlip(iBPatch, theField);
        else
            error([theBCType ' bc not defined']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalPatchType,'inlet')
        if strcmp(theBCType,'fixedValue')
            theField = cfdUpdateFixedValue(iBPatch, theField, theLikeFieldName);
        elseif strcmp(theBCType,'zeroGradient')
            theField = cfdUpdateZeroGradient(iBPatch, theField);
        else
            error('Inlet bc not defined');
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalPatchType,'outlet')
        if strcmp(theBCType,'fixedValue')
            theField = cfdUpdateFixedValue(iBPatch, theField, theLikeFieldName);            
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'outlet')
            theField = cfdUpdateZeroGradient(iBPatch, theField);            
        else
            error([theBCType 'Outlet bc not defined']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalPatchType,'symmetry')
        theField = cfdUpdateZeroGradient(iBPatch, theField);
        %
        % EMPTY
        %        
    elseif strcmp(thePhysicalPatchType,'empty')
        % Do nothing
        %
        % CYCLIC
        %        
    elseif strcmp(thePhysicalPatchType,'cyclic')      
        theField = updateCyclic(iBPatch, theField);     
    else
        error([thePhysicalPatchType '<<<< Physical Condition bc not defined']);
        
    end
                
end

end



%===================================================
% Fixed Value
%===================================================
function theField = cfdUpdateFixedValue(iBPatch, theField, theLikeFieldName)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Apply cfdBoundary condition
theBCValue = cfdValueForBoundaryPatch(theLikeFieldName, iBPatch);
if size(theBCValue, 1)>1
    theField(iBElements) = 2.0*theBCValue - theField(owners_b);
else
    theField(iBElements) = 2.0*theBCValue*ones(length(iBElements),1) - theField(owners_b);
end
end


%===================================================
% No Slip
%===================================================
function theField = cfdUpdateNoSlip(iBPatch, theField)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Apply cfdBoundary condition
theField(iBElements) = -theField(owners_b);
end


%===================================================
% Zero Gradient
%===================================================
function theField = cfdUpdateZeroGradient(iBPatch, theField)

% Get info
iBElements = cfdGetBoundaryElementsSubArrayForBoundaryPatch(iBPatch);
owners_b = cfdGetOwnersSubArrayForBoundaryPatch(iBPatch);

% Copy owner values
theField(iBElements) = theField(owners_b);
end

%===================================================
% Cyclic
%===================================================
function theField = updateCyclic(iBPatch, theField)

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

theField(iBElements) = theField(owners_Nb);
end
