function cfdPerturbateField

global Region;

if Region.foamDictionary.fvSolution.AlguFVM.perturbate
    U = cfdGetField('U');
    innerElem = Region.mesh.innerElemVec;
    Ufield = U(innerElem);
    randVec = rand(size(Ufield));
    Ufield = Ufield.*(1+0.1*randVec);
    U(innerElem) = Ufield;
    cfdSetField(U,'U');
end