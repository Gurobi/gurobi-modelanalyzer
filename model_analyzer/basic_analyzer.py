def process_basic_information(model, data):
    data["MinCoeff"] = model.MinCoeff
    data["MaxCoeff"] = model.MaxCoeff
    data["MinBound"] = model.MinBound
    data["MaxBound"] = model.MaxBound
    data["MinRHS"] = model.MinRHS
    data["MaxRHS"] = model.MaxRHS

    data["MaxQCCoeff"] = model.MaxQCCoeff
    data["MinQCCoeff"] = model.MinQCCoeff
    data["MaxQCLCoeff"] = model.MaxQCLCoeff
    data["MinQCLCoeff"] = model.MinQCLCoeff
    data["MaxQCRHS"] = model.MaxQCRHS
    data["MinQCRHS"] = model.MinQCRHS

    data["numNZs"] = model.NumNZs
    data["NumQNZs"] = model.NumQNZs
    data["NumQCNZs"] = model.NumQCNZs

    data["NumConstrs"] = model.NumConstrs
    data["NumQConstrs"] = model.NumQConstrs
    data["NumSOS"] = model.NumSOS
    data["NumGenConstrs"] = model.NumGenConstrs

    data["MinObjCoeff"] = model.MinObjCoeff
    data["MaxObjCoeff"] = model.MaxObjCoeff
    data["MaxQObjCoeff"] = model.MaxQObjCoeff
    data["MinQObjCoeff"] = model.MinQObjCoeff

    data["NumVars"] = model.NumVars
    data["NumIntVars"] = model.NumIntVars
    data["NumBinVars"] = model.NumBinVars
    data["NumPWLObjVars"] = model.NumPWLObjVars

    data["ModelName"] = model.ModelName
