
function quasiNumericalOperatorConstruction(operators,modelName,models,forceModels,famousEquationType,modelPoints,IneedExternalSources;maskedRegionForFieldInSpace = nothing,maskedRegionForSourcesInSpace=nothing,iExperiment=nothing)

    NpointsUsed=iExperiment.iPointsUsed
    # this is a big wrapper that reads the semi symbolic expressions to give a set of numerical operators (with symbolic expression in time)
    # which will call wrappers of onstructingNumericalDiscretisedEquations(Masked)

    operators=operators["operators"]
    operatorPDE,operatorForce,eqInfo = operators
    exprs,fields,vars,extexprs,extfields,extvars,coordinates=eqInfo

    AjiννᶜU,utilities=operatorPDE
    if IneedExternalSources 
        Γg,utilitiesForce=operatorForce
    end

    lhsConfigurations = @strdict semiSymbolicOpt=AjiννᶜU coordinates modelName models fields vars famousEquationType modelPoints utilities maskedRegion=maskedRegionForFieldInSpace NpointsUsed

    #numOperators,file = @produce_or_load(constructingNumericalDiscretisedEquations,lhsConfigurations,datadir("numOperators",savename(lhsConfigurations));filename = config -> savename(lhsConfigurations; ignores=["vars", "fields"]))
    numOperators = myProduceOrLoad(constructingNumericalDiscretisedEquations,lhsConfigurations,"numOperators","lhs")


    # left-hand side, which is far more recyclable than r.h.s.
    
    costfunctionsLHS,fieldLHS,_=numOperators["numOperators"]
    

    costfunctionsRHS = similar(costfunctionsLHS)
    costfunctionsRHS .= 0.
    fieldRHS = similar(fieldLHS)
    fieldRHS .= 0.
    
    champsLimité = nothing

    if IneedExternalSources 
        rhsConfigurations = @strdict semiSymbolicOpt=Γg coordinates modelName models=forceModels fields=extfields vars=extvars famousEquationType modelPoints utilities=utilitiesForce maskedRegion=maskedRegionForSourcesInSpace NpointsUsed
        numOperators,file=produce_or_load(constructingNumericalDiscretisedEquations,rhsConfigurations,datadir("numOperators",savename(rhsConfigurations));filename = config -> savename("source",rhsConfigurations; ignores=["vars", "fields"]))

        numOperators = myProduceOrLoad(constructingNumericalDiscretisedEquations,rhsConfigurations,"numOperators","rhs")
       
        costfunctionsRHS,fieldRHS,champsLimité=numOperators["numOperators"]

    end
  
    #@show size(costfunctionsLHS),size(costfunctionsRHS)
    #@show costfunctionsRHS[1,430],costfunctionsRHS[1,431],costfunctionsRHS[1,434]
    #costfunctions=0.#costfunctionsLHS[1,1]-costfunctionsLHS[1,1]
    costfunctions = costfunctionsLHS .- costfunctionsRHS
    #numOperators=(costfunctions=costfunctions)
    return costfunctions,fieldLHS,fieldRHS,champsLimité
end
