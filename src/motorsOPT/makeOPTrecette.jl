function TaylorOptions(itplParams,supplementaryOrder)
    options=(YorderBspace=itplParams.YorderBspace,YorderBtime=itplParams.YorderBtime,supplementaryOrder=supplementaryOrder,pointsμInSpace=itplParams.ptsSpace,pointsμInTime=itplParams.ptsTime,offsetμInΔyInSpace=itplParams.offsetSpace,offsetμInΔyInTime=itplParams.offsetTime)
    return options
end

function makeOPTsemiSymbolic(params::Dict)
    @unpack famousEquationType, Δ, orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = params


    # construction of NamedTuples
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,pointsInSpace=pointsInSpace,pointsInTime=pointsInTime)

    # here we can compute the different interpolated Taylor expansion options
    TaylorOptionsμ=TaylorOptions(fieldItpl,supplementaryOrder)
    TaylorOptionsμᶜ=TaylorOptions(materItpl,supplementaryOrder)

    equationCharacteristics=famousEquations(famousEquationType)
    numbersOfTheSystem=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    dependencies,ordersForSplines,configsTaylor=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptions)
    bigα,varM=bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplines)

end