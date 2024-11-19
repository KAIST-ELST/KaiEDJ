using Optimization
using OptimizationOptimJL

export BathDiscHybPH

function BathDiscHybPH( BParam, SParam ; outputlevel=0 )
    @time cost = GetCostFromFlatPH( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam)
    sol = solve(prob, NelderMead())
    if outputlevel > 0 
        @show sol.original
    end

    BParamNew   = [ sol... ]
    return BParamNew
end
