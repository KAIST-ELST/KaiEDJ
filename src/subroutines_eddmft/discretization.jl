using Optimization
using OptimizationOptimJL

export BathParamFlatten
export BathParamReshape
export BathParamReshapePH
export BathParamReshapePHev
export BathDiscHybSpinPH
export BathDiscHybSpinPHev
export BathDiscHybSpin
export BathDiscHybPH
export BathDiscHybPHev
export BathDiscHyb
export BathDiscHybChoose
export BathSymPH
export ShowBathParam
export ReadBathParam!
export WriteBathParam
export InitBathParam!

export WriteBathParamInfo
export WriteBathEnergy
export WriteBathHyb
export WriteBathAll


export BathParams

mutable struct BathParams
    ebathl::Vector{Float64}
    Vil::Matrix{Float64}
    BathParams( norb::Int, nbath::Int )                         = new( zeros(Float64,nbath), zeros(Float64, norb, nbath) )
    BathParams( bp::BathParams )                                = new( bp.ebathl, bp.Vil )
    BathParams( el::Vector{Float64}, vil::Matrix{Float64} )     = new( el, vil )
end


function BathParamFlatten( ebathl::Vector{Float64}, Vil::Matrix{Float64} )
    return [ ebathl... (Vil...)... ]   # Note that it expands in the column-major way.
end

function BathParamReshape( BParam, DimBath )
    ebathl  = BParam[1:DimBath]
    Vil     = reshape( BParam[DimBath+1:end], :, DimBath )
    return ebathl, Vil
end

function BathParamReshapePH( BParamPH, DimBath )
    DimBathHalf = div( DimBath , 2 )
    ebathl      = BParamPH[1:DimBathHalf]
    if DimBath%2 < 1
        ebathl      = [ ebathl... reverse(-ebathl)... ]
    else
        ebathl      = [ ebathl... [0.0]... reverse(-ebathl)... ]
    end
    Vil         = reshape( BParamPH[DimBathHalf+1:end], :, DimBath )
    return ebathl, Vil
end

function BathParamReshapePHev( BParamPH, DimOrb, DimBath )
    DimBathHalf     = div( DimBath , 2 )
    IndBathEnergy   = collect(1:DimBathHalf)
    IndBathHyb      = []
    ebathl          = BParamPH[IndBathEnergy]
    if DimBath%2 < 1
        ebathl      = [ ebathl... reverse(-ebathl)... ]
        IndBathHyb  = [ IndBathEnergy... reverse(IndBathEnergy)... ]
    else
        ebathl      = [ ebathl... [0.0]... reverse(-ebathl)... ]
        IndBathHyb  = [ IndBathEnergy... DimBathHalf+1 reverse(IndBathEnergy)... ]
    end
    Vil     = BParamPH[DimBathHalf+1:end]
    Vil     = reshape( Vil, DimOrb, : )
    Vil     = Vil[:,IndBathHyb]
    return ebathl, Vil
end

@inline function BathDiscHybChoose( ebathl_init, Vil_init, Hybiw, FreqGrid, BathOpt, SpinSep )
    if BathOpt
        if SpinSep
            ebathlnew, Vilnew, nIterOpt, costFinal   = BathDiscHybSpin( ebathl_init, Vil_init, Hybiw, FreqGrid )
            return ebathlnew, Vilnew, nIterOpt, costFinal
        else
            ebathlnew, Vilnew, nIterOpt, costFinal   = BathDiscHybAll( ebathl_init, Vil_init, Hybiw, FreqGrid )
            return ebathlnew, Vilnew, nIterOpt, costFinal
        end
    else
        return ebathl_init, Vil_init, -1, -1
    end
end

@inline function BathDiscHyb( BParam, SParam ; outputlevel=0 )
    return BathDiscHybBFGS(BParam,SParam;outputlevel=outputlevel)
end

function BathDiscHybFunct( BParam, SParam, minimizer ; outputlevel=0 )
    cost = GetCostFromFlat( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    res = optimize( x -> GetCostFromFlat(x,SParam), BParam, minimizer, Optim.Options(iterations=2500, x_tol=1e-14, x_reltol=1e-14, f_abstol=1e-10, f_tol=1e-10, time_limit=90000))
    if outputlevel > 0 
        @show res
    end

    niter       = 0
    cost_final  = -1
    BParamNew   = [ res.minimizer... ]
    try 
        niter       = res.iterations
        cost_final  = res.minimum
    catch
    end
    return BParamNew, niter, cost_final
end

function BathDiscHybBFGS( BParam, SParam ; outputlevel=0 )
    cost = GetCostFromFlat( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    println( "Optimizer : BFGS()" )
    t_opt = @elapsed @time res = optimize( x -> GetCostFromFlat(x,SParam), BParam, BFGS(), Optim.Options(iterations=25000, x_tol=1e-14, x_reltol=1e-14, f_abstol=1e-10, f_tol=1e-10, time_limit=90000))
    if outputlevel > 0 
        @show res
    end
    open( "data/time_optimizer.txt", "a" ) do io
        writedlm( io, t_opt ) 
    end

    niter       = 0
    cost_final  = -1
    BParamNew   = [ res.minimizer... ]
    try 
        niter       = res.iterations
        cost_final  = res.minimum
    catch
    end
    return BParamNew, niter, cost_final
end

function BathDiscHybNM( BParam, SParam ; outputlevel=0 )
    cost = GetCostFromFlat( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    # lbound  = [ -20.0 for i in BParam ]
    # ubound  = [  20.0 for i in BParam ]
    prob = OptimizationProblem(GetCostFromFlat, BParam, SParam) #, lb=lbound, ub=ubound)
    sol = solve(prob, NelderMead() ; maxiters=2500 )
    if outputlevel > 0 
        @show sol.original
    end

    niter       = 0
    cost_final  = -1
    BParamNew   = [ sol... ]
    try 
        niter       = sol.original.iterations
        cost_final  = sol.objective
    catch
    end
    return BParamNew, niter, cost_final
end

function BathDiscHybPH( BParam, SParam ; outputlevel=0 )
    cost = GetCostFromFlatPH( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    # # lbound  = [ -20.0 for i in BParam ]
    # # ubound  = [  20.0 for i in BParam ]
    # prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam) # , lb=lbound, ub=ubound)
    # sol = solve(prob, NelderMead() ; maxiters=2500 )
    # if outputlevel > 0 
    #     @show sol.original
    # end
    # BParamNew   = [ sol... ]

    res = optimize( x -> GetCostFromFlatPH(x,SParam), BParam, LBFGS(), Optim.Options(iterations=2500, x_tol=1e-14, x_reltol=1e-14, f_abstol=1e-10, f_tol=1e-10, time_limit=90000))
    if outputlevel > 0 
        @show res
    end
    BParamNew   = [ res.minimizer... ]

    return BParamNew
end

function BathDiscHybPHev( BParam, SParam ; outputlevel=0 )
    cost = GetCostFromFlatPHev( BParam, SParam )
    if outputlevel > 0 
        println( "Initial cost : $(cost) " ) 
    end

    # # lbound  = [ -20.0 for i in BParam ]
    # # ubound  = [  20.0 for i in BParam ]
    # prob = OptimizationProblem(GetCostFromFlatPHev, BParam, SParam) # , lb=lbound, ub=ubound)
    # sol = solve(prob, NelderMead() ; maxiters=2500 )
    # if outputlevel > 0 
    #     @show sol.original
    # end
    # BParamNew   = [ sol... ]

    res = optimize( x -> GetCostFromFlatPH(x,SParam), BParam, LBFGS(), Optim.Options(iterations=2500, x_tol=1e-14, x_reltol=1e-14, f_abstol=1e-10, f_tol=1e-10, time_limit=90000))
    if outputlevel > 0 
        @show res
    end
    BParamNew   = [ res.minimizer... ]

    return BParamNew
end

function BathDiscHybSpinPH( ebathl, Vil, hybiw_grid, FreqGrid )
    nspinorb, nbath = size( Vil )
    norb            = div( nspinorb, 2 )
    nbathHalf       = div( nbath, 2 )
    IndOrbUp, IndOrbDn = GetOrbUpDn( nspinorb )
    IndBathUp, IndBathDn = GetOrbUpDn( nbath )

    ebathlUp    = ebathl[1:2:end]
    ebathlDn    = ebathl[2:2:end]
    VilUp   = Vil[ IndOrbUp, IndBathUp ]
    VilDn   = Vil[ IndOrbDn, IndBathDn ]

    BParamUpPH  = BathParamFlatten( ebathlUp[1:div(end,2)], VilUp )
    BParamDnPH  = BathParamFlatten( ebathlDn[1:div(end,2)], VilDn )

    NFreq   = length( FreqGrid )
    Hybiwup_init    = [ hybiw_grid[ifreq][IndOrbUp,IndOrbUp] for ifreq in 1:NFreq ]
    Hybiwdn_init    = [ hybiw_grid[ifreq][IndOrbDn,IndOrbDn] for ifreq in 1:NFreq ]

    SParamUp  = ( 
                    Hybiwup_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    SParamDn  = ( 
                    Hybiwdn_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    BParamUpPHNew   = BathDiscHybPH( BParamUpPH, SParamUp ; outputlevel=1)
    BParamDnPHNew   = BathDiscHybPH( BParamDnPH, SParamDn ; outputlevel=1)
    enewup, vnewup  = BathParamReshapePH( BParamUpPHNew, nbathHalf )
    enewdn, vnewdn  = BathParamReshapePH( BParamDnPHNew, nbathHalf )
    println("Parameters (inital) :")
    @show BParamUpPH
    @show BParamDnPH
    println("Parameters (optimized) :")
    @show BParamUpPHNew
    @show BParamDnPHNew
    println("")

    ebathlnew       = zero( ebathl )
    Vilnew          = zero( Vil )
    ebathlnew[1:2:end]              = deepcopy( enewup )
    ebathlnew[2:2:end]              = deepcopy( enewdn )
    Vilnew[ IndOrbUp, IndBathUp ]   = deepcopy( vnewup )
    Vilnew[ IndOrbDn, IndBathDn ]   = deepcopy( vnewdn )
    return ebathlnew, Vilnew
end

function BathDiscHybSpinPHev( ebathl, Vil, hybiw_grid, FreqGrid )
    nspinorb, nbath = size( Vil )
    norb            = div( nspinorb, 2 )
    nbathHalf       = div( nbath, 2 )
    IndOrbUp, IndOrbDn = GetOrbUpDn( nspinorb )
    IndBathUp, IndBathDn = GetOrbUpDn( nbath )

    ebathlUp    = ebathl[1:2:end]
    ebathlDn    = ebathl[2:2:end]
    VilUp   = Vil[ IndOrbUp, IndBathUp ]
    VilDn   = Vil[ IndOrbDn, IndBathDn ]

    ebathlUpPH  = ebathlUp[1:div(end,2)]
    ebathlDnPH  = ebathlDn[1:div(end,2)]
    VilUpPH     = VilUp[:,1:div(end+1,2)]
    VilDnPH     = VilDn[:,1:div(end+1,2)]

    BParamUpPH  = BathParamFlatten( ebathlUpPH, VilUpPH )
    BParamDnPH  = BathParamFlatten( ebathlDnPH, VilDnPH )

    NFreq   = length( FreqGrid )
    Hybiwup_init    = [ hybiw_grid[ifreq][IndOrbUp,IndOrbUp] for ifreq in 1:NFreq ]
    Hybiwdn_init    = [ hybiw_grid[ifreq][IndOrbDn,IndOrbDn] for ifreq in 1:NFreq ]

    SParamUp  = ( 
                    Hybiwup_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    SParamDn  = ( 
                    Hybiwdn_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    BParamUpPHNew   = BathDiscHybPHev( BParamUpPH, SParamUp ; outputlevel=1)
    BParamDnPHNew   = BathDiscHybPHev( BParamDnPH, SParamDn ; outputlevel=1)
    enewup, vnewup  = BathParamReshapePHev( BParamUpPHNew, norb, nbathHalf )
    enewdn, vnewdn  = BathParamReshapePHev( BParamDnPHNew, norb, nbathHalf )
    println("Parameters (inital) :")
    @show BParamUpPH
    @show BParamDnPH
    println("Parameters (optimized) :")
    @show BParamUpPHNew
    @show BParamDnPHNew
    println("")

    ebathlnew       = zero( ebathl )
    Vilnew          = zero( Vil )
    ebathlnew[1:2:end]              = deepcopy( enewup )
    ebathlnew[2:2:end]              = deepcopy( enewdn )
    Vilnew[ IndOrbUp, IndBathUp ]   = deepcopy( vnewup )
    Vilnew[ IndOrbDn, IndBathDn ]   = deepcopy( vnewdn )
    return ebathlnew, Vilnew
end

function BathDiscHybSpin( ebathl, Vil, hybiw_grid, FreqGrid )
    nspinorb, nbath = size( Vil )
    norb            = div( nspinorb, 2 )
    nbathHalf       = div( nbath, 2 )
    IndOrbUp, IndOrbDn = GetOrbUpDn( nspinorb )
    IndBathUp, IndBathDn = GetOrbUpDn( nbath )

    ebathlUp    = ebathl[1:2:end]
    ebathlDn    = ebathl[2:2:end]
    VilUp   = Vil[ IndOrbUp, IndBathUp ]
    VilDn   = Vil[ IndOrbDn, IndBathDn ]

    BParamUp  = BathParamFlatten( ebathlUp, VilUp )
    BParamDn  = BathParamFlatten( ebathlDn, VilDn )

    NFreq   = length( FreqGrid )
    Hybiwup_init    = [ hybiw_grid[ifreq][IndOrbUp,IndOrbUp] for ifreq in 1:NFreq ]
    Hybiwdn_init    = [ hybiw_grid[ifreq][IndOrbDn,IndOrbDn] for ifreq in 1:NFreq ]

    SParamUp  = ( 
                    Hybiwup_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    SParamDn  = ( 
                    Hybiwdn_init, 
                    FreqGrid, 
                    norb,
                    nbathHalf
                    )
    BParamUpNew, niterUp, costUp   = BathDiscHyb( BParamUp, SParamUp ; outputlevel=1)
    BParamDnNew, niterDn, costDn   = BathDiscHyb( BParamDn, SParamDn ; outputlevel=1)
    enewup, vnewup  = BathParamReshape( BParamUpNew, nbathHalf )
    enewdn, vnewdn  = BathParamReshape( BParamDnNew, nbathHalf )
    println("Parameters (inital) :")
    @show BParamUp
    @show BParamDn
    println("Parameters (optimized) :")
    @show BParamUpNew
    @show BParamDnNew
    println("")

    ebathlnew       = zero( ebathl )
    Vilnew          = zero( Vil )
    ebathlnew[1:2:end]              = deepcopy( enewup )
    ebathlnew[2:2:end]              = deepcopy( enewdn )
    Vilnew[ IndOrbUp, IndBathUp ]   = deepcopy( vnewup )
    Vilnew[ IndOrbDn, IndBathDn ]   = deepcopy( vnewdn )
    return ebathlnew, Vilnew, niterUp+niterDn, costUp+costDn
end

function BathDiscHybAll( ebathl, Vil, hybiw_grid, FreqGrid )
    nspinorb, nbath = size( Vil )
    BParam  = BathParamFlatten( ebathl, Vil )
    SParam    = ( 
                    hybiw_grid,
                    FreqGrid, 
                    nspinorb,
                    nbath
                    )
    BParamNew, niter, cost   = BathDiscHyb( BParam, SParam ; outputlevel=1)
    ebathlnew, Vilnew  = BathParamReshape( BParamNew, nbath )
    println("")
    println("Parameters (inital) :")
    @show BParam
    println("")
    println("Parameters (optimized) :")
    @show BParamNew
    println("")

    return deepcopy(ebathlnew), deepcopy(Vilnew), niter, cost
end

function BathSymPH( ebathl, Vil ; outputlevel=0 ) 
    nbath   = size(Vil)[2]
    BParamPH  = BathParamFlatten( ebathl[1:div(end,2)], Vil )
    enew, vnew  = BathParamReshapePH( BParamUpPH, nbath )
    if outputlevel > 0 
        println( "" )
        println( "ebathl before PH-sym : $(ebathl)" ) 
        println( "ebathl after  PH-sym : $(enew)" ) 
        # println( "Vil    before PH-sym : " ) ; writedlm(stdout, Vil, " ")
        # println( "Vil    after  PH-sym : " ) ; writedlm(stdout, vnew, " ")
        println( "" )
    end
    return enew, vnew
end

@inline function ShowBathParam( ebathl, Vil, SpinSep::Bool ; prefix="", affix="" )
    if SpinSep
        ShowBathParamSpin( ebathl, Vil )
    else
        ShowBathParam( ebathl, Vil )
    end
end

function ShowBathParam( ebathl, Vil ; prefix="", affix="" )
    println( "$(prefix)" )
    println( "ebathl : $(ebathl)" ) 
    println( "Vil    : " ) ; writedlm(stdout, Vil, " ")
    println( "$(affix)" )
end

@inline ShowBathParam( BP::BathParams ; prefix="", affix="" ) = ShowBathParam(BP.ebathl, BP.Vil; prefix=prefix, affix=affix)

function ShowBathParam( BPV::Vector{BathParams} ; prefix="", affix="" )
    for (ibp, bp_i) in Iterators.enumerate(BPV)
        ShowBathParam( bp_i; prefix="[ind_imp=$(ibp)]",affix=affix)
    end
end

function ShowBathParamSpin( ebathl, Vil )
    nbath       = length(ebathl)
    IndBathUp, IndBathDn = GetOrbUpDn( nbath )
    ebathlUp    = ebathl[IndBathUp]
    ebathlDn    = ebathl[IndBathDn]
    VilUp       = Vil[:,IndBathUp]
    VilDn       = Vil[:,IndBathDn]
    ShowBathParam( ebathlUp, VilUp, prefix="BathUp::" )
    ShowBathParam( ebathlDn, VilDn, prefix="BathDn::" )
end

function ReadBathParam!( ebathl, Vil ; fname="./BathParams_in.dat" )
    bp_raw  = readdlm( fname )
    ebathl[:]   .= deepcopy(bp_raw[:,1])
    Vil[:,:]    .= transpose(deepcopy(bp_raw[:,2:end]))
    println( "ebathl and Vil are loaded from '$(fname)'." )
    return 0 
end

function ReadBathParam!( BPV::Vector{BathParams} ; fname="./BathParams_in.dat", outputlevel=1 )
    bp_raw  = readdlm( fname )
    Nimp    = length(BPV)
    norbimp, nbathimp   = size(BPV[1].Vil)
    nrow, ncol  = size(bp_raw)
    if (nrow != (Nimp*nbathimp)) & ((ncol-1)!=norbimp)
        error("BathParams : dimensions are not compatible with the input-file. ($(Nimp)*$(nbathimp), $(norbimp)+1) != $(size(bp_raw)) ")
    end

    IndBathImp  = Iterators.partition( 1:nrow, nbathimp )
    @show IndBathImp
    for (indimp, indbath) in Iterators.enumerate(IndBathImp)
        BPV[indimp].ebathl[:]   .= deepcopy(bp_raw[indbath,1])
        BPV[indimp].Vil[:,:]    .= transpose(deepcopy(bp_raw[indbath,2:end]))
    end
    if outputlevel>0
        println( "ebathl and Vil are loaded from '$(fname)'." )
    end
    return 0 
end

function WriteBathParam( ebathl, Vil ; fname="./BathParams" )
    fopen   = open(fname, "w")
    nspinorb, nbath = size( Vil )
    for i in 1:nbath
        writedlm( fopen, [ ebathl[i] Vil[:,i]... ], " " )
    end
    close(fopen)
end

function WriteBathParam( BPV::Vector{BathParams} ; fname="./BathParams" )
    fopen   = open(fname, "w")
    for bp in BPV
        nspinorb, nbath = size( bp.Vil )
        for i in 1:nbath
            writedlm( fopen, [ bp.ebathl[i] bp.Vil[:,i]... ], " " )
        end
    end
    close(fopen)
end

function WriteBathParamInfo( ebathl, Vil, fname ; prefix="", affix="", iiter=0 )
    open(fname, "a" ) do io
        println(io, "$(prefix)" )
        println(io, "iter = $(iiter)" )
        println(io, "ebathl : $(ebathl)" ) 
        println(io, "Vil    : " ) ; writedlm(io, Vil, " ")
        println(io, "" )
    end
end

function WriteBathEnergy( ebathl, fname ; iiter=0 )
    open(fname, "a" ) do io
        write( io, format("{:3d} ",iiter) )
        writedlm( io, transpose(ebathl), " " )
    end
end

function WriteBathHyb( Vil, fname ; iiter=0)
    open(fname, "a" ) do io
        write( io, format("{:3d} ",iiter) )
        writedlm( io, transpose([Vil...]), " " )
    end
end

function WriteBathAll( BPV::Vector{BathParams} ; iiter=0 )
    fnameBathInfo   = "bath_info.dat"
    fnameBathE      = "bath_energy.dat"
    fnameBathV      = "bath_hyb.dat"

    ebathlAll   = vcat(  [BP.ebathl for BP in BPV]... )
    VilAll      = vcat(  [BP.Vil    for BP in BPV]... )
    WriteBathParamInfo( ebathlAll, VilAll,  fnameBathInfo ; iiter=iiter )
    WriteBathEnergy(    ebathlAll,          fnameBathE ; iiter=iiter )
    WriteBathHyb(       VilAll,             fnameBathV ; iiter=iiter )
end

function WriteBathAll( ebathl, Vil ; iiter=0 )
    fnameBathInfo   = "bath_info.dat"
    fnameBathE      = "bath_energy.dat"
    fnameBathV      = "bath_hyb.dat"

    WriteBathParamInfo( ebathl, Vil,    fnameBathInfo ; iiter=iiter )
    WriteBathEnergy(    ebathl,         fnameBathE ; iiter=iiter )
    WriteBathHyb(       Vil,            fnameBathV ; iiter=iiter )
end


function InitBathParam!( ebathl, Vil ; outputlevel=0 )
    nspinorb, nbath = size(Vil)
    nbathHalf   = div( nbath, 2 )
    nspinorbHalf    = div( nspinorb, 2 )
    @show size(Vil)
    ebathl[1:2:end]     = collect(LinRange( -0.8, 0.9, nbathHalf ))
    ebathl[2:2:end]     = collect(LinRange( -0.8, 0.9, nbathHalf ))
    for i in 1:nspinorbHalf
        for j in 1:nbathHalf
            Vil[2*i-1,2*j-1]    = 1 - ebathl[2*i-1]^2
            Vil[2*i  ,2*j  ]    = 1 - ebathl[2*i  ]^2
        end
    end
    if outputlevel>0
        println( "ebathl and Vil are manually chosen." )
    end
end

function InitBathParam!( BPV::Vector{BathParams} ; outputlevel=1 )
    for (ibp, bp_i) in Iterators.enumerate(BPV)
        InitBathParam!( bp_i.ebathl, bp_i.Vil )
    end
    if outputlevel>0
        println( "ebathl and Vil are manually chosen." )
    end
end
