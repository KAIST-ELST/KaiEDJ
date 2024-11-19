using JLD2

export SolveEDBethe
export SolveEDBetheDegen
export SolveEDTbDegen
export SolveEDTbDegenNimp

function ShowGSSectorInfo( GSSectorInfo )
    println("")
    println("GSSectorInfo : ")
    for GSSectorInfo_i in GSSectorInfo
        iGSSector, Ei, Zi, veci = GSSectorInfo_i
        @show (iGSSector, Ei, Zi)
    end
    println("")
end

function SolveEDBetheDegen( Hybiw, ebathl_init, Vil_init, FreqGrid, ntot, nspinorb, chem_int, chem, ECorrOrb, opccaavec ;
                                D=1 , outputlevel=1, BathOpt=true, beta_boltz=0, TolEGS=1e-8, TolEGSsec=1e-3, TolBoltz=1e-5, GSSectorInfo=[],
                                IterTag=0  )
    if BathOpt
        ebathlnew, Vilnew   = BathDiscHybSpinPHev( ebathl_init, Vil_init, Hybiw, FreqGrid )
    else
        ebathlnew, Vilnew   = (ebathl_init, Vil_init)
    end
    if outputlevel > 0 
        #ShowBathParam( ebathlnew, Vilnew )
        ShowBathParamSpin( ebathlnew, Vilnew )
    end

    opcavec     = [ GetOpBathParam(ebathlnew, Vilnew, ibath->ibath+nspinorb),
                    GetOpMatrix(ECorrOrb-I*chem-I*chem_int, nspinorb) ,
                    ]
    ESecMin_arr    = SearchGSSector( ntot, opcavec, opccaavec ; outputlevel=1 )

    Nsector     = length(ESecMin_arr)
    iGSSector0  = argmin(ESecMin_arr)
    E0          = ESecMin_arr[iGSSector0]
    BoolSearchSector    = map( x -> x<TolEGSsec , ESecMin_arr .- E0 )
    IndSearchSector = collect(1:Nsector)[BoolSearchSector]
    println("Ground states are in the following sectors.")
    @show IndSearchSector
    @show E0
    println("")

    println("Find more ground states in the choosen sectors ...")
    if length(GSSectorInfo) == 0
        for iSearchSector in IndSearchSector
            esyssec_AR  = GetGSFromSectorLanczos( iSearchSector, ntot, opcavec, opccaavec ; outputlevel=1 )
            E_arr       = real(esyssec_AR[1])
            if outputlevel > 0 
                @show (iSearchSector, E_arr)
            end
            if beta_boltz == 0 
                BoolGSSectorInside  = map( x -> x<TolEGS , E_arr .- E0 )
                BoltzWeight         = [ 1 for Ei in E_arr ][BoolGSSectorInside]
                @show E_arr .- E0
                @show BoolGSSectorInside  
                println("")
            else
                BoolGSSectorInside  = map( x -> exp(-beta_boltz*x)>TolBoltz , E_arr .- E0 )
                BoltzWeight         = [ exp(-beta_boltz*(Ei-E0)) for Ei in E_arr ][BoolGSSectorInside]
                @show E_arr .- E0
                @show BoolGSSectorInside  
                println("")
            end
            IndGSSectorInside   = collect(1:length(E_arr))[BoolGSSectorInside]
            println("")
            # if length(IndGSSectorInside)==0
            #     IndGSSectorInside   = [1  ]
            #     BoltzWeight         = [1.0]
            # end
            for iGSSector in IndGSSectorInside
                push!( GSSectorInfo, [ (iSearchSector,iGSSector), E_arr[iGSSector], BoltzWeight[iGSSector], deepcopy(esyssec_AR[2][:,iGSSector]) ] )
            end
            println("GSSectorInfo (updated)::" )
            ShowGSSectorInfo( GSSectorInfo )
        end

        if outputlevel > 0
            if beta_boltz == 0
                TolEGSFromBoltz = nothing
            else 
                TolEGSFromBoltz = -log( TolBoltz ) / beta_boltz 
            end
            @show (TolEGS, TolBoltz, beta_boltz, TolEGSFromBoltz)
            println("GSSectorInfo (final)::" )
            ShowGSSectorInfo( GSSectorInfo )
        end

        if length(GSSectorInfo)==0
            error("GSSectorInfo is empty.")
        end
    end

    GSSectorBoltzArr    = [ GSSectorInfo_i[3] for GSSectorInfo_i in GSSectorInfo ]
    Zpart               = sum(GSSectorBoltzArr)
    GSSectorBoltzArr_norm   = GSSectorBoltzArr / Zpart
    @show GSSectorBoltzArr
    @show Zpart
    @show GSSectorBoltzArr_norm

    GimpDegen   = Matrix{ComplexF64}[ zeros(ComplexF64,nspinorb,nspinorb) for freq in FreqGrid ]
    for (iGSSectorInfo, GSSectorInfo_i) in Iterators.enumerate(GSSectorInfo)
        Gimp    = GetGreenImpurityFromGS(  nspinorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ) * GSSectorInfo_i[3]
        GimpDegen   .+= Gimp
    end
    gimpup  = GetijarrayFromVecMat( GimpDegen, 1, 1 )
    gimpdn  = GetijarrayFromVecMat( GimpDegen, 2, 2 )
    
    G0imp   = GetGreenDiscGrid( ebathlnew, Vilnew, FreqGrid, ECorrOrb-I*chem )
    g0impup = GetijarrayFromVecMat( G0imp, 1, 1 )
    g0impdn = GetijarrayFromVecMat( G0imp, 2, 2 )

    selfup  = GetSelf( g0impup, gimpup )
    selfdn  = GetSelf( g0impdn, gimpdn )

    
    # G0newimp    = GetGreenImpGrid( zeros(nspinorb,nspinorb), (D*D / 2. / 2. ) * GimpDegen, FreqGrid )
    GbethenewUp = GetGzBetheFromSelf( FreqGrid, selfup )
    GbethenewDn = GetGzBetheFromSelf( FreqGrid, selfdn )

    Gbethenew   = GetVecMatFromDiag( [GbethenewUp, GbethenewDn] )
    Selfnew     = GetVecMatFromDiag( [selfup, selfdn] )

    return GimpDegen, G0imp, Gbethenew, (D*D / 2. / 2. )*Gbethenew, Selfnew, ebathlnew, Vilnew, GSSectorInfo
end


function SolveEDBethe( Hybiw, ebathl_init, Vil_init, FreqGrid, ntot, nspinorb, chem_int, chem, ECorrOrb, opccaavec ;
                        D=1 , outputlevel=1, BathOpt=true, IterTag=0 )
    if BathOpt
        ebathlnew, Vilnew   = BathDiscHybSpinPHev( ebathl_init, Vil_init, Hybiw, FreqGrid )
    else
        ebathlnew, Vilnew   = (ebathl_init, Vil_init)
    end
    if outputlevel > 0 
        #ShowBathParam( ebathlnew, Vilnew )
        ShowBathParamSpin( ebathlnew, Vilnew )
    end

    opcavec     = [ GetOpBathParam(ebathlnew, Vilnew, ibath->ibath+nspinorb),
                    GetOpMatrix(ECorrOrb-I*chem-I*chem_int, nspinorb) ,
                    ]
    Emin_arr    = SearchGSSector( ntot, opcavec, opccaavec ; outputlevel=1 )
    iGSSector   = argmin(Emin_arr)
    esyssec_AR  = GetGSFromSector( iGSSector, ntot, opcavec, opccaavec ; outputlevel=1 )

    ieval        = 1
    IndGSSector = [ [ iGSSector, esyssec_AR[1][ieval], 1.0, esyssec_AR[2][:,ieval] ] ]

    Gimp    = GetGreenImpurityFromGS(  nspinorb, IndGSSector[1], ntot, opcavec, opccaavec, FreqGrid )
    gimpup  = GetijarrayFromVecMat( Gimp, 1, 1 )
    gimpdn  = GetijarrayFromVecMat( Gimp, 2, 2 )
    
    G0imp   = GetGreenDiscGrid( ebathlnew, Vilnew, FreqGrid, ECorrOrb-I*chem )
    g0impup = GetijarrayFromVecMat( G0imp, 1, 1 )
    g0impdn = GetijarrayFromVecMat( G0imp, 2, 2 )

    selfup  = GetSelf( g0impup, gimpup )
    selfdn  = GetSelf( g0impdn, gimpdn )

    
    # G0newimp    = GetGreenImpGrid( zeros(nspinorb,nspinorb), (D*D / 2. / 2. ) * Gimp, FreqGrid )
    GbethenewUp = GetGzBetheFromSelf( FreqGrid, selfup )
    GbethenewDn = GetGzBetheFromSelf( FreqGrid, selfdn )

    Gbethenew   = GetVecMatFromDiag( [GbethenewUp, GbethenewDn] )
    Selfnew     = GetVecMatFromDiag( [selfup, selfdn] )

    return Gimp, G0imp, Gbethenew, (D*D / 2. / 2. )*Gbethenew, Selfnew, ebathlnew, Vilnew
end


function SolveEDTbDegen( hkarr, karr, wkarr, Hybiw, ebathl_init, Vil_init, FreqGrid, ntot, nspinorb, chem_int, chem, ECorrOrb, opccaavec ;
                            outputlevel=1, BathOpt=true, beta_boltz=0, TolGS=1e-9, TolEGS=1e-7, TolEGSsec=1e-2, TolBoltz=1e-5, GSSectorInfo=[], 
                            SpinSep=true, SelfTrev=true, SelfDiag=true, GimpDiag=true, IterTag=0,
                            IndCorrOrb=[], 
                            show_obs=true  )
    TOSolver = TimerOutput()
    nIterOpt    = -1 
    costFinal   = -1
    @timeit TOSolver "BathOpt" if BathOpt
        if SpinSep
            ebathlnew, Vilnew, nIterOpt, costFinal   = BathDiscHybSpin( ebathl_init, Vil_init, Hybiw, FreqGrid )
        else
            ebathlnew, Vilnew, nIterOpt, costFinal   = BathDiscHyb_grid( ebathl_init, Vil_init, Hybiw, FreqGrid )
        end
    else
        ebathlnew, Vilnew   = (ebathl_init, Vil_init)
    end
    println("BathOpt-done")
    WriteBathParam( ebathlnew, Vilnew ; fname="BathParams_out.dat")
    if outputlevel > 0 
        if SpinSep
            ShowBathParamSpin( ebathlnew, Vilnew )
        else
            ShowBathParam( ebathlnew, Vilnew )
        end
    end

    opcavec     = [ GetOpBathParam(ebathlnew, Vilnew, ibath->ibath+nspinorb),
                    GetOpMatrix(ECorrOrb-I*chem-I*chem_int, nspinorb) ,
                    ]
    @timeit TOSolver "SearchGSSector" ESecMin_arr    = SearchGSSector( ntot, opcavec, opccaavec ; outputlevel=outputlevel, nev=14 )
    println("SearchGSSector-done")

    Nsector     = length(ESecMin_arr)
    iGSSector0  = argmin(ESecMin_arr)
    E0          = ESecMin_arr[iGSSector0]
    BoolSearchSector    = map( x -> x<TolEGSsec , ESecMin_arr .- E0 )
    IndSearchSector = collect(1:Nsector)[BoolSearchSector]
    println("Ground states are in the following sectors.")
    @show IndSearchSector
    @show E0
    println("")

    println("Find more ground states in the choosen sectors ...")
    println("TolGS = $(TolGS)")
    @timeit TOSolver "GetGS" if length(GSSectorInfo) == 0
        for iSearchSector in IndSearchSector
            @time esyssec_AR  = GetGSFromSector( iSearchSector, ntot, opcavec, opccaavec ; outputlevel=outputlevel, nev=10, TolGS=TolGS )
            E_arr       = real(esyssec_AR[1])
            if outputlevel > 0 
                @show (iSearchSector, E_arr)
            end
            if beta_boltz == 0 
                BoolGSSectorInside  = map( x -> x<TolEGS , E_arr .- E0 )
                BoltzWeight         = [ 1 for Ei in E_arr ][BoolGSSectorInside]
                @show E_arr .- E0
                @show BoolGSSectorInside  
                println("")
            else
                BoolGSSectorInside  = map( x -> exp(-beta_boltz*x)>TolBoltz , E_arr .- E0 )
                BoltzWeight         = [ exp(-beta_boltz*(Ei-E0)) for Ei in E_arr ][BoolGSSectorInside]
                @show E_arr .- E0
                @show BoolGSSectorInside  
                println("")
            end
            IndGSSectorInside   = collect(1:length(E_arr))[BoolGSSectorInside]
            println("")
            # if length(IndGSSectorInside)==0
            #     IndGSSectorInside   = [1  ]
            #     BoltzWeight         = [1.0]
            # end
            for iGSSector in IndGSSectorInside
                push!( GSSectorInfo, [ (iSearchSector,iGSSector), E_arr[iGSSector], BoltzWeight[iGSSector], deepcopy(esyssec_AR[2][:,iGSSector]) ] )
            end
            println("GSSectorInfo (updated)::" )
            ShowGSSectorInfo( GSSectorInfo )
        end

        if outputlevel > 0
            if beta_boltz == 0
                TolEGSFromBoltz = nothing
            else 
                TolEGSFromBoltz = -log( TolBoltz ) / beta_boltz 
            end
            @show (TolEGS, TolBoltz, beta_boltz, TolEGSFromBoltz)
            println("GSSectorInfo (final)::" )
            ShowGSSectorInfo( GSSectorInfo )
        end

        if length(GSSectorInfo)==0
            error("GSSectorInfo is empty.")
        end
    end
    println("GetGS-done")

    GSSectorBoltzArr    = [ GSSectorInfo_i[3] for GSSectorInfo_i in GSSectorInfo ]
    Zpart               = sum(GSSectorBoltzArr)
    GSSectorBoltzArr_norm   = GSSectorBoltzArr / Zpart
    @show GSSectorBoltzArr
    @show Zpart
    @show GSSectorBoltzArr_norm

    GimpDegen   = Matrix{ComplexF64}[ zeros(ComplexF64,nspinorb,nspinorb) for freq in FreqGrid ]
    @timeit TOSolver "GetGimp" for (iGSSectorInfo, GSSectorInfo_i) in Iterators.enumerate(GSSectorInfo)
        Wboltz  = GSSectorBoltzArr_norm[iGSSectorInfo]
        Gimp    = GetGreenImpurityFromGS(  nspinorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ) * Wboltz
        GimpDegen   .+= Gimp
    end
    println("GetGimp-done")
    
    @timeit TOSolver "GetObservables" begin
    DMDegen     = zeros(ComplexF64,nspinorb,nspinorb)
    nGS         = length(GSSectorInfo)
    obs_info    = [[] for i=1:nGS+1]
    DMEach      = [zero(DMDegen) for i=1:nGS]
    nImp_info   = [[] for i=1:nGS]
    for (iGSSectorInfo, GSSectorInfo_i) in Iterators.enumerate(GSSectorInfo)
        Wboltz  = GSSectorBoltzArr_norm[iGSSectorInfo]
        DM      = GetDensityMatrix( GSSectorInfo_i, ntot, nspinorb ) 
        DMEach[iGSSectorInfo]   = deepcopy(DM)
        DMDegen .+= DM * Wboltz
        push!( obs_info[iGSSectorInfo], iGSSectorInfo, Wboltz, GetNtotNdiagFromDM( DM ), GetSzTotSzdiagFromDM( DM ) )
        nImpArr = GetNtot( GSSectorInfo_i, ntot, ntot )
        push!( nImp_info[iGSSectorInfo], iGSSectorInfo, Wboltz, sum(nImpArr),  nImpArr )
    end
    push!( obs_info[nGS+1], "tot", "1", GetNtotNdiagFromDM( DMDegen ), GetSzTotSzdiagFromDM( DMDegen ) )
    end
    observables = [ DMDegen, DMEach ]
    println("GetObservables-done")
    @timeit TOSolver "GimpDegenDiag" if GimpDiag
        ConstrainDiagAll!( GimpDegen, nspinorb )
    end
    
    @timeit TOSolver "Dyson_freq" begin
    G0imp       = GetGreenDiscGrid( ebathlnew, Vilnew, FreqGrid, ECorrOrb-I*chem )
    Selfimp     = GetSelf( G0imp, GimpDegen, chem_int )
    Hybimpnew   = GetHybFromGreenLocalGrid( G0imp, FreqGrid, ECorrOrb-I*chem )
    end
    println("Dyson_freq-done")

    @timeit TOSolver "Constraint" begin
    @timeit TOSolver "SelfTrev" if SelfTrev
        ConstrainTrevAllSz!( Selfimp, nspinorb )
    end
    @timeit TOSolver "SelfDiag" if SelfDiag
        ConstrainDiagAll!( Selfimp, nspinorb )
    end
    end
    
    @timeit TOSolver "Gkw_Dyson" begin
    @timeit TOSolver "Glattnew"  Glattnew    = GetGreenLocalFromHkSelfGrid( hkarr, wkarr, FreqGrid, Selfimp, collect(I(nspinorb)*(chem+chem_int)) )
    @timeit TOSolver "G0lattnew" G0lattnew   = GetGreenLocalFromDyson( Glattnew, Selfimp, chem_int ) 
    @timeit TOSolver "Hybnew"    Hybnew      = GetHybFromGreenLocalGrid( G0lattnew, FreqGrid, ECorrOrb - I*chem )
    end
    println("Gkw_Dyson-done")


    if show_obs 
        rho_imp     = GetDensityMatrixFromGreenImFreq( GimpDegen, FreqGrid )
        nimp_arr    = real(diag(rho_imp))
        nimptot     = sum(nimp_arr)

        rho_imp0    = GetDensityMatrixFromGreenImFreq( G0imp,     FreqGrid )
        nimp0_arr   = real(diag(rho_imp0))
        nimp0tot    = sum(nimp0_arr)

        rho_latt    = GetDensityMatrixFromGreenImFreq( Glattnew,  FreqGrid )
        nlatt_arr   = real(diag(rho_latt))
        nlatttot    = sum(nlatt_arr)

        rho_latt0   = GetDensityMatrixFromGreenImFreq( G0lattnew, FreqGrid )
        nlatt0_arr  = real(diag(rho_latt0))
        nlatt0tot   = sum(nlatt0_arr)

        n_info      = [ [nimptot, nimp_arr...], [nimp0tot, nimp0_arr], [nlatttot, nlatt_arr...], [nlatt0tot, nlatt0_arr...] ]

        println("Observables inside ground-states :")
        writedlm( stdout, obs_info )
        println("Occupations inside ground-states :")
        writedlm( stdout, nImp_info )
        println("Occupations from GF (Gimp, G0imp, Glattnew, G0lattnew) :")
        @show n_info
    end

    @show TOSolver
    return GimpDegen, G0imp, Hybimpnew, Glattnew, Hybnew, Selfimp, ebathlnew, Vilnew, GSSectorInfo, nIterOpt, costFinal, observables
end


function SolveEDTbDegenNimp( hkarr, karr, wkarr, Hybiw, BPV_init::Vector{BathParams}, FreqGrid, chem_int, chem, ECorrOrbNimp, opccaavec ;
                            D=1 , outputlevel=1, BathOpt=true, beta_boltz=0, TolGS=1e-9, TolEGS=1e-7, TolEGSsec=1e-2, TolBoltz=1e-5, GSSectorInfoNimp=[], 
                            SpinSep=true, SelfTrev=false, SelfDiag=false, GimpDiag=true, IterTag=0,
                            IndCorrOrbNimp=[],
                            show_obs=true  )
    TOSolver = TimerOutput()

    Nimp    = length(BPV_init)
    @assert length(IndCorrOrbNimp)!=0
    @assert length(IndCorrOrbNimp)==Nimp

    nspinorb, nbath = size(BPV_init[1].Vil)
    BPV_new                 = BathParams[]
    GSSectorInfoNimp_new    = []
    GSSectorInfoNimp_new    = []
    nIterOptNimp            = []
    costFinalNimp           = []
    status_occu_sz_all      = []
    status_occu_sz_imp      = []
    imp_infoeach            = []
    if GSSectorInfoNimp == []
        GSSectorInfoNimp    = [ [] for i=1:Nimp ]
    end

    dimhk       = size(hkarr[1])[1]
    SelfCluster = [ zero(hkarr[1]) for i in FreqGrid ]

    Hybiw_imparr    = [ [ Hyb_j[icorrorb,icorrorb] for Hyb_j in Hybiw ] for icorrorb in IndCorrOrbNimp ]
    for iimp in 1:Nimp
        init_println( "Impurity ($(iimp)/$(Nimp))" )
        ftail_imp   = "_imp$(iimp)"
        IndCorrOrb_imp  = IndCorrOrbNimp[iimp]
        Hybiw_imp       = Hybiw_imparr[iimp]
        BP_init         = BPV_init[iimp]
        nspinorb, nbath = size(BP_init.Vil)
        ntot            = nspinorb+nbath
        ECorrOrb_imp    = ECorrOrbNimp[IndCorrOrb_imp, IndCorrOrb_imp]

        ebathlnew, Vilnew, nIterOpt, costFinal  = BathDiscHybChoose( BP_init.ebathl, BP_init.Vil, Hybiw_imp, FreqGrid, BathOpt, SpinSep )
        if BathOpt 
            WriteBathParam( ebathlnew, Vilnew ; fname="BathParams_out$(ftail_imp).dat")
        end
        if outputlevel > 0 
            ShowBathParam( ebathlnew, Vilnew, SpinSep )
        end
        push!( BPV_new, BathParams(ebathlnew, Vilnew) )
        push!( nIterOptNimp, nIterOpt )
        push!( costFinalNimp, costFinal )

        Gimp, DMDegen, DMEach, GSSectorInfo    = SolveEDTbDegenOnlyGimp( ebathlnew, Vilnew, FreqGrid, ntot, nspinorb, chem_int, chem, ECorrOrb_imp, opccaavec ;
                            outputlevel=1, beta_boltz=beta_boltz, TolGS=TolGS, TolEGS=TolEGS, TolEGSsec=TolEGSsec, TolBoltz=TolBoltz, GSSectorInfo=GSSectorInfoNimp[iimp], 
                            GimpDiag=GimpDiag, IterTag=IterTag, ImpTag=iimp  )
        push!( GSSectorInfoNimp_new, GSSectorInfo )
        G0imp       = GetGreenDiscGrid( ebathlnew, Vilnew, FreqGrid, ECorrOrb_imp-I*chem )
        Selfimp     = GetSelf( G0imp, Gimp, chem_int )
        Hybimpnew   = GetHybFromGreenLocalGrid( G0imp, FreqGrid, ECorrOrb_imp-I*chem )

        if SelfTrev
            ConstrainTrevAllSz!( Selfimp, nspinorb )
        end
        if SelfDiag
            ConstrainDiagAll!( Selfimp, nspinorb )
        end

        map( ((sc,si),) -> sc[IndCorrOrb_imp,IndCorrOrb_imp] .= si , Iterators.zip( SelfCluster, Selfimp) )
        # SelfCluster[IndCorrOrb_imp,IndCorrOrb_imp] .= Selfimp

        if isnan( sum(sum(Selfimp)) )
            WriteVecMat( imag(FreqGrid),    GetijarrayFromVecMat( Selfimp,  1, 1 ), "./data/"*"Selfiw_imp$(iimp)_1_1_test.dat" )
            WriteVecMat( imag(FreqGrid),    GetijarrayFromVecMat( Selfimp,  1, 2 ), "./data/"*"Selfiw_imp$(iimp)_1_2_test.dat" )
            WriteVecMat( imag(FreqGrid),    GetijarrayFromVecMat( Selfimp,  2, 1 ), "./data/"*"Selfiw_imp$(iimp)_2_1_test.dat" )
            WriteVecMat( imag(FreqGrid),    GetijarrayFromVecMat( Selfimp,  2, 2 ), "./data/"*"Selfiw_imp$(iimp)_2_2_test.dat" )
            WriteMat( G0imp, "./data/"*"g0impiw_imp$(iimp)_test.dat" )
            WriteMat( Gimp,  "./data/"*"gimpiw_imp$(iimp)_test.dat" )
            error( " isnan( sum(sum(Selfimp)) ) = true" )
        end


        nimp_info   = GetDensityInfoFromGreenImFreq(Gimp,      FreqGrid)
        n0imp_info  = GetDensityInfoFromGreenImFreq(G0imp,     FreqGrid)
        szimpinfo   = GetSzTotSzdiagFromDM( DMDegen )
        if show_obs 
            println("Occupations from GF (Gimp, G0imp) and Sz [iimp=$(iimp)] :")
            @show (nimp_info, n0imp_info, szimpinfo)
        end
        push!( imp_infoeach, [nimp_info... n0imp_info... szimpinfo... ] )
        push!( status_occu_sz_imp, [nimp_info[1], n0imp_info[1], szimpinfo[1]] )

        GC.gc()
    end
    println("=== All impurities are solved. ===")
    println("")

    Glattnew    = GetGreenLocalFromHkSelfGrid( hkarr, wkarr, FreqGrid, SelfCluster, collect(I(dimhk)*(chem+chem_int)) )
    G0lattnew   = GetGreenLocalFromDyson( Glattnew, SelfCluster, chem_int ) 
    Hybnewiw    = GetHybFromGreenLocalGrid( G0lattnew, FreqGrid, ECorrOrbNimp - I*chem )

    nlatt_info  = GetDensityInfoFromGreenImFreq(Glattnew,  FreqGrid)
    n0latt_info = GetDensityInfoFromGreenImFreq(G0lattnew, FreqGrid)
    szlatt_arr  = nlatt_info[2:end] .* [ 0.5 * (-1)^(i%2+1) for i in 1:dimhk ]
    szlatt_info = [ sum(szlatt_arr) szlatt_arr... ]
    push!( status_occu_sz_all, nlatt_info[1], n0latt_info[1], szlatt_info[1],  sum(status_occu_sz_imp,dims=1)... )
    latt_info   = [ nlatt_info... n0latt_info... szlatt_info... ]
    if show_obs 
        println("Occupations from GF (Glattnew, G0lattnew) and Sz :")
        @show (nlatt_info, n0latt_info, szlatt_info)
    end

    return Glattnew, Hybnewiw, SelfCluster, BPV_new, GSSectorInfoNimp_new, nIterOptNimp, costFinalNimp, [ status_occu_sz_all, imp_infoeach, latt_info ]
end

function SolveEDTbDegenOnlyGimp( ebathl, Vil, FreqGrid, ntot, nspinorb, chem_int, chem, ECorrOrb, opccaavec ;
                            outputlevel=1, BathOpt=true, beta_boltz=0, TolGS=1e-9, TolEGS=1e-7, TolEGSsec=1e-2, TolBoltz=1e-5, GSSectorInfo=[], 
                            GimpDiag=true, IterTag=0, ImpTag=0,
                            IndCorrOrb=[] )
    TOSolver = TimerOutput()


    opcavec     = [ GetOpBathParam(ebathl, Vil, ibath->ibath+nspinorb),
                    GetOpMatrix(ECorrOrb-I*chem-I*chem_int, nspinorb) ,
                    ]

    @timeit TOSolver "GetGS" begin
        if length(GSSectorInfo) == 0
            GSSectorInfo    = GetGSAll( ntot, opcavec, opccaavec ; beta_boltz=beta_boltz, TolGS=TolGS, TolEGS=TolEGS, TolBoltz=TolBoltz, outputlevel=0)
        end
    end
    println("GetGS-done")

    # save("raw/gsinfo_imp$(ImpTag).jld2", "GSSectorInfo", GSSectorInfo)
    # println("SaveGS-done")
    # exit(1)

    @timeit TOSolver "GetGimp" begin
        GimpDegen   = GetGimpDegen(  GSSectorInfo, nspinorb, ntot, opcavec, opccaavec, FreqGrid )
    end
    println("GetGimp-done")
    
    @timeit TOSolver "GetObservables" begin
        DMDegen, DMEach, obs_info, nImp_info    = GetObservables( GSSectorInfo, nspinorb, ntot)
    end
    println("GetObservables-done")
    
    println("Observables inside ground-states :")
    writedlm( stdout, obs_info )
    println("Occupations inside ground-states :")
    writedlm( stdout, nImp_info )

    @timeit TOSolver "GimpDegenDiag" if GimpDiag
        ConstrainDiagAll!( GimpDegen, nspinorb )
    end

    return GimpDegen, DMDegen, DMEach, GSSectorInfo
end
