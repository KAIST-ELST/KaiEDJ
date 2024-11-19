using SparseArrays
using LinearAlgebra
using Arpack
using ThreadedSparseCSR
using SparseMatricesCSR
using Base.Threads

export Hamil
export ConstructHamil_ij!
export ConstructHamilAll_ij!
export MakeSparseComplexMatrix
export MakeSparseComplexMatrix!
# export CheckConvertSymmetric!
export SearchGSSector
export GetGSFromSector
export GetGSFromSectorLanczos
export GetGreenImpurityOrbFromGS
export GetGreenImpurityFromGS

export GetSaveCFcoeffFromGS

mutable struct Hamil
    MatSparse::SparseMatrixCSC
    Hamil( dim::Int ) = new( spzeros(ComplexF64, dim,dim) )
end

using TimerOutputs

function ConstructHamil_ij!( H::Hamil, opca::OpCA, HashF::Vector{Int64}, HashFInv::T ; outputlevel=0 , TO=nothing, ratbuffer=15 )  where T # HashF is equiv to the hashing function mapping i to binst
    dim = length( HashF )
    v       = opca.v
    l, m    = opca.i
    nthreads    = Threads.nthreads()
    nbuffer     = div( ratbuffer*dim , nthreads )
    I   = [ ones(Int, nbuffer) for i in 1:nthreads ]
    J   = [ ones(Int, nbuffer) for i in 1:nthreads ]
    V   = [ zeros(ComplexF64, nbuffer) for i in 1:nthreads ]
    nnz   = zeros(Int,nthreads)
    Threads.@threads for indBin_i in 1:dim 
        binst_i             = HashF[indBin_i]-1 
        FerSign, indBin_j   = OpCAIndBinst( opca, indBin_i, HashF, HashFInv )
        if (indBin_j > 0)
            id_th   = Threads.threadid()
            nnz[id_th] += 1
            nnz_th  = nnz[id_th]
            I[id_th][nnz_th]    = indBin_i
            J[id_th][nnz_th]    = indBin_j
            V[id_th][nnz_th]    = v * FerSign
            # @timeit TO "add_Hspth" Hspth[Threads.threadid()][indBin_i, indBin_j] += v * FerSign
            # if outputlevel > 0
            #     println( "(CdagC:$(opca.i))", " i=$(indBin_i) j=$(indBin_j) ; ", "binst_i : ", getbittaili(binst_i), " , binst_j : ", getbittaili(HashF[indBin_j]), " , FerSign : $(FerSign)" )
            # end
        end
    end
    for i in 1:nthreads
        resize!( I[i], nnz[i] )
        resize!( J[i], nnz[i] )
        resize!( V[i], nnz[i] )
    end
    I   = collect(Iterators.flatten(I)) # [ (I...)... ]
    J   = collect(Iterators.flatten(J)) # [ (J...)... ]
    V   = collect(Iterators.flatten(V)) # [ (V...)... ]
    nnz = collect(Iterators.flatten(nnz)) # [ (nnz...)... ]
    Hspth    = sparse(I, J, V, dim, dim)
    H.MatSparse += Hspth
    # for i in 1:nthreads
    #     resize!( I[i], nnz[i] )
    #     resize!( J[i], nnz[i] )
    #     resize!( V[i], nnz[i] )
    #     Hspth    = sparse(I[i],J[i],V[i], dim, dim)
    #     @timeit TO "add_H" H.MatSparse += Hspth
    # end
end

function ConstructHamil_ij!( H::Hamil, opccaa::OpCCAA, HashF::Vector{Int64}, HashFInv::T ; outputlevel=0 , TO=nothing, ratbuffer=15 ) where T
    dim = length( HashF )
    v               = opccaa.v
    l1, l2, l3, l4  = opccaa.i
    nthreads    = Threads.nthreads()
    nbuffer     = div( ratbuffer*dim , nthreads )
    I   = [ ones(Int, nbuffer) for i in 1:nthreads ]
    J   = [ ones(Int, nbuffer) for i in 1:nthreads ]
    V   = [ zeros(ComplexF64, nbuffer) for i in 1:nthreads ]
    nnz   = zeros(Int,nthreads)
    Threads.@threads for indBin_i in 1:dim 
        binst_i             = HashF[indBin_i]-1 
        FerSign, indBin_j   = OpCCAAIndBinst( opccaa, indBin_i, HashF, HashFInv )
        # if outputlevel > 0 
        #     try 
        #         println( "(CdagCdagCC:$(opccaa.i))", " i=$(indBin_i) j=$(indBin_j) ; ", "binst_i : ", getbittaili(binst_i), " , binst_j : ", getbittaili(HashF[indBin_j]), " , FerSign : $(FerSign)" )
        #     catch
        #         println( "(CdagCdagCC:$(opccaa.i))", " i=$(indBin_i) j=$(indBin_j) ; ", "binst_i : ", getbittaili(binst_i) )
        #     end
        # end
        if (indBin_j > 0)
            id_th   = Threads.threadid()
            nnz[id_th] += 1
            nnz_th  = nnz[id_th]
            I[id_th][nnz_th] = indBin_i
            J[id_th][nnz_th] = indBin_j
            V[id_th][nnz_th] = v * FerSign
            # @timeit TO "add_Hspth" Hspth[Threads.threadid()][indBin_i, indBin_j] += v * FerSign
            # println( "(CdagCdagCC:$(opccaa.i))", " i=$(indBin_i) j=$(indBin_j) ; ", "binst_i : ", getbittaili(binst_i), " , binst_j : ", getbittaili(HashF[indBin_j]))
        end
    end
    for i in 1:nthreads
        resize!( I[i], nnz[i] )
        resize!( J[i], nnz[i] )
        resize!( V[i], nnz[i] )
    end
    I   = collect(Iterators.flatten(I)) # [ (I...)... ]
    J   = collect(Iterators.flatten(J)) # [ (J...)... ]
    V   = collect(Iterators.flatten(V)) # [ (V...)... ]
    nnz = collect(Iterators.flatten(nnz)) # [ (nnz...)... ]
    Hspth    = sparse(I, J, V, dim, dim)
    H.MatSparse += Hspth
    # for i in 1:nthreads
    #     resize!( I[i], nnz[i] )
    #     resize!( J[i], nnz[i] )
    #     resize!( V[i], nnz[i] )
    #     Hspth    = sparse(I[i],J[i],V[i], dim, dim)
    #     @timeit TO "add_H" H.MatSparse += Hspth
    # end
end

function MakeSparseComplexMatrix( M ) 
    mr  = sparse(real(M))
    mi  = sparse(imag(M))
    return sparse( mr + mi *im )
end

function MakeSparseComplexMatrix!( M ) 
    mr  = sparse(real(M))
    mi  = sparse(imag(M))
    M   = sparse( mr + mi *im )
end

## It is NOT working somehow ##
# function CheckConvertSymmetric!( M ) 
#     if issymmetric(M)
#         mr  = sparse(real(M))
#         M   = deepcopy(sparse(mr) * 0)
#         @show count(!iszero, M)
#         println( "CheckConvertSymmetric :: okay done.")
#     end
# end

## Usingle : on-site interaction tensor for a single-band ##
# Uijkl   = OpCCAA[]
# push!( Uijkl, OpCCAA( U, [1,2,2,1] ) )

## Ut2g : on-site interaction tensor for t2g shell ##
# Uijkl   = OpCCAA[]
# # pushing example : push!( Uijkl, OpCCAA( 1, [1,2,2,1] ) )
# push!( Uijkl, OpCCAA( U,             [2*i-1, 2*i  , 2*i,     2*i-1] ) ) ## t2g : intra-orbital density-density  ##
# push!( Uijkl, OpCCAA( U-2.0*JHund,   [2*i-1, 2*j  , 2*j,     2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i != j , opposite-spin
# push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i-1, 2*j-1, 2*j-1,   2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i <  j , up-spin
# push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i  , 2*j  , 2*j  ,   2*i  ] ) ) ## t2g : inter-orbital density-density  ## # i <  j , dn-spin
# push!( Uijkl, OpCCAA( -JHund,        [2*i-1, 2*j,   2*j-1,   2*i  ] ) ) ## t2g : pair-exchange                  ## # i != j , up|dn => dn|up
# push!( Uijkl, OpCCAA(  JHund,        [2*i-1, 2*i,   2*j  ,   2*j-1] ) ) ## t2g : pair-hopping                   ## # i != j , i,up|i,dn => j,up|j,dn
## Actual code ##
# Uijkl   = OpCCAA[]
# U       = 10.0
# JHund   = 1.0
# for i in 1:norb 
#     push!( Uijkl, OpCCAA( U,             [2*i-1, 2*i  , 2*i,     2*i-1] ) ) ## t2g : intra-orbital density-density  ##
#     for j in 1:norb 
#         if i!=j
#             push!( Uijkl, OpCCAA( U-2.0*JHund,   [2*i-1, 2*j  , 2*j,     2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i != j , opposite-spin
#             push!( Uijkl, OpCCAA( -JHund,        [2*i-1, 2*j,   2*j-1,   2*i  ] ) ) ## t2g : pair-exchange                  ## # i != j , up|dn => dn|up
#             push!( Uijkl, OpCCAA(  JHund,        [2*i-1, 2*i,   2*j  ,   2*j-1] ) ) ## t2g : pair-hopping                   ## # i != j , i,up|i,dn => j,up|j,dn
#         end
#         if j<i
#             push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i-1, 2*j-1, 2*j-1,   2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i <  j , up-spin
#             push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i  , 2*j  , 2*j  ,   2*i  ] ) ) ## t2g : inter-orbital density-density  ## # i <  j , dn-spin
#         end
#     end
# end

function ConstructHamilAll_ij!( H::Hamil, opcaVV::Vector{Vector{OpCA}}, opccaaVV::Vector{Vector{OpCCAA}}, HashF, HashFInv ; outputlevel=0 )  # HashF is equiv to the hashing function mapping i to binst
    for opcaV in opcaVV
        for opca in opcaV
            if outputlevel > 0 
                println(" Constructing opca $(opca) ..." )
            end
            ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=outputlevel ) 
        end
    end
    for opccaaV in opccaaVV
        for opccaa in opccaaV
            if outputlevel > 0 
                println(" Constructing opccaa $(opccaa) ..." )
            end
            ConstructHamil_ij!( H, opccaa, HashF, HashFInv ; outputlevel=outputlevel ) 
        end
    end
end

function ConstructHamilAll_ij!( H::Hamil, opcaVec::Vector{OpCA}, opccaaVec::Vector{OpCCAA}, HashF, HashFInv ; outputlevel=0 )  # HashF is equiv to the hashing function mapping i to binst
    for opca in opcaVec
        if outputlevel > 0 
            println(" Constructing opca $(opca) ..." )
        end
        ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=outputlevel ) 
    end
    for opccaa in opccaaVec
        if outputlevel > 0 
            println(" Constructing opccaa $(opccaa) ..." )
        end
        ConstructHamil_ij!( H, opccaa, HashF, HashFInv ; outputlevel=outputlevel ) 
    end
end

function SearchGSSector( ntot, opcavec, opccaavec ; TolGS=1e-8, outputlevel=0, nev=12 ) 
    hashfsec    = HashingIntToSecarrBin(;norb=ntot, outputlevel=outputlevel-1)
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec ; outputlevel=outputlevel-1)
    hashfinv    = HashingInvSecarrBinToInt( hashfsec ; outputlevel=outputlevel-1)

    nsector     = length(hashfsec)
    Emin_arr    = zeros(nsector)
    println("GSSector Info : ")
    dimsec_arr  = [ length(hashfdum) for hashfdum in hashfsec ]
    for idum in Iterators.enumerate( dimsec_arr ) 
        @show idum
    end
    println("")

    if outputlevel > 0 
        println("Searching GSSector ...")
    end
    for isector in 1:nsector
        HashF       = hashfsec[isector]
        HashFInv    = hashfinv
        wfargs      = [ HashFInv ]
        DimSec      = length(HashF)
        if outputlevel > 0 
            println("$(isector)-th sector ; dim = $(DimSec)" )
        end 
        ######## Setting Up Hamiltonian Matrix of Each Sector from Operator #########
        HSec   = Hamil( DimSec )
        HSec.MatSparse  = zero(HSec.MatSparse)
        ConstructHamilAll_ij!( HSec, opcavec, opccaavec, HashF, HashFInv ; outputlevel=outputlevel )
    
        if DimSec < 2 
            Emin_arr[isector] = HSec.MatSparse[1,1]
        else 
            ######## Look-over Hamiltonian Diagonalization #########
            nvecLanc = DimSec > nev ? nev : DimSec-1
            try 
                # @time global esys_AR = eigs( HSec.MatSparse ; which=:SR , nev = nvecLanc, tol=TolGS )
                # println("eigs::Arpack-done")
                vec0    = normalize!(rand(Complex{Float64}, DimSec))
                global esys_AR = LanczosIter(HSec.MatSparse, vec0, nvecLanc)
                println("LanczosIter(H;isector=$(isector))-done")
            catch y
                @show y 
                @show isector
                if outputlevel > 1 
                    @show HashF
                end
                @show nev
                @show nvecLanc
                @show DimSec
                @show size(HSec.MatSparse)
                if outputlevel > 1 
                    @show HSec.MatSparse
                end
                if isa( y, ArgumentError )
                    @time global esys_eig = eigen( collect(HSec.MatSparse) )
                    global esys_AR = [esys_eig.values transpose(esys_eig.vectors)]
                    println("eigen(H)-done")
                else
                    println("WARNNING:: Cannot diagonalize $(isector)-th sector.")
                end
            end
            eval_sec    = esys_AR[1]
            Emin_arr[isector]   = real(eval_sec[1])
        end
        flush(stdout)
    end
    if outputlevel == 1 
        println( "" )
        println( "Minimum energy eigenvalue for each sector :" )
        for isector in 1:nsector
            println( "$(isector)-th sector Emin = $(Emin_arr[isector])" )
        end
        println( "" )
    end
    return Emin_arr
end

function GetGSFromSector( isector, ntot, opcavec, opccaavec ; TolGS=1e-8, outputlevel=0, nev=12 ) 
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )

    nsector     = length(hashfsec)

    HashF       = hashfsec[isector]
    HashFInv    = hashfinv
    wfargs      = [ HashFInv ]
    DimSec      = length(HashF)
    
    ######## Setting Up Hamiltonian Matrix of Each Sector from Operator #########
    HSec   = Hamil( DimSec )
    HSec.MatSparse  = zero(HSec.MatSparse)
    ConstructHamilAll_ij!( HSec, opcavec, opccaavec, HashF, HashFInv ; outputlevel=outputlevel )
    
    nvec    = nev > DimSec ? max(1,DimSec-2) : nev
    ncv     = min( 4*nev + 1, DimSec )
    @show (isector, nev, DimSec)
    global esys_AR = eigs( HSec.MatSparse ; which=:SR , nev = nvec, ncv=ncv, tol=TolGS, maxiter=7000 ) ## ith-eval => esAr[1][i] ; ith-evec => esAr[2][:,i]
    println("eigs::Arpack-done")
    evals    = esys_AR[1]
    if outputlevel > 0 
        neval   = length(evals)
        println( "" )
        println( "$(isector)-th sector energy eigenvalues :" )
        for ieval in 1:neval
            println( "$(isector)-th sector, $(ieval)-th energy = $(evals[ieval])" )
        end
        println( "" )
    end
    flush(stdout)
    return esys_AR
end

function LanczosIter( HMat, v0, nev ; nlanc=15, TolLancGS=1e-14, nvecComp=2 ) 
    vi  = deepcopy( v0 )
    viset_comp  = [ v0 for i=1:nvecComp ]
    evals_v = []
    evecs_v = []
    ThreadedSparseCSR.multithread_matmul(BaseThreads())
    Hcsr    = SparseMatrixCSR( Transpose(HMat) )
    for ilanc in 1:nlanc 
        HL, VL  = lanczos_algorithm( Hcsr, vi, nev )
        esys    = eigen(HL)
        evals_v = deepcopy( esys.values )
        evecs_v = deepcopy( VL[:,1:nev] * esys.vectors )
        evecs_v = mapslices( x->norm(x)>1e-14 ? normalize(x) : zero(x), evecs_v, dims=1 )
        vf      = evecs_v[:,1]
        vfset_comp   = [ evecs_v[:,i] for i=1:nvecComp ]
        overlap_set = [ abs2(dot(vix,vfx)) for (vix,vfx) in Iterators.zip(viset_comp,vfset_comp) ]

        if ilanc%5 ==0
            println( "iter_lanc ($(ilanc)/$(nlanc)) : eig = $(real(esys.values))" )
            flush(stdout)
        end

        if sum(map( x-> x > 1-TolLancGS, overlap_set)) > nvecComp
            @show overlap_set
            break
        end
        # if abs2(dot(vi,vf))> 1-TolLancGS
        #     #@show (conj(transpose(VL)) * VL )[1,:]
        #     break
        # end
        vi  = deepcopy(vf)
        viset_comp = deepcopy(vfset_comp)
    end
    return evals_v, evecs_v
end

function GetGSFromSectorLanczos( isector, ntot, opcavec, opccaavec ; TolGS=1e-12, outputlevel=0, nev=12 ) 
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )

    nsector     = length(hashfsec)

    HashF       = hashfsec[isector]
    HashFInv    = hashfinv
    wfargs      = [ HashFInv ]
    DimSec      = length(HashF)
    
    ######## Setting Up Hamiltonian Matrix of Each Sector from Operator #########
    HSec   = Hamil( DimSec )
    HSec.MatSparse  = zero(HSec.MatSparse)
    ConstructHamilAll_ij!( HSec, opcavec, opccaavec, HashF, HashFInv )
    
    vec0    = normalize!(rand(Complex{Float64}, DimSec))
    nvec    = DimSec > nev ? nev : DimSec
    global esys_AR = LanczosIter(HSec.MatSparse, vec0, nvec )
    evals    = esys_AR[1]
    if outputlevel > 0 
        neval   = length(evals)
        println( "" )
        println( "$(isector)-th sector energy eigenvalues :" )
        for ieval in 1:neval
            println( "$(isector)-th sector, $(ieval)-th energy = $(evals[ieval])" )
        end
        println( "" )
    end
    flush(stdout)
    return esys_AR
end

function GetGreenImpurityOrbFromGS(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ; NvecLanc=700, outputlevel=0 )
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    HashFAni    = hashfsec[isector-1]
    HashFCre    = hashfsec[isector+1]
    DimSec      = length(HashF)
    DimSecAni   = length(HashFAni)
    DimSecCre   = length(HashFCre)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )
    ciwf    = WF(DimSecAni, HashFAni, wfargs... )
    cdagjwf = WF(DimSecCre, HashFCre, wfargs... )

    ci      = OpA( 1.0, [iorb] ) 
    OpAWFadd( ci, wf, ciwf ) 
    cdagj   = OpC( 1.0, [jorb] ) 
    OpCWFadd( cdagj, wf, cdagjwf ) 
    # println( "$(iorb) $(jorb) ciwf    norm : $(norm(    ciwf.Probamp ))" )
    # println( "$(iorb) $(jorb) cdagjwf norm : $(norm( cdagjwf.Probamp ))" )


    HSecAni   = Hamil( DimSecAni )
    HSecCre   = Hamil( DimSecCre )
    HSecAni.MatSparse  = deepcopy(zero(HSecAni.MatSparse))
    HSecCre.MatSparse  = deepcopy(zero(HSecCre.MatSparse))

    ConstructHamilAll_ij!( HSecAni, opcavec, opccaavec, HashFAni, hashfinv ; outputlevel=outputlevel)
    ConstructHamilAll_ij!( HSecCre, opcavec, opccaavec, HashFCre, hashfinv ; outputlevel=outputlevel)

    normciwf    = norm(ciwf.Probamp)
    normcdagjwf = norm(cdagjwf.Probamp)
    ydat1   = GetGreenFromHCFcoeffA(    HSecAni.MatSparse, ciwf.Probamp,    FreqGrid, Eimp0 ; NvecLanc=NvecLanc)
    ydat2   = GetGreenFromHCFcoeffAdag( HSecCre.MatSparse, cdagjwf.Probamp, FreqGrid, Eimp0 ; NvecLanc=NvecLanc)
    if abs(normciwf)<1e-12
        ydat1   = [ zero(y_i) for y_i in ydat1 ] 
    end
    if abs(normcdagjwf)<1e-12
        ydat2   = [ zero(y_i) for y_i in ydat2 ] 
    end
    return ydat1 + ydat2
end

function GetGreenImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ; NvecLanc=700, cicoeff=1.0, cjcoeff=1.0, cdagicoeff=1.0, cdagjcoeff=1.0 , outputlevel=0)
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    HashFAni    = hashfsec[isector-1]
    HashFCre    = hashfsec[isector+1]
    DimSec      = length(HashF)
    DimSecAni   = length(HashFAni)
    DimSecCre   = length(HashFCre)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )
    ciwf    = WF(DimSecAni, HashFAni, wfargs... )
    cjwf    = WF(DimSecAni, HashFAni, wfargs... )
    cdagiwf = WF(DimSecCre, HashFCre, wfargs... )
    cdagjwf = WF(DimSecCre, HashFCre, wfargs... )

    ci      = OpA( cicoeff,     [iorb] ) ; OpAWFadd( ci, wf, ciwf ) 
    cj      = OpA( cjcoeff,     [jorb] ) ; OpAWFadd( cj, wf, cjwf ) 
    cdagi   = OpC( cdagicoeff,  [iorb] ) ; OpCWFadd( cdagi, wf, cdagiwf ) 
    cdagj   = OpC( cdagjcoeff,  [jorb] ) ; OpCWFadd( cdagj, wf, cdagjwf ) 
    # println( "$(iorb) $(jorb) ciwf    norm : $(norm(    ciwf.Probamp ))" )
    # println( "$(iorb) $(jorb) cdagiwf norm : $(norm( cdagiwf.Probamp ))" )



    HSecAni   = Hamil( DimSecAni )
    HSecCre   = Hamil( DimSecCre )
    HSecAni.MatSparse  = deepcopy(zero(HSecAni.MatSparse))
    HSecCre.MatSparse  = deepcopy(zero(HSecCre.MatSparse))

    ConstructHamilAll_ij!( HSecAni, opcavec, opccaavec, HashFAni, hashfinv ; outputlevel=outputlevel)
    ConstructHamilAll_ij!( HSecCre, opcavec, opccaavec, HashFCre, hashfinv ; outputlevel=outputlevel)

    ydat1   = GetGreenFromHCFcoeffA(    HSecAni.MatSparse, ciwf.Probamp + cjwf.Probamp,         FreqGrid, Eimp0 ; NvecLanc=NvecLanc)
    ydat2   = GetGreenFromHCFcoeffAdag( HSecCre.MatSparse, cdagiwf.Probamp + cdagjwf.Probamp,   FreqGrid, Eimp0 ; NvecLanc=NvecLanc)
    return ydat1 + ydat2
end

function GetGreenImpurityFromGS(  norb::Int, GSSectorInfo_i, ntot::Int, opcavec, opccaavec, FreqGrid ; NvecLanc=700 )
    gorb    = [ [[] for i in 1:norb] for j in 1:norb]
    for iorb in 1:norb
        # delta_ij    = iorb==jorb ? 1.0 : 0.0
        gorb[iorb][iorb] = GetGreenImpurityOrbFromGS(  iorb, iorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ; NvecLanc=NvecLanc )
    end
    for iorb=1:norb, jorb=iorb+1:norb
        Gij_1_1_1_1     =  GetGreenImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ; NvecLanc=NvecLanc)
        Gij_i_i_1_m1    =  GetGreenImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ; NvecLanc=NvecLanc, cjcoeff=-1.0*im, cdagjcoeff=1.0*im ) * im
        gorb[iorb][jorb]    = ( Gij_1_1_1_1 + Gij_i_i_1_m1 - (1.0+im)*gorb[iorb][iorb] - (1.0+im)*gorb[jorb][jorb] ) * 0.5
        gorb[jorb][iorb]    = ( Gij_1_1_1_1 - Gij_i_i_1_m1 - (1.0-im)*gorb[iorb][iorb] - (1.0-im)*gorb[jorb][jorb] ) * 0.5
    end
    nfreq   = length(FreqGrid)
    GorbFormatted    = Matrix{ComplexF64}[ zeros(ComplexF64,norb,norb) for ifreq in 1:nfreq ]
    for (iorb, jorb, ifreq) in Iterators.product( 1:norb, 1:norb, 1:nfreq )
        GorbFormatted[ifreq][iorb,jorb] = deepcopy(gorb[iorb][jorb][ifreq])
    end
    return GorbFormatted
end

function GetCFcoeffHAPsi( H, AVec ; NvecLanc=800 )
    AbsAVec     = norm(AVec)^2
    TriMat, Vecs = lanczos_algorithm( H, normalize(AVec), NvecLanc)
    return [ AbsAVec, diag(TriMat), diag(TriMat,1) ]
end

function GetCFcoeffImpurityOrbFromGS(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec ; NvecLanc=700, outputlevel=0 )
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    HashFAni    = hashfsec[isector-1]
    HashFCre    = hashfsec[isector+1]
    DimSec      = length(HashF)
    DimSecAni   = length(HashFAni)
    DimSecCre   = length(HashFCre)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )
    ciwf    = WF(DimSecAni, HashFAni, wfargs... )
    cdagiwf = WF(DimSecCre, HashFCre, wfargs... )

    ci      = OpA( 1.0, [iorb] ) 
    OpAWFadd( ci, wf, ciwf ) 
    cdagi   = OpC( 1.0, [jorb] ) 
    OpCWFadd( cdagi, wf, cdagiwf ) 
    # println( "$(iorb) $(jorb) ciwf    norm : $(norm(    ciwf.Probamp ))" )
    # println( "$(iorb) $(jorb) cdagiwf norm : $(norm( cdagiwf.Probamp ))" )


    HSecAni   = Hamil( DimSecAni )
    HSecCre   = Hamil( DimSecCre )
    HSecAni.MatSparse  = deepcopy(zero(HSecAni.MatSparse))
    HSecCre.MatSparse  = deepcopy(zero(HSecCre.MatSparse))

    ConstructHamilAll_ij!( HSecAni, opcavec, opccaavec, HashFAni, hashfinv ; outputlevel=outputlevel)
    ConstructHamilAll_ij!( HSecCre, opcavec, opccaavec, HashFCre, hashfinv ; outputlevel=outputlevel)

    normciwf    = norm(ciwf.Probamp)
    normcdagiwf = norm(cdagiwf.Probamp)

    CFDataAni   = GetCFcoeffHAPsi( HSecAni.MatSparse, ciwf.Probamp,   ; NvecLanc=NvecLanc)
    CFDataCre   = GetCFcoeffHAPsi( HSecCre.MatSparse, cdagiwf.Probamp ; NvecLanc=NvecLanc)
    if abs(normciwf)<1e-12
        CFDataAni[1]    = 0.0
    end
    if abs(normcdagiwf)<1e-12
        CFDataCre[1]    = 0.0
    end
    return CFDataAni, CFDataCre
end

function GetCFcoeffImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec ; NvecLanc=700, cicoeff=1.0, cjcoeff=1.0, cdagicoeff=1.0, cdagjcoeff=1.0 , outputlevel=0)
    hashfsec    = HashingIntToSecarrBin(;norb=ntot )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    HashFAni    = hashfsec[isector-1]
    HashFCre    = hashfsec[isector+1]
    DimSec      = length(HashF)
    DimSecAni   = length(HashFAni)
    DimSecCre   = length(HashFCre)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )
    ciwf    = WF(DimSecAni, HashFAni, wfargs... )
    cjwf    = WF(DimSecAni, HashFAni, wfargs... )
    cdagiwf = WF(DimSecCre, HashFCre, wfargs... )
    cdagjwf = WF(DimSecCre, HashFCre, wfargs... )

    ci      = OpA( cicoeff,     [iorb] ) ; OpAWFadd( ci, wf, ciwf ) 
    cj      = OpA( cjcoeff,     [jorb] ) ; OpAWFadd( cj, wf, cjwf ) 
    cdagi   = OpC( cdagicoeff,  [iorb] ) ; OpCWFadd( cdagi, wf, cdagiwf ) 
    cdagj   = OpC( cdagjcoeff,  [jorb] ) ; OpCWFadd( cdagj, wf, cdagjwf ) 
    # println( "$(iorb) $(jorb) ciwf    norm : $(norm(    ciwf.Probamp ))" )
    # println( "$(iorb) $(jorb) cdagiwf norm : $(norm( cdagiwf.Probamp ))" )



    HSecAni   = Hamil( DimSecAni )
    HSecCre   = Hamil( DimSecCre )
    HSecAni.MatSparse  = deepcopy(zero(HSecAni.MatSparse))
    HSecCre.MatSparse  = deepcopy(zero(HSecCre.MatSparse))

    ConstructHamilAll_ij!( HSecAni, opcavec, opccaavec, HashFAni, hashfinv ; outputlevel=outputlevel)
    ConstructHamilAll_ij!( HSecCre, opcavec, opccaavec, HashFCre, hashfinv ; outputlevel=outputlevel)

    CFDataAni   = GetCFcoeffHAPsi( HSecAni.MatSparse, ciwf.Probamp + cjwf.Probamp,      ; NvecLanc=NvecLanc)
    CFDataCre   = GetCFcoeffHAPsi( HSecCre.MatSparse, cdagiwf.Probamp + cdagjwf.Probamp ; NvecLanc=NvecLanc)
    return CFDataAni, CFDataCre
end


function GetSaveCFcoeffFromGS(  norb::Int, GSSectorInfo_i, ntot::Int, opcavec, opccaavec ; NvecLanc=700, IterTag=0 )
    CF_fdir     = "data/CF_coeff/"
    run( `mkdir $(CF_fdir)` )
    CF_fname    = CF_fdir * "coeff_i$(IterTag).h5"

    WriteHDF5( CF_fname, "norb", norb )
    Eimp0   = GSSectorInfo_i[2]
    WriteHDF5( CF_fname, "Eimp0", Eimp0 )

    CFcoeffDiagOrb       = [[] for i in 1:norb]
    for iorb=1:norb
        CFcoeffDiagOrb[iorb]    = GetCFcoeffImpurityOrbFromGS(  iorb, iorb, GSSectorInfo_i, ntot, opcavec, opccaavec ; NvecLanc=NvecLanc )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(iorb)/CFDataAni", CFcoeffDiagOrb[iorb][1] )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(iorb)/CFDataCre", CFcoeffDiagOrb[iorb][2] )
    end
    CFcoeff_11               = [ [[] for i in 1:norb] for j in 1:norb]
    CFcoeff_1ibar            = [ [[] for i in 1:norb] for j in 1:norb]
    for iorb=1:norb, jorb=iorb+1:norb
        CFcoeff_11[iorb][jorb]       =  GetCFcoeffImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec ; NvecLanc=NvecLanc)
        CFcoeff_1ibar[iorb][jorb]    =  GetCFcoeffImpurityOrbFromGSLinearcomb(  iorb, jorb, GSSectorInfo_i, ntot, opcavec, opccaavec ; NvecLanc=NvecLanc, cjcoeff=-1.0*im, cdagjcoeff=1.0*im )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_11/CFDataAni", CFcoeff_11[iorb][1] )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_11/CFDataCre", CFcoeff_11[iorb][2] )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_1ibar/CFDataAni", CFcoeff_1ibar[iorb][1] )
        WriteHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_1ibar/CFDataCre", CFcoeff_1ibar[iorb][2] )
    end
end

function GetGreenImpurityFromGSCFLoaded(  norb::Int, GSSectorInfo_i, ntot::Int, opcavec, opccaavec, FreqGrid ; NvecLanc=700, IterTag=0 )
    CF_fdir     = "data/CF_coeff/"
    CF_fname    = CF_fdir * "coeff_i$(IterTag).h5"
    if !isfile(CF_fname)
        error( "File existence :: $(CF_fname)" )
    end

    Eimp0   = ReadHDF5( CF_fname, "Eimp0")

    gorb    = [ [[] for i in 1:norb] for j in 1:norb]
    for iorb=1:norb
        CFDAni  = ReadHDF5_CFData( CF_fname, "$(iorb)_$(iorb)/CFDataAni")
        CFDCre  = ReadHDF5_CFData( CF_fname, "$(iorb)_$(iorb)/CFDataCre")
        GA = GetGreenFromHCFLoadedAni( CFDAni, FreqGrid, Eimp0 )
        GC = GetGreenFromHCFLoadedCre( CFDCre, FreqGrid, Eimp0 )
        gorb[iorb][iorb] = GA + GC
    end
    CFcoeff_11               = [ [[] for i in 1:norb] for j in 1:norb]
    CFcoeff_1ibar            = [ [[] for i in 1:norb] for j in 1:norb]
    for iorb=1:norb, jorb=iorb+1:norb
        CFDAni_11       = ReadHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_11/CFDataAni")
        CFDCre_11       = ReadHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_11/CFDataCre")
        CFDAni_1ibar    = ReadHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_1ibar/CFDataAni")
        CFDCre_1ibar    = ReadHDF5_CFData( CF_fname, "$(iorb)_$(jorb)_1ibar/CFDataCre")
        GA11    = GetGreenFromHCFLoadedAni( CFDAni_11, FreqGrid, Eimp0 )
        GC11    = GetGreenFromHCFLoadedCre( CFDCre_11, FreqGrid, Eimp0 )
        GA1ibar = GetGreenFromHCFLoadedAni( CFDAni_1ibar, FreqGrid, Eimp0 )
        GC1ibar = GetGreenFromHCFLoadedCre( CFDCre_1ibar, FreqGrid, Eimp0 )
        G11     = GA11 + GC11
        G1ibar  = GA1ibar + GC1ibar 
        gorb[iorb][jorb]    = ( G11 + im*G1ibar - (1.0+im)*gorb[iorb][iorb] - (1.0+im)*gorb[jorb][jorb] ) * 0.5
        gorb[jorb][iorb]    = ( G11 - im*G1ibar - (1.0-im)*gorb[iorb][iorb] - (1.0-im)*gorb[jorb][jorb] ) * 0.5
    end
    nfreq   = length(FreqGrid)
    GorbFormatted    = Matrix{ComplexF64}[ zeros(ComplexF64,norb,norb) for ifreq in 1:nfreq ]
    for (iorb, jorb, ifreq) in Iterators.product( 1:norb, 1:norb, 1:nfreq )
        GorbFormatted[ifreq][iorb,jorb] = deepcopy(gorb[iorb][jorb][ifreq])
    end
    return GorbFormatted
end


function FindGSSector( ntot, opcavec, opccaavec ; TolEGSsec=1e-2, outputlevel=0, nev=20 )
    ESecMin_arr    = SearchGSSector( ntot, opcavec, opccaavec ; outputlevel=outputlevel, nev=nev )
    Nsector     = length(ESecMin_arr)
    iGSSector0  = argmin(ESecMin_arr)
    E0          = ESecMin_arr[iGSSector0]
    BoolSearchSector    = map( x -> x<TolEGSsec , ESecMin_arr .- E0 )
    IndSearchSector = collect(1:Nsector)[BoolSearchSector]

    println("Ground states are in the following sectors.")
    @show IndSearchSector
    @show E0
    println("")
    return E0, IndSearchSector
end

function FindGSFromSector( IndSearchSector, E0, ntot, opcavec, opccaavec ; beta_boltz=0, TolGS=1e-8, TolEGS=1e-7, TolBoltz=1e-5, outputlevel=0, nev=12 ) 
    GSSectorInfo    = []
    for iSearchSector in IndSearchSector
        @time esyssec_AR  = GetGSFromSector( iSearchSector, ntot, opcavec, opccaavec ; outputlevel=outputlevel, nev=nev, TolGS=TolGS )
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

    GSSectorBoltzArr    = [ GSSectorInfo_i[3] for GSSectorInfo_i in GSSectorInfo ]
    Zpart               = sum(GSSectorBoltzArr)
    GSSectorBoltzArr_norm   = GSSectorBoltzArr / Zpart
    @show GSSectorBoltzArr
    @show Zpart
    @show GSSectorBoltzArr_norm

    
    for (GSSectorInfo_i, wboltz_norm_i) in Iterators.zip( GSSectorInfo, GSSectorBoltzArr_norm )
        push!( GSSectorInfo_i, wboltz_norm_i, ) 
    end

    return GSSectorInfo
end

function GetGSAll( ntot, opcavec, opccaavec ; TolEGSsec=1e-2, beta_boltz=beta_boltz, TolGS=TolGS, TolEGS=TolEGS, TolBoltz=TolBoltz, outputlevel=0)
    E0, IndSearchSector = FindGSSector( ntot, opcavec, opccaavec ; TolEGSsec=TolEGSsec, outputlevel=outputlevel, nev=20 )
    println("SearchGSSector-done")

    println("Find more ground states in the choosen sectors ...")
    GSSectorInfo    = FindGSFromSector( IndSearchSector, E0, ntot, opcavec, opccaavec ; 
                        beta_boltz=beta_boltz, TolGS=TolGS, TolEGS=TolEGS, TolBoltz=TolBoltz, outputlevel=0, nev=12 ) 
    return GSSectorInfo
end
