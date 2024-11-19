using LinearAlgebra
using DelimitedFiles
using FastGaussQuadrature
using Distributed
using TimerOutputs

export GetDeltaHybDiscGrid
export GetDeltaFromG0Grid
export GetVecMatFromDiag
export GetijarrayFromVecMat
export GetGreenImpGrid
export GetGreenDiscGrid
export GetImFreqValGrid
export GetDeltaHybDiscGridFromFlat
export GetDeltaHybDiscGridFromFlatPH
export GetDeltaHybDiscGridFromFlatPHev
export GetSelf
export GetDyson
export GetCost
export GetCostFromFlat
export GetCostFromFlatPH
export GetCostFromFlatPHev
export GetCostGlocinvFromFlat
export GetCostGlocinvFromFlatPH
export GetGreenFromHCFcoeff
export GetGreenFromHCFcoeffA
export GetGreenFromHCFcoeffAdag

export GetGreenLatticeTrFromDiagGrid
export GetGreenLatticeDiagFromHkSelfGrid
export GetGreenLatticeFromHkSelfGrid
export GetGreenLocalFromHkGrid
export GetGreenLocalFromHkSelfGrid
export GetGreenLocalFromHybGrid
export GetHybFromGreenLocalFromHkGrid
export GetHybFromGreenLocalGrid

export GetGimpDegen

export ConstrainTrev!
export ConstrainDiag!
export ConstrainDiagAll!

export GetDensityMatrixFromGreenImFreq
export GetDensityInfoFromGreenImFreq

function GetImFreqValGrid( beta, ngrid ) 
    return collect(range( pi, step=2*pi, length=ngrid )) / beta
end

function GetVecMatFromDiag( Mdiagarr_grid ) 
    norb    = length( Mdiagarr_grid )
    nfreq   = length( Mdiagarr_grid[1] )
    M       = Matrix{ComplexF64}[ zeros(norb,norb) for i in 1:nfreq ]
    for ifreq in 1:nfreq
        for iorb in 1:norb
            M[ifreq][iorb,iorb]    = Mdiagarr_grid[iorb][ifreq]
        end
    end
    return M
end

function GetijarrayFromVecMat( VMat, i ,j ) 
    return [ m[i,j] for m in VMat ]
end

@inline function GetDeltaFromG0Grid( G0, wgrid )
    nfreq   = length(wgrid) 
    return [ wgrid[ifreq]*I - inv(G0[ifreq]) for ifreq in 1:nfreq ]
end

function GetDeltaHybDiscGrid( ebathl, Vil, wgrid )
    return [ GetDeltaHybDisc( ebathl, Vil, w ) for w in wgrid ]
end

function GetDeltaHybDisc( ebathl, Vil, w )
    DimOrb, DimBath = size(Vil)
    Mij             = zeros(ComplexF64, DimOrb, DimOrb )
    for (i,j) in Iterators.product( 1:DimOrb, 1:DimOrb )
        # println( "$i, $j :",  sum([ ( conj(Vil[i,ibath]) * Vil[j,ibath] )/( w - ebathl[ibath] ) for ibath in 1:DimBath]) )
        Mij[i,j]    = sum([ ( conj(Vil[i,ibath]) * Vil[j,ibath] )/( w - ebathl[ibath] ) for ibath in 1:DimBath])
    end
    return Mij 
end

function GetGreenDiscGrid( ebathl, Vil, wgrid, Eorb )
    return [ inv(w*I -Eorb - GetDeltaHybDisc( ebathl, Vil, w )) for w in wgrid ]
end

function GetGreenImpGrid( eij, Hybij_w, wgrid )
    ngrid   = length(wgrid)
    return [ GetGreenImp( eij, Hybij_w[ind_w], wgrid[ind_w] ) for ind_w in 1:ngrid ]
end

@inline GetGreenImp( eij, Hybij, w )    = inv(-eij + I * w - Hybij)


@inline function GetCost( GijOrig_w, GijFit_w )
    diff_w = GijOrig_w - GijFit_w 
    return sum([ norm(diff)^2 for diff in diff_w ])
end

@inline function GetCost( GijOrig_w::Vector{Matrix{ComplexF64}}, GijFit_w::Vector{Matrix{ComplexF64}} )
    diff_w = GijOrig_w - GijFit_w 
    return sum([ norm(diff)^2 for diff in diff_w ])
end

@inline function GetCost( GijOrig_w, GijFit_w, weight_w, ngrid )
    diff_w = GijOrig_w - GijFit_w 
    return sum([  norm(diff_w[ind_w])^2 * weight_w[ind_w]  for ind_w in 1:ngrid ])
end

@inline function GetSelfDyson( g0g::Tuple{Matrix{ComplexF64},Matrix{ComplexF64}} )
    return inv(g0g[1]) - inv(g0g[2])
end

@inline function GetSelfDyson( g0g::Tuple{ComplexF64,ComplexF64} )
    return inv(g0g[1]) - inv(g0g[2])
end

# @inline function GetSelf( g0_grid::Vector{T}, g_grid::Vector{T} ) where T
#     return GetSelfDyson.( collect(Iterators.zip( g0_grid, g_grid )) )
# end
@inline function GetSelf( g0_grid, g_grid )
    nfreq   = length(g0_grid)
    return [ inv(g0_grid[ifreq]) - inv(g_grid[ifreq])  for ifreq in 1:nfreq ]
end

@inline function GetSelf( g0_grid, g_grid, chem_int::ComplexF64 )
    nfreq   = length(g0_grid)
    return [ inv(g0_grid[ifreq]) - inv(g_grid[ifreq]) + I*chem_int  for ifreq in 1:nfreq ]
end
@inline GetSelf( g0_grid, g_grid, chem_int::Float64 )   = GetSelf( g0_grid, g_grid, chem_int + 0.0*im )

@inline function GetGreenLocalFromDyson( g_grid, self_grid, chem_int )
    nfreq   = length(g_grid)
    return [ inv( inv(g_grid[ifreq]) + self_grid[ifreq] - I*chem_int ) for ifreq in 1:nfreq ]
end


@inline function GetDeltaHybDiscGridFromFlat( bathparam, wgrid, DimOrb, DimBath )
    ebathl  = bathparam[1:DimBath]
    Vil     = bathparam[DimBath+1:end]
    Vil     = reshape( Vil, DimOrb, DimBath )
    return [ GetDeltaHybDiscFromFlat( ebathl, Vil, w, DimOrb, DimBath ) for w in wgrid ]
end

function GetDeltaHybDiscFromFlat( ebathl, Vil, w, DimOrb, DimBath )
    Mij             = zeros(ComplexF64, DimOrb, DimOrb )
    for (i,j) in Iterators.product( 1:DimOrb, 1:DimOrb )
        Mij[i,j]    += sum([ ( conj(Vil[i,ibath]) * Vil[j,ibath] )/( w - ebathl[ibath] ) for ibath in 1:DimBath])
    end
    return Mij 
end

@inline function GetDeltaHybDiscGridFromFlatPH( bathparamPH, wgrid, DimOrb, DimBath )
    DimBathHalf = div( DimBath , 2 )
    ebathl      = bathparamPH[1:DimBathHalf]
    if DimBath%2 < 1
        ebathl      = [ ebathl... reverse(-ebathl)... ]
    else
        ebathl      = [ ebathl... 0.0 reverse(-ebathl)... ]
    end
    Vil         = bathparamPH[DimBathHalf+1:end]
    Vil         = reshape( Vil, DimOrb, DimBath )
    return [ GetDeltaHybDiscFromFlat( ebathl, Vil, w, DimOrb, DimBath ) for w in wgrid ]
end


@inline function GetDeltaHybDiscGridFromFlatPHev( bathparamPHev, wgrid, DimOrb, DimBath )
    DimBathHalf     = div( DimBath , 2 )
    IndBathEnergy   = collect(1:DimBathHalf)
    IndBathHyb      = []
    ebathl          = bathparamPHev[IndBathEnergy]
    if DimBath%2 < 1
        ebathl      = [ ebathl... reverse(-ebathl)... ]
        IndBathHyb  = [ IndBathEnergy..., reverse(IndBathEnergy)... ]
    else
        ebathl      = [ ebathl... 0.0 reverse(-ebathl)... ]
        IndBathHyb  = [ IndBathEnergy..., DimBathHalf+1, reverse(IndBathEnergy)... ]
    end
    Vil     = bathparamPHev[DimBathHalf+1:end]
    Vil     = reshape( Vil, DimOrb, : )
    Vil     = Vil[:,IndBathHyb]
    return [ GetDeltaHybDiscFromFlat( ebathl, Vil, w, DimOrb, DimBath ) for w in wgrid ]
end


@inline function GetCostFromFlat( bathparam, systemparam ) 
    HybwOrig    = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    HybwParam   = GetDeltaHybDiscGridFromFlat( bathparam, wgrid, dimorb, dimbath )
    return GetCost( HybwOrig, HybwParam )
end

@inline function GetCostFromFlat( bathparam::Matrix{Float64}, systemparam::Tuple{Vector{Matrix{ComplexF64}}, Vector{ComplexF64}, Int64, Int64} ) 
    HybwOrig    = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    HybwParam   = GetDeltaHybDiscGridFromFlat( bathparam, wgrid, dimorb, dimbath )
    return GetCost( HybwOrig, HybwParam )
end

@inline function GetCostFromFlatPH( bathparamPH, systemparam )
    HybwOrig    = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    HybwParam   = GetDeltaHybDiscGridFromFlatPH( bathparamPH, wgrid, dimorb, dimbath )
    return GetCost( HybwOrig, HybwParam )
end

@inline function GetCostFromFlatPHev( bathparamPH, systemparam )
    HybwOrig    = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    HybwParam   = GetDeltaHybDiscGridFromFlatPHev( bathparamPH, wgrid, dimorb, dimbath )
    return GetCost( HybwOrig, HybwParam )
end


@inline function GetCostGlocinvFromFlat( bathparam, systemparam ) 
    GinvwOrig   = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    GinvwParam  = GetDeltaHybDiscGridFromFlat( bathparam, wgrid, dimorb, dimbath )
    return GetCost( GinvwOrig, GinvwParam )
end

@inline function GetCostGlocinvFromFlatPH( bathparamPH, systemparam )
    GinvwOrig   = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    GinvwParam  = GetDeltaHybDiscGridFromFlatPH( bathparamPH, wgrid, dimorb, dimbath )
    return GetCost( GinvwOrig, GinvwParam )
end


@inline function GetGreenFromHCFcoeff( H, AdagVec, AVec, Freq, E0 ; NvecLanc=800, NCFcoeff=0 )
    return GetGreenFromHCFcoeffFer( H, AdagVec, AVec, Freq, E0 ; NvecLanc=NvecLanc, NCFcoeff=NCFcoeff )
end

@inline function GetGreenFromHCFcoeffFer( H, AdagVec, AVec, Freq, E0 ; NvecLanc=800, NCFcoeff=0 )
    G_2 = GetGreenFromHCFcoeffAdag( H, AdagVec, Freq, E0 ; NvecLanc=NvecLanc, NCFcoeff=NCFcoeff )
    G_1 = GetGreenFromHCFcoeffA(    H, AVec,    Freq, E0 ; NvecLanc=NvecLanc, NCFcoeff=NCFcoeff )
    return G_1 + G_2 
end

@inline function GetGreenFromHCFcoeffFerBos( H, AdagVec, AVec, Freq, E0 ; NvecLanc=800, NCFcoeff=0, SgnFerBos=1 )
    G_2 = GetGreenFromHCFcoeffAdag( H, AdagVec, Freq, E0 ; NvecLanc=NvecLanc, NCFcoeff=NCFcoeff )
    G_1 = GetGreenFromHCFcoeffA(    H, AVec,    Freq, E0 ; NvecLanc=NvecLanc, NCFcoeff=NCFcoeff )
    return G_1 + SgnFerBos * G_2 
end

function GetGreenFromHCFcoeffAdag( H, AdagVec, Freq, E0 ; NvecLanc=800, NCFcoeff=0 )
    AbsAdagVec  = norm(AdagVec)^2
    TriMat, Vecs = lanczos_algorithm( H, normalize(AdagVec), NvecLanc)

    dimM        = size(TriMat)[1]
    NCFcoeff    = NCFcoeff<1 ? dimM : NCFcoeff
    z           = Freq .+ E0  # Absorption
    GAdag       = [ continued_fraction(-TriMat, NCFcoeff, w ) for w in z ]  ### inv(z - H)
    return GAdag * AbsAdagVec
end

function GetGreenFromHCFcoeffA( H, AVec, Freq, E0 ; NvecLanc=800, NCFcoeff=0 )
    AbsAVec     = norm(AVec)^2
    TriMat, Vecs = lanczos_algorithm( H, normalize(AVec), NvecLanc)

    dimM        = size(TriMat)[1]
    NCFcoeff    = NCFcoeff<1 ? dimM : NCFcoeff
    z           = Freq .- E0  # Emission
    GA          = [ continued_fraction(TriMat, NCFcoeff, w ) for w in z ]  ### inv(z + H)
    return GA * AbsAVec

end

function GetGreenFromHCFLoadedAni( CFData, Freq, E0 ; NCFcoeff=0 )
    AbsAVec     = CFData[1]
    TriMatDiag  = CFData[2]
    TriMatUpper = CFData[3]
    TriDiag     = Tridiagonal( conj(TriMatUpper), TriMatDiag, TriMatUpper )

    dimM        = size(TriMat)[1]
    TriMat      = Matrix{ComplexF64}(undef,dimM,dimM)
    TriMat      = zero(TriMat) + TriDiag

    NCFcoeff    = NCFcoeff<1 ? dimM : NCFcoeff
    z           = Freq .- E0  # Emission
    GA          = [ continued_fraction(TriMat, NCFcoeff, w ) for w in z ]  ### inv(z + H)
    return GA * AbsAVec
end

function GetGreenFromHCFLoadedCre( CFData, Freq, E0 ; NCFcoeff=0 )
    AbsAVec     = CFData[1]
    TriMatDiag  = CFData[2]
    TriMatUpper = CFData[3]
    TriDiag     = Tridiagonal( conj(TriMatUpper), TriMatDiag, TriMatUpper )

    dimM        = size(TriMat)[1]
    TriMat      = Matrix{ComplexF64}(undef,dimM,dimM)
    TriMat      = zero(TriMat) + TriDiag

    NCFcoeff    = NCFcoeff<1 ? dimM : NCFcoeff
    z           = Freq .+ E0  # Absorption
    GA          = [ continued_fraction(-TriMat, NCFcoeff, w ) for w in z ]  ### inv(z - H)
    return GA * AbsAVec
end

@inline function GetGreenLocalFromHk( Hk, weight_karr, Freq , chem )
    # @show size([ inv( -Hk_i + I * Freq ) for Hk_i in Hk ])
    # @show size( weight_karr )

    nk  = length(Hk)
    # @show nk

    zk      = [ zero(Hk[1]) + I * Freq + chem  for ik=1:nk ]
    gzkinv  = -Hk + zk
    M       = weight_karr .* inv.(gzkinv)
    # @show size(M)
    # @show typeof(M)
    # @show M
    # writedlm(stdout,M )
    # @show size(sum(M, dims=1))
    # @show sum(M, dims=1)
    # @show size(sum(M, dims=1)[1])
    # @show sum(M, dims=1)[1]
    # exit(1)
    return sum(M, dims=1)[1]
end

@inline function GetGreenLocalFromHkSelf( Hk, weight_karr, Freq, self, chem )
    nk      = length(Hk)
    zk      = [ zero(Hk[1]) + I * Freq - self + chem for ik=1:nk]
    gzkinv  = -Hk + zk
    M       = weight_karr .* inv.(gzkinv)
    return sum(M, dims=1)[1]
end

@inline function GetGreenLatticeDiagFromHkSelf( Hk, Freq, self, chem )
    nk  = length(Hk)
    M   = [ diag(inv( -Hk[ik] + I * Freq - self + chem)) for ik in 1:nk ]
    return M
end

@inline function GkSegment( (Hki, wki), freqi, chem )
    return inv( -Hki + I*freqi + chem ) * wki
end
@inline function GkSegment( Hki::Matrix{ComplexF64}, wki::Float64, freqi::ComplexF64, chem )
    return inv( -Hki + I*freqi + chem ) * wki
end

@inline function GetGreenLocalFromHkThread( Hk, weight_karr, Freq, chem )
    nk  = length(Hk)
    nth = Threads.nthreads()

    M   = [ zero(Hk[1]) for ith in 1:nth ]
    chunks = collect(Iterators.partition( 1:nk, div(nk,nth) ))

    Threads.@threads for ith in 1:nth
        chunk   = chunks[ith]
        wkchunk = weight_karr[chunk]
        Hkchunk = Hk[chunk]
        nkchunk = length(chunk)
        zk      = [ zero(Hkchunk[1]) + I * Freq - chem for ik=1:nkchunk ]
        Mchunk  = wkchunk .* ( inv.(Hkchunk + zk) )
        # Mchunk  = GkSegment.( collect(Iterators.zip( Hkchunk, wkchunk )), Freq, chem )
        # Mchunk  = [ GkSegment( hk, wk, Freq, chem ) for (hk,wk) in collect(Iterators.zip( Hkchunk, wkchunk )) ]
        
        # @show size(Mchunk)
        # @show typeof(Mchunk)
        # @show size(sum(Mchunk,dims=1))
        # @show typeof(sum(Mchunk,dims=1))
        M[ith]  = sum(Mchunk, dims=1)[1]
    end
    return sum(M, dims=1)[1]
end

@inline function GetGreenLocalFromHkGrid( Hk, weight_karr, Freq, chem)
    return pmap( freq_i -> GetGreenLocalFromHk(Hk,weight_karr,freq_i, chem), Freq )
    # return [ GetGreenLocalFromHk(Hk,weight_karr,freq_i, chem) for freq_i in Freq ]
    # return pmap( freq_i -> GetGreenLocalFromHkThread(Hk,weight_karr,freq_i, chem), Freq )
    # return [ GetGreenLocalFromHkThread(Hk,weight_karr,freq_i, chem) for freq_i in Freq ]
end

@inline function GetGreenLocalFromHkSelfGrid( Hk, weight_karr, Freq, Self_freq, chem)
    return pmap( ((freq_i,self_i),) -> GetGreenLocalFromHkSelf(Hk,weight_karr,freq_i,self_i,chem), Iterators.zip(Freq,Self_freq) )
end

@inline function GetGreenLocalFromHybGrid( Hyb, Freq, Eorb )
    nfreq   = length(Freq) 
    return pmap( ifreq -> inv( -Hyb[ifreq] + I * Freq[ifreq] - Eorb ), 1:nfreq )
    # return [ inv( -Hyb[ifreq] + I * Freq[ifreq] - Eorb ) for ifreq in 1:nfreq ]
end

@inline function GetGreenLatticeDiagFromHkSelfGrid( Hk, Freq, Self_freq, chem_int)
    return pmap( ((freq_i,self_i),) -> GetGreenLatticeDiagFromHkSelf(Hk,freq_i,self_i,chem_int), Iterators.zip(Freq,Self_freq) )
end

@inline function GetGreenLatticeTrFromDiag( Gk_diag )
    nk  = length(Gk_diag)
    G   = [ sum(Gk_diag_i) for Gk_diag_i in Gk_diag ]
    return G
end

@inline function GetGreenLatticeTrFromDiagGrid( Gkw_diag )
    return pmap( Gk_diag -> GetGreenLatticeTrFromDiag( Gk_diag ), Gkw_diag )
end


@inline function GetHybFromGreenLocalFromHk( Hk, weight_karr, Freq, chem, Eorb )
    nk  = length(Hk)
    G   = [ weight_karr[ik] * inv( -Hk[ik] + I * Freq + chem) for ik in 1:nk ]
    Gloc= sum(G, dims=1)[1]
    Hyb = - inv(Gloc) + I * Freq - Eorb + chem
    return Hyb
end

@inline function GetHybFromGreenLocalFromHkGrid( Hk, weight_karr, Freq, chem, Eorb )
    return pmap( freq_i -> GetHybFromGreenLocalFromHk(Hk,weight_karr,freq_i,chem,Eorb), Freq)
end

@inline function GetHybFromGreenLocalGrid( gloc, Freq, Eorb )
    nfreq   = length(Freq) 
    return pmap( indfreq -> -inv(gloc[indfreq]) + I * Freq[indfreq] - Eorb, 1:nfreq )
end

@inline function invMgrid( Ggrid::Vector{Matrix{ComplexF64}} )
    return [ inv(Gcomp) for Gcomp in Ggrid ]
end

function ConstrainTrev!( g_grid, pair_ind ) 
    ind_up, ind_dn  = pair_ind
    for g in g_grid
        g[ind_dn,ind_dn]    = g[ind_up, ind_up]
    end
end

function ConstrainTrevAllSz!( g_grid, norb ) 
    ind_ups = collect(1:2:norb)
    ind_dns = collect(2:2:norb)
    for pair_i in Iterators.zip( ind_ups, ind_dns )
        ConstrainTrev!( g_grid, pair_i )
    end
end

function ConstrainDiag!( g_grid, pair_ind ) 
    ind_up, ind_dn  = pair_ind
    for g in g_grid
        g[ind_up,ind_dn]    = zero(g[ind_up, ind_dn])
    end
end

function ConstrainDiagAll!( g_grid, norb ) 
    ind_pairs = [ (i,j) for i in 1:norb for j in (i+1):norb ]
    for pair_i in ind_pairs
        @show pair_i
        ConstrainDiag!( g_grid, pair_i )
        ConstrainDiag!( g_grid, reverse(pair_i) )
    end
end

@inline function GetDensityMatrixFromGreenImFreq( giw, imfreq ; epsilon=1e-6 )
    beta    = pi / imag(imfreq[1])
    return sum( [ giw_i * exp(imfreq_i*epsilon) + conj(giw_i) * exp(-imfreq_i*epsilon) for (giw_i,imfreq_i) in Iterators.zip(giw,imfreq) ], dims=1)[1] / beta + 0.5*I
end

@inline function GetDensityInfoFromGreenImFreq( giw, imfreq ; epsilon=1e-6 )
    rho     = GetDensityMatrixFromGreenImFreq( giw, imfreq ; epsilon=epsilon ) 
    n_arr   = real(diag(rho))
    ntot    = sum(n_arr)
    return [ntot, n_arr... ]
end

function GetGimpDegen( GSSectorInfo, nspinorb, ntot, opcavec, opccaavec, FreqGrid )
    GimpDegen   = Matrix{ComplexF64}[ zeros(ComplexF64,nspinorb,nspinorb) for freq in FreqGrid ]
    for (iGSSectorInfo, GSSectorInfo_i) in Iterators.enumerate(GSSectorInfo)
        Wboltz  = GSSectorInfo_i[5]
        Gimp    = GetGreenImpurityFromGS(  nspinorb, GSSectorInfo_i, ntot, opcavec, opccaavec, FreqGrid ) * Wboltz
        GimpDegen   .+= Gimp
    end
    return GimpDegen
end
