export GetNtotNdiagFromDM
export GetSzTotSzdiagFromDM
export GetSzTotSzdiagFromArr
export GetObservables

function GetObservables(GSSectorInfo, nspinorb, ntot)
    DMDegen     = zeros(ComplexF64,nspinorb,nspinorb)
    nGS         = length(GSSectorInfo)
    obs_info    = [[] for i=1:nGS+1]
    DMEach      = [zero(DMDegen) for i=1:nGS]
    nImp_info   = [[] for i=1:nGS]
    for (iGSSectorInfo, GSSectorInfo_i) in Iterators.enumerate(GSSectorInfo)
        Wboltz  = GSSectorInfo_i[5]
        DM      = GetDensityMatrix( GSSectorInfo_i, ntot, nspinorb ) 
        DMEach[iGSSectorInfo]   = deepcopy(DM)
        DMDegen .+= DM * Wboltz
        push!( obs_info[iGSSectorInfo], iGSSectorInfo, Wboltz, GetNtotNdiagFromDM( DM ), GetSzTotSzdiagFromDM( DM ) )
        nImpArr = GetNtot( GSSectorInfo_i, ntot, ntot )
        push!( nImp_info[iGSSectorInfo], iGSSectorInfo, Wboltz, sum(nImpArr),  nImpArr )
    end
    push!( obs_info[nGS+1], "tot", "1", GetNtotNdiagFromDM( DMDegen ), GetSzTotSzdiagFromDM( DMDegen ) )

    return DMDegen, DMEach, obs_info, nImp_info
end

function GetNtot( GSSectorInfo_i, norbsys, nspinorb )
    hashfsec    = HashingIntToSecarrBin(;norb=norbsys )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    DimSec      = length(HashF)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )

    OpNtot  = [ OpCA(1.0, [ i , i ]) for i=1:nspinorb ]
    Ntotval = [ 0.0*im for i=1:nspinorb ]

    for i=1:nspinorb
        Ntotwf      = WF(DimSec, HashF, wfargs... )
        OpCAWFadd( OpNtot[i], wf, Ntotwf )
        Ntotval[i]  = dot( wf.Probamp, Ntotwf.Probamp )
    end
    return real(Ntotval)
end

function GetDensityMatrix( GSSectorInfo_i, norbsys, nspinorb )
    hashfsec    = HashingIntToSecarrBin(;norb=norbsys )
    hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec )
    hashfinv    = HashingInvSecarrBinToInt( hashfsec )
    wfargs      = [ hashfinv ]

    isector = GSSectorInfo_i[1][1]
    Eimp0   = GSSectorInfo_i[2]
    wboltz  = GSSectorInfo_i[3]
    gsamp   = GSSectorInfo_i[4]

    HashF       = hashfsec[isector]
    DimSec      = length(HashF)

    wf      = WF(DimSec, HashF, gsamp, wfargs... )

    OpRho       = [ OpCA(1.0, [ i , j ]) for i=1:nspinorb, j=1:nspinorb ]
    DensityM    = [ 0.0*im for i=1:nspinorb, j=1:nspinorb ]

    for i=1:nspinorb, j=1:nspinorb 
        rhowf   = WF(DimSec, HashF, wfargs... )
        OpCAWFadd( OpRho[i,j], wf, rhowf )
        DensityM[i,j]   = dot( wf.Probamp, rhowf.Probamp )
    end
    return DensityM
end

function GetNtotNdiagFromDM( DM::Matrix{ComplexF64} )
    ndiag   = real( diag(DM) )
    ntot    = sum(ndiag)
    return [ ntot ndiag... ]
end

function GetSzTotSzdiagFromDM( DM::Matrix{ComplexF64} )
    ndiag   = real( diag(DM) )
    szmatdiag   = [ 0.5*(-1)^(i-1) for i = 1:size(DM)[1] ]
    szdiag  = ndiag .* szmatdiag
    sztot   = sum(szdiag)
    return [ sztot szdiag... ]
end

function GetSzTotSzdiagFromArr( ndiag::Vector{ComplexF64} )
    szmatdiag   = [ 0.5*(-1)^(i-1) for i = 1:size(DM)[1] ]
    szdiag  = ndiag .* szmatdiag
    sztot   = sum(szdiag)
    return [ sztot szdiag... ]
end
