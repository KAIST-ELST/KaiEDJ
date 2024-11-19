export GetOpSz
export GetOpUSlaterKanamori
export GetOpBathParam
export GetOpChemPot

function GetOpSz( ntot ) 
    OpSz        = OpCA[]
    ntotHalf    = div(ntot,2)
    @assert ntot%2 < 1 
    for i in 1:ntotHalf
        push!( OpSz, OpCA(  1. / 2 , [ 2*i-1, 2*i-1 ] ) )
        push!( OpSz, OpCA( -1. / 2 , [ 2*i  , 2*i   ] ) )
    end
    return OpSz
end

@inline GetOpSzOrb( norb )  = GetOpSz( norb )


function GetOpUSlaterKanamori( ; U=1, JHund=0.167, norb=3 )
    Uijkl   = OpCCAA[]
    for i in 1:norb 
        push!( Uijkl, OpCCAA( U,             [2*i-1, 2*i  , 2*i,     2*i-1] ) )         ## t2g : intra-orbital density-density  ##
        for j in 1:norb 
            if i!=j
                push!( Uijkl, OpCCAA( U-2.0*JHund,   [2*i-1, 2*j  , 2*j,     2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i != j , opposite-spin
                push!( Uijkl, OpCCAA( -JHund,        [2*i-1, 2*j,   2*j-1,   2*i  ] ) ) ## t2g : pair-exchange                  ## # i != j , up|dn => dn|up
                push!( Uijkl, OpCCAA(  JHund,        [2*i-1, 2*i,   2*j  ,   2*j-1] ) ) ## t2g : pair-hopping                   ## # i != j , i,up|i,dn => j,up|j,dn
            end
            if j<i
                push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i-1, 2*j-1, 2*j-1,   2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i <  j , up-spin
                push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i  , 2*j  , 2*j  ,   2*i  ] ) ) ## t2g : inter-orbital density-density  ## # i <  j , dn-spin
            end
        end
    end
    return Uijkl 
end

function GetOpBathParam( ebathl, Vil, IndHamilBath ) 
    norb, nbath = size( Vil )
    Op_BathAll  = OpCA[]
    for j in 1:nbath
        jbath   = IndHamilBath( j )
        for i in 1:norb
            if abs(Vil[i,j])>1e-15 
                push!( Op_BathAll, OpCA(      Vil[i,j],  [i,     jbath] ) )
                push!( Op_BathAll, OpCA( conj(Vil[i,j]), [jbath, i    ] ) )
            end
        end
        if abs(ebathl[j])>1e-15 
            push!( Op_BathAll, OpCA( ebathl[j], [jbath,jbath] ) )
        end
    end
    return Op_BathAll
end

function GetOpChemPot( val::ComplexF64, rangeorb::UnitRange{Int64} )
    tij = OpCA[]
    for i in rangeorb
        push!( tij, OpCA( -val , [ i, i ] ) )  
    end
    return tij
end

@inline function GetOpChemPot( val::ComplexF64, norb::Int )
    return GetOpChemPot( val, 1:norb )
end

@inline function GetOpMatrix( M::Matrix{ComplexF64}, norb::Int )
    return [[ OpCA( M[i,j], [i,j] ) for i=1:norb, j=1:norb ]...]
end

@inline function GetOpChemPot( M::Matrix{ComplexF64}, norb::Int )
    return GetOpMatrix( -M, norb ) 
end

@inline GetOpChemPot( val::Float64, rangeorb::UnitRange{Int64} )    = GetOpChemPot( convert(ComplexF64,val), rangeorb )
@inline GetOpChemPot( val::Float64, norb::Int )                     = GetOpChemPot( convert(ComplexF64,val), norb )
@inline GetOpMatrix( M::Matrix{Float64}, norb::Int )    = GetOpMatrix( M .+ 0.0*im, norb )
@inline GetOpChemPot( M::Matrix{Float64}, norb::Int )   = GetOpChemPot( M .+ 0.0*im, norb )
