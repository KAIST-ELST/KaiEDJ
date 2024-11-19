export OpCA
export OpCCAA
export ArrOpCA
export ArrOpCCAA
export act_vij
export act_vij_sign
export act_vijkl
export act_vijkl_sign
export checkact_vij_sign
export checkact_vijkl_sign
export OpCABinst
export OpCAIndBinst
export OpCCAAIndBinst
export OpCAWFadd
export OpC
export OpA
export OpCWFadd
export OpAWFadd

struct OpA
    v::ComplexF64
    i::Array{Int,1}
end

struct OpC
    v::ComplexF64
    i::Array{Int,1}
end

struct OpCA
    v::ComplexF64
    i::Array{Int,1}
end

struct ArrOpCA
    Vij::Array{OpCA,1}
    ArrOpCA()   = new()
end

struct OpCCAA
    v::ComplexF64
    i::Array{Int,1}
end

struct ArrOpCCAA
    Vijkl::Array{OpCCAA,1}
end

function OpCABinst( op::OpCA, stbin::Int)
    v       = op.v
    i, j    = op.i 
    if CheckOpCABin( stbin, i, j ) 
        FerSign, stbin_new  = act_vij_sign(stbin, i, j)
        return v*FerSign, stbin_new
    else 
        return 0, 0
    end
end

function OpCAIndBinst( op::OpCA, istbin::Int, hashf::Vector{Int64}, hashfinv::Vector{Vector{Int64}} ; outputlevel=0 )
    v       = op.v
    i, j    = op.i 
    stbin   = hashf[istbin]-1
    FerSign, stbin_new  = checkact_vij_sign(stbin, i, j)
    if stbin_new > 0
        # if outputlevel > 0
        #     @show hashf
        #     @show hashfinv
        #     @show stbin
        #     @show stbin_new
        #     @show getbittaili(stbin)
        #     @show getbittaili(stbin_new)
        # end
        jstbin              = hashfinv[stbin_new+1][2]
        return FerSign, jstbin
    else 
        return 0, 0
    end
end

function OpCAIndBinst( op::OpCA, istbin::Int, hashf::Vector{Int64}, hashfinv::Vector{Int64} ; outputlevel=0 )
    v       = op.v
    i, j    = op.i 
    stbin   = hashf[istbin]-1
    FerSign, stbin_new  = checkact_vij_sign(stbin, i, j)
    if stbin_new > 0
        # if outputlevel > 0
        #     @show hashf
        #     @show hashfinv
        #     @show stbin
        #     @show stbin_new
        #     @show getbittaili(stbin)
        #     @show getbittaili(stbin_new)
        # end
        jstbin              = hashfinv[stbin_new+1]
        return FerSign, jstbin
    else 
        return 0, 0
    end
end

function OpCCAAIndBinst( op::OpCCAA, istbin::Int, hashf::Vector{Int64}, hashfinv::Vector{Vector{Int64}} ; outputlevel=0 )
    v           = op.v
    i, j, k, l  = op.i 
    stbin   = hashf[istbin]-1
    FerSign, stbin_new  = checkact_vijkl_sign(stbin, i, j, k, l)
    # if outputlevel > 0
    #     @show hashf
    #     @show hashfinv
    #     @show stbin
    #     @show stbin_new
    #     @show getbittaili(stbin)
    #     @show getbittaili(stbin_new)
    # end
    if stbin_new > 0
        jstbin              = hashfinv[stbin_new+1][2]
        return FerSign, jstbin
    else 
        return 0, 0
    end
end

function OpCCAAIndBinst( op::OpCCAA, istbin::Int, hashf::Vector{Int64}, hashfinv::Vector{Int64} ; outputlevel=0 )
    v           = op.v
    i, j, k, l  = op.i 
    stbin   = hashf[istbin]-1
    FerSign, stbin_new  = checkact_vijkl_sign(stbin, i, j, k, l)
    # if outputlevel > 0
    #     @show hashf
    #     @show hashfinv
    #     @show stbin
    #     @show stbin_new
    #     @show getbittaili(stbin)
    #     @show getbittaili(stbin_new)
    # end
    if stbin_new > 0
        jstbin              = hashfinv[stbin_new+1]
        return FerSign, jstbin
    else 
        return 0, 0
    end
end

function OpCAWFadd( op::OpCA, wf::WF, wfnew::WF ; outputlevel=0 , tol=1e-15) 
    v       = op.v
    i, j    = op.i 
    dim     = wf.Dim
    for indBin in 1:dim
        stbin   = wf.BasisBinstate[indBin]-1
        coeff   = wf.Probamp[indBin]
        if abs(coeff) > tol
            FerSign, stbin_new  = checkact_vij_sign( stbin, i, j ) 
            if (stbin_new > 0)&( iszero(wf.Noccubasis-wfnew.Noccubasis) )   ## Number-conservation 
                indBin_new          = wfnew.HashfInv[stbin_new+1]
                wfnew.Probamp[indBin_new]  += coeff * v * FerSign
                # if outputlevel > 0
                #     println( "bin1=$(indBin) bin2=$(indBin_new) [$(getbittaili(stbin))] [$(getbittaili(stbin_new))] ", format("{:.8f} {:.8f} ",abs(wfnew.Probamp[indBin_new]), abs(coeff)), ", $v , $FerSign" )
                # end
            end
        end
    end
    ## Probably, this cdagi cj matrix can be stored in WF and so CCAA is.
    ## For instance, vij_mat[ st, ij ]  = st_new 
    ## For instance, vijkl_mat[ st, ijkl ]  = st_new 
end

function OpCWFadd( op::OpC, wf::WF, wfnew::WF ) 
    v       = op.v
    i       = op.i[1] 
    dim     = wf.Dim
    for indBin in 1:dim
        stbin   = wf.BasisBinstate[indBin]-1
        coeff   = wf.Probamp[indBin]
        FerSign, stbin_new  = checkact_cdagi_sign( stbin, i )
        if (stbin_new > 0)&( (wfnew.Noccubasis<0) | iszero(wf.Noccubasis+1-wfnew.Noccubasis) )  ## Non-conserved-number basis or Number-increase 
            indBin_new          = wfnew.HashfInv[stbin_new+1]
            wfnew.Probamp[indBin_new]  += coeff * v * FerSign
            # println( "$(indBin) $(indBin_new) [$(getbittaili(stbin))] [$(getbittaili(stbin_new))] $(wfnew.Probamp[indBin_new]) ; $(coeff), $v , $FerSign" )
        end
    end
    ## Probably, this cdagi cj matrix can be stored in WF and so CCAA is.
    ## For instance, vij_mat[ st, ij ]  = st_new 
    ## For instance, vijkl_mat[ st, ijkl ]  = st_new 
end

function OpAWFadd( op::OpA, wf::WF, wfnew::WF ) 
    v       = op.v
    i       = op.i[1] 
    dim     = wf.Dim
    for indBin in 1:dim
        stbin   = wf.BasisBinstate[indBin]-1
        coeff   = wf.Probamp[indBin]
        FerSign, stbin_new  = checkact_ci_sign( stbin, i )
        if (stbin_new > 0)&( (wfnew.Noccubasis<0) | iszero(wf.Noccubasis-1-wfnew.Noccubasis) )  ## Non-conserved-number basis or Number-decrease
            indBin_new          = wfnew.HashfInv[stbin_new+1] # HashfInv[stbin_new+1][2]
            wfnew.Probamp[indBin_new]  += coeff * v * FerSign
            # println( "$(indBin) $(indBin_new) [$(getbittaili(stbin))] [$(getbittaili(stbin_new))] $(wfnew.Probamp[indBin_new]) ; $(coeff), $v , $FerSign" )
        end
    end
    ## Probably, this cdagi cj matrix can be stored in WF and so CCAA is.
    ## For instance, vij_mat[ st, ij ]  = st_new 
    ## For instance, vijkl_mat[ st, ijkl ]  = st_new 
end

@inline CheckOpCABin( st, i, j )            = isemp_zbase(st,i-1) & isocc_zbase(st,j-1)
@inline CheckOpCCAABin( st, i, j, k, l )    = isemp_zbase(st,i-1) & isemp_zbase(st,j-1) & isocc_zbase(st,k-1) & isocc_zbase(st,l-1)


@inline act_vij( st , i , j ) = cdagi_zbase( ci_zbase(deepcopy(st),j-1), i-1 )
# function act_vij( st , i , j )
#     st_new  = deepcopy(st)
#     st_new  = ci_zbase( st_new, j )
#     st_new  = cdagi_zbase( st_new, i )
#     return st_new 
# end

function act_vij_sign( st , i , j )
    ii, jj      = [ i-1, j-1 ]
    st_new      = deepcopy(st)
    nferswap    = NFerSwap(st_new,jj)
    # println( "[$(jj)=$(isocc_zbase(st_new,jj))] $(getbittaili(st_new)) $(nferswap)" )
    st_new      = ci_zbase( st_new, jj )
    nferswap    += NFerSwap(st_new,ii)
    # println( "[$(ii)=$(isemp_zbase(st_new,ii))] $(getbittaili(st_new)) $(nferswap)" )
    st_new      = cdagi_zbase( st_new, ii )
    # println(" returning : $( (-1)^nferswap) , $(st_new) " )
    return (-1)^nferswap, st_new 
end


function checkact_vij_sign( st::Int , i::Int , j::Int )
    ii = i-1
    jj = j-1
    st_new      = deepcopy(st)
    nferswap    = 0
    # println( "[$(jj)=$(isocc_zbase(st_new,jj))] $(getbittaili(st_new)) $(nferswap)" )

    if isocc_zbase( st_new, jj ) 
        nferswap    = NFerSwap(st_new,jj)
        st_new      = ci_zbase( st_new, jj )
        # println( "[$(ii)=$(isemp_zbase(st_new,ii))] $(getbittaili(st_new)) $(nferswap)" )

        if isemp_zbase( st_new, ii )
            nferswap    += NFerSwap(st_new,ii)
            st_new      = cdagi_zbase( st_new, ii )
            # println( "[completed] $(getbittaili(st_new)) $(nferswap)" )
            return (-1)^nferswap, st_new
        else
            return 0, 0
        end
    else
        return 0, 0
    end
end



@inline act_vijkl( st , i , j, k, l ) = cdagi_zbase( cdagi_zbase( ci_zbase( ci_zbase(deepcopy(st),l-1), k-1 ), j-1), i-1 )
# function act_vijkl( st , i , j, k, l )
#     st_new  = deepcopy(st)
#     st_new  = ci_zbase( st_new, l )
#     st_new  = ci_zbase( st_new, k )
#     st_new  = cdagi_zbase( st_new, j )
#     st_new  = cdagi_zbase( st_new, i )
#     return st_new 
# end
function act_vijkl_sign( st , i , j, k, l )
    ii, jj, kk, ll  = [ i-1, j-1, k-1, l-1 ]
    st_new      = deepcopy(st)
    nferswap    = NFerSwap(st_new,ll)
    st_new      = ci_zbase( st_new, ll )
    nferswap    += NFerSwap(st_new,kk)
    st_new      = ci_zbase( st_new, kk )
    nferswap    += NFerSwap(st_new,jj)
    st_new      = cdagi_zbase( st_new, jj )
    nferswap    += NFerSwap(st_new,ii)
    st_new      = cdagi_zbase( st_new, ii )
    return (-1)^nferswap, st_new 
end
function checkact_vijkl_sign( st::Int , i::Int , j::Int, k::Int, l::Int )
    ii = i-1
    jj = j-1
    kk = k-1
    ll = l-1
    st_new  = deepcopy(st)
    nferswap = 0
    # println( "[$(ll)=$(isocc_zbase(st_new,ll))] $(getbittaili(st_new)) $(nferswap)" )

    if isocc_zbase( st_new, ll ) 
        nferswap    = NFerSwap(st_new,ll)
        st_new      = ci_zbase( st_new, ll )
        # println( "[$(kk)=$(isocc_zbase(st_new,kk))] $(getbittaili(st_new)) $(nferswap)" )

        if isocc_zbase( st_new, kk )
            nferswap    += NFerSwap(st_new,kk)
            st_new      = ci_zbase( st_new, kk )
            # println( "[$(jj)=$(isemp_zbase(st_new,jj))] $(getbittaili(st_new)) $(nferswap)" )

            if isemp_zbase( st_new, jj )
                nferswap    += NFerSwap(st_new,jj)
                st_new      = cdagi_zbase( st_new, jj )
                # println( "[$(ii)=$(isemp_zbase(st_new,ii))] $(getbittaili(st_new)) $(nferswap)" )

                if isemp_zbase( st_new, ii )
                    nferswap    += NFerSwap(st_new,ii)
                    st_new      = cdagi_zbase( st_new, ii )
                    # println( "[completed] $(getbittaili(st_new)) $(nferswap)" )
                    return (-1)^nferswap, st_new
                else
                    return 0, 0
                end
            else
                return 0, 0
            end
        else
            return 0, 0
        end
    else 
        return 0, 0
    end
end


function checkact_cdagi_sign( st , i )
    ii  = i-1
    st_new      = deepcopy(st)
    nferswap    = 0
    # println( "[$(ii)=$(isemp_zbase(st_new,ii))] $(getbittaili(st_new)) $(nferswap)" )
    if isemp_zbase( st_new, ii )
        nferswap    = NFerSwap(st_new,ii)
        st_new      = cdagi_zbase( st_new, ii )
        # println( "[completed] $(getbittaili(st_new)) $(nferswap)" )
        return (-1)^nferswap, st_new
    else
        return 0, 0
    end
end

function checkact_ci_sign( st , i )
    ii  = i-1
    st_new      = deepcopy(st)
    nferswap    = 0
    # println( "[$(ii)=$(isocc_zbase(st_new,ii))] $(getbittaili(st_new)) $(nferswap)" )
    if isocc_zbase( st_new, ii )
        nferswap    = NFerSwap(st_new,ii)
        st_new      = ci_zbase( st_new, ii )
        # println( "[completed] $(getbittaili(st_new)) $(nferswap)" )
        return (-1)^nferswap, st_new
    else
        return 0, 0
    end
end

