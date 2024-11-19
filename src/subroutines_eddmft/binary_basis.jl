using Formatting

export cdagi_zbase
export ci_zbase
export isocc_zbase
export isemp_zbase
export NFerSwap
export OpNtot
export getbittaili

Norb     = 1
Nspinorb    = 2*Norb

Nsite   = 3
Ntotalbasis  = Nspinorb^Nsite

stateArrInt = collect(1:Ntotalbasis)
stateArrOcc = Vector{Int64}

# @inline cdagi_zbase(st,i) = st|(1<<i)
# @inline ci_zbase(st,i)    = st&(~(1<<i))
# @inline isocc_zbase(st,i) = Bool((st>>i)&1)
# @inline isemp_zbase(st,i) = Bool(((~st)>>i)&1)

@inline cdagi_zbase(st::Int,i::Int) = st|(1<<i)
@inline ci_zbase(st::Int,i::Int)    = st&(~(1<<i))
@inline isocc_zbase(st::Int,i::Int) = Bool((st>>i)&1)
@inline isemp_zbase(st::Int,i::Int) = Bool(((~st)>>i)&1)

@inline function NFerSwap( st::Int, iorb::Int ) 
    st_clear = (st >> iorb) << iorb
    return count_ones(st - st_clear)
end

# function NFerSwap( st, iorb ) 
#     nswap   = 0 
#     for jorb in 0:iorb-1
#         # println( "          ", jorb, isocc_zbase(st,jorb) )
#         if isocc_zbase(st,jorb)
#             nswap += 1
#         end
#     end
#     return nswap
# end


@inline OpNtot(st)  = count_ones(st)

# function OpNtot( st, norb )
#     Ntot   = 0 
#     for jorb in 0:norb-1
#         # println( "          ", jorb, isocc_zbase(st,jorb) )
#         if isocc_zbase(st,jorb)
#             Ntot += 1
#         end
#     end
#     return Ntot
# end


@inline getbittaili(binst)  = format("{}",bitstring(binst)[end-16:end])

## TEST EXAMPLE  ##
# println( "(x, i), format('{:s} {:s} {:s}',getbittaili(x), getbittaili(cdagi_zbase(x,i)), getbittaili(ci_zbase(x,i)))" )
# for x in 0:3
#     println("")
#     for i in 0:3
#         # println( (x, i), format("{:s} {:s} {:s}",getbittaili(x), getbittaili(cdagi_zbase(x,i)), getbittaili(ci_zbase(x,i))) )
#         cdx = cdagi_zbase(x,i)
#         cx  = ci_zbase(x,i)
#         println( ("state $x", "$i-th orbital") )
#         println( format("  |{:s}> {:s} {:s} ; nswap={:d} ",getbittaili(x),   isocc_zbase(x,i)  , isemp_zbase(x,i)  , NFerSwap(x,i)   ))
#         println( format("  |{:s}> {:s} {:s} ; nswap={:d} ",getbittaili(cdx), isocc_zbase(cdx,i), isemp_zbase(cdx,i), NFerSwap(cdx,i) ))
#         println( format("  |{:s}> {:s} {:s} ; nswap={:d} ",getbittaili(cx),  isocc_zbase(cx,i) , isemp_zbase(cx,i) , NFerSwap(cx,i)  ))
#     end
# end
