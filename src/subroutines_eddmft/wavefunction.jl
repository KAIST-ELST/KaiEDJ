export WF
export HashingIntToSecarrBin
export HashingIntToBinall
export HashingIntToBinallOrdered
export HashingInvSecarrBinToAllNparInt
export HashingInvSecarrBinToInt
export HashingInvBinallToInt
export ShowHashfsec
export ShowHashfsecDim
export ShowHashfInv
export ShowWFBasis
export ShowWFProb
export WFInner
export WFNorm
export WFCompare
export WFCompareArr
export WFCompareArrFtnbin
export WFCompareAbs
export WFCompareAbsArr
export WFCompareAbsArrFtnbin
export WFCompareAbsdiff
export ProbCompareArr
export ProbCompareAbsArr
export WFCopyDiffBasis
export ShowProbArrNonzero

mutable struct WF
    Dim::Int
    BasisBinstate::Vector{Int64}
    Probamp::Vector{ComplexF64}

    HashfInv::Vector{Int}
    IndSector::Int
    Noccubasis::Int

    WF( wf::WF )    = new( wf.Dim, wf.BasisBinstate, zeros(wf.Dim)*im, wf.HashfInv, wf.IndSector, wf.Noccubasis )
    WF( d::Int64 )  = new( d, collect(0:d-1), zeros(d)*im, [], -1, -1 )
    WF( d::Int64, bb::Vector{Int64} )                                   = new( d, bb, zeros(d)*im, [], -1, -1 )
    WF( d::Int64, bb::Vector{Int64}, bbInv::Vector{Int64} )             = new( d, bb, zeros(d)*im, bbInv, -1, -1 )
    WF( d::Int64, bb::Vector{Int64}, bbInv::Vector{Int64}, isec::Int )  = new( d, bb, zeros(d)*im, bbInv, isec, -1 )
    WF( d::Int64, bb::Vector{Int64}, bbInv::Vector{Int64}, isec::Int, noccu::Int )  = new( d, bb, zeros(d)*im, bbInv, isec, noccu )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{Float64}, bbInv::Vector{Int64} )                        = new( d, bb, pamp, bbInv, -1, -1 )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{Float64}, bbInv::Vector{Int64}, isec::Int )             = new( d, bb, pamp, bbInv, isec, -1 )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{Float64}, bbInv::Vector{Int64}, isec::Int, noccu::Int ) = new( d, bb, pamp, bbInv, isec, noccu )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{ComplexF64}, bbInv::Vector{Int64} )                        = new( d, bb, pamp, bbInv, -1, -1 )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{ComplexF64}, bbInv::Vector{Int64}, isec::Int )             = new( d, bb, pamp, bbInv, isec, -1 )
    WF( d::Int64, bb::Vector{Int64}, pamp::Vector{ComplexF64}, bbInv::Vector{Int64}, isec::Int, noccu::Int ) = new( d, bb, pamp, bbInv, isec, noccu )
end

@inline function WFInner( wf::WF, wf2::WF )
    return wf.Probamp' * wf2.Probamp
end

@inline function WFadd( wf::WF, wf2::WF )
    return wf.Probamp += wf2.Probamp
end

@inline function WFNorm( wf::WF )
    return norm(wf.Probamp)
end

function ShowWFProb( wf::WF )
    for i in 1:wf.Dim
        binst   = wf.BasisBinstate[i]-1
        println( "WF-basis[ $i ] : $(binst) [$(getbittaili(binst))] prob=$(abs(wf.Probamp[i]))" )
    end
end


function WFCopyDiffBasis( wf::WF, wf2::WF )
    for i in 1:wf.Dim
        binst_i = wf.BasisBinstate[i]-1
        j_wf2   = wf2.HashfInv[binst_i+1]
        wf2.Probamp[j_wf2]  = wf.Probamp[i]
    end
end

function WFCompare( wf::WF, wf2::WF )
    for i in 1:wf.Dim
        binst   = wf.BasisBinstate[i]-1
        binst2  = wf2.BasisBinstate[i]-1
        println( "WF-basis[ $i ] : $(binst+1) [$(getbittaili(binst))] [$(getbittaili(binst2))] prob1=$(wf.Probamp[i]) prob2=$(wf2.Probamp[i])" )
    end
end

function WFCompareArr( wfarr::Vector{WF} ; tol=1e-5)
    wfdum   = wfarr[1]
    for i in 1:wfdum.Dim
        binst   = wfdum.BasisBinstate[i]
        BoolAllzero = prod([ abs(wf.Probamp[i])<tol for wf in wfarr ])
        if !BoolAllzero
            probarr = [ abs(wf.Probamp[i])<tol ? format("{:10d}",0) :  format("{:10.3e} +i {:10.3e}", real(wf.Probamp[i]), imag(wf.Probamp[i])) for wf in wfarr ]
            println( "WF-basis[ $i ] : $(binst) [$(getbittaili(binst))] $(probarr...)" )
        end
    end
end

function WFCompareArrFtnbin( wfarr::Vector{WF}, ftnbin ; tol=1e-5)
    wfdum   = wfarr[1]
    for i in 1:wfdum.Dim
        binst   = wfdum.BasisBinstate[i]
        BoolAllzero = prod([ abs(wf.Probamp[i])<tol for wf in wfarr ])
        if !BoolAllzero
            probarr = [ abs(wf.Probamp[i])<tol ? format("{:10d}",0) :  format("{:10.3e} +i {:10.3e}", real(wf.Probamp[i]), imag(wf.Probamp[i])) for wf in wfarr ]
            println( format("WF-basis[ {:4d} ] : {:4d} [$(getbittaili(binst))][$(ftnbin(binst))] $(probarr...)", i, binst) )
        end
    end
end

function WFCompareAbs( wf::WF, wf2::WF )
    for i in 1:wf.Dim
        binst   = wf.BasisBinstate[i]
        binst2  = wf2.BasisBinstate[i]
        println( "WF-basis[ $i ] : $(binst) [$(getbittaili(binst))] [$(getbittaili(binst2))] prob1=$(abs(wf.Probamp[i])) prob2=$(abs(wf2.Probamp[i]))" )
    end
end

function WFCompareAbsArr( wfarr::Vector{WF} ; tol=1e-5)
    wfdum   = wfarr[1]
    for i in 1:wfdum.Dim
        binst   = wfdum.BasisBinstate[i]
        BoolAllzero = prod([ abs(wf.Probamp[i])<tol for wf in wfarr ])
        if !BoolAllzero
            probarr = [ abs(wf.Probamp[i])<tol ? format("{:10d}",0) :  format("{:10.3e}", abs(wf.Probamp[i])) for wf in wfarr ]
            println( "WF-basis[ $i ] : $(binst) [$(getbittaili(binst))] $(probarr...)" )
        end
    end
end

function WFCompareAbsArrFtnbin( wfarr::Vector{WF}, ftnbin ; tol=1e-5)
    wfdum   = wfarr[1]
    for i in 1:wfdum.Dim
        binst   = wfdum.BasisBinstate[i]
        BoolAllzero = prod([ abs(wf.Probamp[i])<tol for wf in wfarr ])
        if !BoolAllzero
            probarr = [ abs(wf.Probamp[i])<tol ? format("{:10d}",0) :  format("{:10.3e}", abs(wf.Probamp[i])) for wf in wfarr ]
            println( format("WF-basis[ {:4d} ] : {:4d} [$(getbittaili(binst))][$(ftnbin(binst))] $(probarr...)", i, binst) )
        end
    end
end

function ShowProbArrNonzero( wfprobarr::Vector ; tol=1e-5)
    for i in 1:length(wfprobarr[1])
        BoolAllzero = prod([ abs(wfprob[i])<tol for wfprob in wfprobarr ])
        if !BoolAllzero
            val_wfprobarr = [ abs(wfprob[i])<tol ? format("{:10d}",0) :  format("{:10.3e}", abs(wfprob[i])) for wfprob in wfprobarr ]
            println( "[ $i ] : $(val_wfprobarr...)" )
        end
    end
end

function ProbCompareArr( wfprobarr::Vector ; tol=1e-5)
    for i in 1:length(wfprobarr[1])
        BoolAllzero = prod([ abs(wfprob[i])<tol for wfprob in wfprobarr ])
        if !BoolAllzero
            val_wfprobarr = [ abs(wfprob[i])<tol ? format("{:10d}",0) :  format("{:10.3e} +i {:10.3e}", real(wfprob[i]), imag(wfprob[i]) ) for wfprob in wfprobarr ]
            println( "[ $i ] : $(val_wfprobarr...)" )
        end
    end
end

function ProbCompareAbsArr( wfprobarr::Vector ; tol=1e-5)
    for i in 1:length(wfprobarr[1])
        BoolAllzero = prod([ abs(wfprob[i])<tol for wfprob in wfprobarr ])
        if !BoolAllzero
            val_wfprobarr = [ abs(wfprob[i])<tol ? format("{:10d}",0) :  format("{:10.3e}", abs(wfprob[i])) for wfprob in wfprobarr ]
            println( "[ $i ] : $(val_wfprobarr...)" )
        end
    end
end

@inline function WFCompareAbsdiff( wf::WF, wf2::WF )
    return abs.(wf.Probamp) - abs.(wf.Probamp) 
end

function ShowWFBasis( wf::WF)
    for i in 1:wf.Dim
        binst   = wf.BasisBinstate[i]-1
        println( "WF-basis[ $i ] : $(binst+1) [$(getbittaili(binst))]" )
    end
end

function HashingIntToBinall( ; norb  = 4, outputlevel=0)
    dim     = 2^norb
    nsector = norb
    HashfIntToBin   = Int[]
    for binst in 0:dim-1
        ibinst  = binst+1 
        push!( HashfIntToBin, ibinst )
        if outputlevel > 1
            # println( getbittaili(binst) )
            println( format("ibinst=$(ibinst) -> HashfIntToBin[ind={}] = $(ibinst) [{}]", ibinst, getbittaili(binst)) )
        end
    end
    if outputlevel > 0 
        @show typeof(HashfIntToBin)
        @show HashfIntToBin
    end
    return HashfIntToBin
end

function HashingInvBinallToInt( hashf::Vector{Int} ; outputlevel=0)
    dim     = length(hashf)
    HashfBinallToInt   = zeros(Int,dim)
    nbasis  = dim
    for ib in 1:nbasis
        ist  = hashf[ib]
        HashfBinallToInt[ist]    = ib
        # push!(HashfBinallToInt, ib)
        if outputlevel > 1 
            println( format("HashfBinallToInt[ibinst={}] : {} [{}]", ist, ib, getbittaili(ib)) )
        end
    end
    if outputlevel > 0 
        @show typeof(HashfBinallToInt)
        @show HashfBinallToInt
    end
    return HashfBinallToInt
end


function HashingIntToSecarrBin( ; norb  = 4, outputlevel=0)
    dim     = 2^norb
    nsector = norb
    HashfIntToSecarrBin   = Vector{Int}[ [] for i in 0:nsector]
    for binst in 0:dim-1
        ntot    = AddIntIndex!( HashfIntToSecarrBin, binst )    # add ibinst=binst+1  in [ntot+1]-component
        if outputlevel > 1
            # println( getbittaili(binst) )
            println( format("ind={} pushed to HashfIntToSecarrBin[$(ntot+1)] (N=$(ntot)) : {}", binst+1, HashfIntToSecarrBin[ntot+1]) )
        end
    end
    if outputlevel > 0 
        @show typeof(HashfIntToSecarrBin)
        @show HashfIntToSecarrBin
    end
    return HashfIntToSecarrBin
end

function HashingIntToBinallOrdered( ; norb  = 4, outputlevel=0)
    hashf   =  vcat( HashingIntToSecarrBin(;norb=norb)... )
    if outputlevel > 0 
        @show typeof(hashf)
        @show hashf
    end
    return hashf
end

function HashingInvSecarrBinToAllNparInt( hashf::Vector{Vector{Int}} ; outputlevel=0)
    dim = length(collect(Iterators.flatten(hashf)))
    HashfBinToNinfoIndbin   = Vector{Int}[ [] for i in 1:dim]
    Npar    = length(hashf)
    for i in 1:Npar
        nbasis  = length(hashf[i])
        if outputlevel > 1
            println( format("N={} :",i) )
        end
        for ib in 1:nbasis
            ist  = hashf[i][ib]
            push!(HashfBinToNinfoIndbin[ist], i, ib)
            if outputlevel > 1
                println( format("ibinst=$(ist) [$(getbittaili(ist-1))] -> HashfBinToNinfoIndbin[ibinst={}] = [{} {}]", ist, i, ib) )
            end
        end
    end
    if outputlevel > 0 
        @show typeof(HashfBinToNinfoIndbin)
        @show HashfBinToNinfoIndbin
    end
    return HashfBinToNinfoIndbin
end

function HashingInvSecarrBinToInt( hashfsec::Vector{Vector{Int}} ; outputlevel=0)
    dim = length(collect(Iterators.flatten(hashfsec)))
    HashfInvSecarrBinToInt  = zeros(Int,dim)
    Nsec    = length(hashfsec)
    for i in 1:Nsec
        nbasis  = length(hashfsec[i])
        if outputlevel > 1
            println( format("N={} :",i-1) )
        end
        for ib in 1:nbasis
            ist  = hashfsec[i][ib]
            HashfInvSecarrBinToInt[ist]  = ib
            if outputlevel > 0
                println( format("Set HashfInvSecarrBinToInt[ibinst=$(ist)] = $(ib) ") )
            end
        end
    end
    if outputlevel > 0 
        @show typeof(HashfInvSecarrBinToInt)
        @show HashfInvSecarrBinToInt
    end
    return HashfInvSecarrBinToInt
end

@inline function AddIntIndex!( arr::Vector{Vector{Int}}, st)
    ntot    = OpNtot(st)
    push!( arr[ntot+1], st+1 )
    return ntot
end

function ShowHashfsec( hashfsec::Vector{Vector{Int}} )
    Npar    = length(hashfsec)
    println("hashfsec[Npar=$(Npar)][dim_ipar]")
    println("<ipar,ibinst> -> <ibinstate> = hashfsec[ipar][istate] [occupation/binary representation of binstate (=ibinstate-1)]")
    for i in 1:Npar
        nbasis  = length(hashfsec[i])
        println( format("N={} (ipar={}) :",i-1,i) )
        for ib in 1:nbasis
            ist  = hashfsec[i][ib]
            println( format("{} -> {} [binrep={}]", (i,ib), ist, getbittaili(ist-1)) )
        end
    end
end

function ShowHashfsecDim( hashfsec::Vector{Vector{Int}} )
    dimsec_arr  = [ length(hashfdum) for hashfdum in hashfsec ]
    for (isector,dimsec_i) in Iterators.enumerate( dimsec_arr ) 
        @show (isector,dimsec_i)
    end
end

function ShowHashfInv( hashfinv::Vector{Int} )
    dim    = length(hashfinv)
    println("hashfinv[dim_all]")
    println("<ibinstate> [occupation/binary representation of binstate(=ibinstate-1)] ->  <istate_on_ipar> of <ipar/Ntot> ")
    println("HashfInv : ")
    for i in 1:dim
        println( format("ibinst={} [{}] -> {} ", i, getbittaili(i-1), hashfinv[i] ) )
    end
end

