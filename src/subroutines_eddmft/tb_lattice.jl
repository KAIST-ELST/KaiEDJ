using LinearAlgebra

export read_wann90
export hr_to_hk
export hr_dftforge_to_hk
export get_hlocal
export ReadKPOINTS
export LinRangeVector
export LinRangeVector!
export SetKPOINTSDisp
export GetKarrAllFromKpSymLines

export calc_cart_distance_kpsymlines
export get_cum_index_cart

export TBWann

struct TBWann
    fname_hr
    nbasis_hr
    hr
    hk
    nk1grid::Int
    k1arr
    wk1arr
    nkgrid3D::Int
    karr3D
    wkarr3D
end

function TBWann( ; fname_hr="./wannier90_hr.dat", nk1grid=8, nbasis_hr=2, 
                            ISPIN=1,
                            local_phase=false,
                            outputlevel=0 )
    println("")
    println("Reading $(fname_hr)")
    hr = read_wann90(fname_hr)
    if outputlevel > 0
        @show size(hr)
        @show typeof(hr)
    end
    xarr, wkarr = gausslegendre( nk1grid )
    k1arr   = [ (x + 1) * pi for x in xarr ]
    w1arr   = wkarr / 2.

    karr3D  = [[ [kx,ky,kz] for (kx,ky,kz) in Iterators.product(k1arr,k1arr,k1arr) ]...]
    wkarr3D = [[ wkx*wky*wkz for (wkx,wky,wkz) in Iterators.product(w1arr,w1arr,w1arr)]...]
    nkgrid3D    = length(karr3D)
    @show nkgrid3D
    @show nbasis_hr
    @time hk    = hr_to_hk( hr, nbasis_hr, karr3D )
    hk      = hr_to_hk( hr, nbasis_hr, karr3D ; ISPIN=ISPIN, local_phase=local_phase )
    return TBWann( fname_hr, nbasis_hr, hr, hk, 
                    nk1grid, k1arr, w1arr, 
                    nkgrid3D, karr3D, wkarr3D)
end

function ParseW90lines( dat_line ) 
    item_arr    = split(dat_line)
    return (    parse( Int, item_arr[1] )  ,
                parse( Int, item_arr[2] )  ,
                parse( Int, item_arr[3] )  ,
                parse( Int, item_arr[4] )  ,
                parse( Int, item_arr[5] )  ,
                parse( Float64, item_arr[6] )  ,
                parse( Float64, item_arr[7] ) 
                )
end

function read_pos( pos_file )
    readdlm( posfile ) 
end

function read_wann90( hr_dat_file )
    ## hri ~ ( na1, na2, na3, iorb, jorb, tre, tim )
    fopen   = open(hr_dat_file)
    frr = readlines(fopen)
    nlines  = length(frr) 
    hr  = Tuple{Int,Int,Int,Int,Int,Float64,Float64}[]
    for i in 1:nlines
        push!( hr, ParseW90lines( frr[i] ) )
    end
    close(fopen)
    return hr
end

function hr_to_hk( hr, norb, karr ; ISPIN=1, local_phase=false )
    ## hr ~ t * exp( - im * k * rij ) 
    nk  = length(karr)

    # if local_phase 
    #   pos_orb = readdlm( fnamePOSCORRORBITAL )
    #   norb_pos    = size(pos_orb)[1]
    #   if (norb!=norb_pos)  
    #       error( "Input-norb and POSCORRORBITAL-norb are different. $(norb) != $(norb_pos) " )
    #   end
    # end
    if ISPIN > 1 
        hk  = Matrix{ComplexF64}[ zeros(ComplexF64,norb,norb) for ik in 1:nk ]
        for (indk,ki) in Iterators.enumerate(karr)
            for (indhr, hri)  in Iterators.enumerate(hr)
                hki = hk[indk]
                navec = [ hri[1:3]... ]
                i, j        = hri[4:5]
                tr, ti      = hri[6:7]
                # rij         = pos_orb[j,:] - pos_orb[i,:]
                hki[i,j]    += (tr + im * ti ) * exp( - im * sum(ki .* (navec)) ) 
                # @show ( hri, "ki=",ki, "eth=", exp( - im * sum(ki .* (navec)) ) , "hkij=",hki[i,j] )
            end
        end
    else
        println("Spin-partner is being created with doubling orbital components.")
        hk  = Matrix{ComplexF64}[ zeros(ComplexF64,2*norb,2*norb) for ik in 1:nk ]
        for (indk,ki) in Iterators.enumerate(karr)
            for (indhr, hri)  in Iterators.enumerate(hr)
                hki = hk[indk]
                navec = [ hri[1:3]... ]
                i, j        = hri[4:5]
                tr, ti      = hri[6:7]
                # rij         = pos_orb[j,:] - pos_orb[i,:]
                hki[2*i-1,2*j-1]    += (tr + im * ti ) * exp( - im * sum(ki .* (navec)) ) 
                hki[2*i  ,2*j  ]    += (tr + im * ti ) * exp( - im * sum(ki .* (navec)) ) 
                # @show ( hri, "ki=",ki, "eth=", exp( - im * sum(ki .* (navec)) ) , "hkij=",hki[i,j] )
            end
        end
    end
    return hk
end


function hr_dftforge_to_hk( hr_info_dftforge, norb, karr ; ISPIN=1, local_phase=false )
    ## hr ~ t * exp( - im * k * rij ) 
    nk  = length(karr)

    hr_dftforge = hr_info_dftforge.scf_r.Hks_R[1]          # Array of norb*norb matrix
    rvec_as_mat = hr_info_dftforge.scf_r.R_vector_mat[1]   # Array of 3-comp vector

    if ISPIN > 1 
        hk  = Matrix{ComplexF64}[ zeros(ComplexF64,norb,norb) for ik in 1:nk ]
        for (indk,ki) in Iterators.enumerate(karr)
            hk[indk]    += sum( exp.( -im * rvec_as_mat * ki ) .* hr_dftforge )
            hk[indk]    += sum( exp.( -im * rvec_as_mat * ki ) .* hr_dftforge )
        end
    else
        println("Spin-partner is being created with doubling orbital components.")
        hk  = Matrix{ComplexF64}[ zeros(ComplexF64,2*norb,2*norb) for ik in 1:nk ]
        ind_odd  = 1:2:2*norb
        ind_even = 2:2:2*norb
        for (indk,ki) in Iterators.enumerate(karr)
            hk[indk][ind_odd,ind_odd]   += sum( exp.( -im * rvec_as_mat * ki ) .* hr_dftforge )
            hk[indk][ind_even,ind_even] += sum( exp.( -im * rvec_as_mat * ki ) .* hr_dftforge )
        end
    end
    return hk
end


@inline function HkSegmentFromHr( (a1, a2, a3, i, j, tr, ti), ki, norb )
    hki = zeros(ComplexF64,norb,norb)
    hki[i,j]    = (tr + im * ti ) * exp( - im * dot(ki ,[a1, a2, a3]) ) 
    return hki
end

function hr_to_hk_thread( hr, norb, karr )
    ## hr ~ t * exp( - im * k * rij ) 
    nk  = length(karr)
    hk  = Matrix{ComplexF64}[ zeros(ComplexF64,norb,norb) for ik in 1:nk ]

    nth = Threads.nthreads()
    chunks = collect(Iterators.partition( 1:nk, div(nk,nth) ))

    Threads.@threads for ith in 1:nth
        chunk   = chunks[ith]
        hkchunk = hk[chunk]
        kchunk  = karr[chunk]
        for (hki,ki) in Iterators.zip(hkchunk,kchunk)
            hki = sum(  map( x -> HkSegmentFromHr( x, ki, norb ), hr ), dims=1 )[1]
        end

    end
    return hk
end

function get_hlocal( hk, wkarr ) 
    return sum( [ hki * wkarri for (hki,wkarri) in Iterators.zip( hk, wkarr ) ] )
end

@inline function ParseKPOINTSlines( fline ) 
    kpos        = split(split(fline, "!")[1])
    kposval     = [ parse(Float64,ki) for ki in kpos ]
    kname       = split(fline, "!")[2]
    kname       = filter(x -> !isspace(x), kname)
    @show ( kposval, kname )
    return kposval, kname
end

function LinRangeVector( v1, v2, npoint ) 
    k1, k2, k3 = [ collect(a) for a in LinRange.( v1, v2, npoint) ]
    karr = Vector{Float64}[]

    for ik in 1:npoint
        push!( karr, [ k1[ik],k2[ik],k3[ik] ] )
    end
    return karr 
end

function LinRangeVector!( karr, v1, v2, npoint ) 
    k1, k2, k3 = [ collect(a) for a in LinRange.( v1, v2, npoint) ]
    for ik in 1:npoint
        push!( karr, [ k1[ik],k2[ik],k3[ik] ] )
    end
end

mutable struct KpSymLine
    nkpsym::Int64
    kpsym1pos::Vector{Float64}
    kpsym1lab::String
    kpsym2pos::Vector{Float64}
    kpsym2lab::String
    karrsymline::Vector{Vector{Float64}}
    mode::String
    # KpSymLine( (symkpos, symklab) ) =  new( 20, symkpos[1], symklab[1], symkpos[2], symklab[2], LinRangeVector(symkpos[1], symkpos[2]), "f" )
    # KpSymLine( symkpos1, symklab1, symkpos2, symklab2 ) =  new( 20, symkpos1, symklab1, symkpos2, symklab2, LinRangeVector(symkpos1, symkpos2), "f" )
    # KpSymLine( (symkpos1, symklab1), (symkpos2, symklab2) ) =  new( 20, symkpos1, symklab1, symkpos2, symklab2, LinRangeVector(symkpos1, symkpos2), "f" )
    KpSymLine( (symkposlab1, symkposlab2) )         =  new( 20,  symkposlab1[1], symkposlab1[2], symkposlab2[1], symkposlab2[2], LinRangeVector(symkposlab1[1], symkposlab2[1],20), "f" )
    KpSymLine( nkp, (symkposlab1, symkposlab2) )    =  new( nkp, symkposlab1[1], symkposlab1[2], symkposlab2[1], symkposlab2[2], LinRangeVector(symkposlab1[1], symkposlab2[1],nkp), "f" )
end

function ReadKPOINTS( fname="./KPOINTS", nk=0 )
    fopen   = open(fname)
    flines = readlines(fopen)
    nlines  = length(flines) 
    #karr    = Tuple{Float64,Float64,Float64}[]
    symkparr    = Tuple{Vector{Float64},String}[]

    if nk == 0 
        println( flines[2] )
        nk  = parse( Int, split( flines[2] )[1] )
    end

    for i in 5:nlines
        println( "fline $i len=$(length(flines[i])): ", flines[i] )
        if length(flines[i]) > 3
            push!( symkparr, ParseKPOINTSlines( flines[i] ) )
        end
    end

    println("KPOINTS : " ) 
    for kp in symkparr
        println( kp ) 
    end
    println("")
    
    close(fopen)
    return symkparr, nk
end

function SetKPOINTSDisp( ; fname="./KPOINTS", nkpsymline=0, favec="./AVEC", avec=[], fbvec="./BVEC", bvec=[] )  
    kpline_raw, nk_raw  = ReadKPOINTS(fname)
    if nkpsymline == 0
        nkpsymline = nk_raw
    end
    nkpline_raw = length(kpline_raw)
    nsymline    = div( nkpline_raw, 2 )
    @show nsymline

    KParr  = KpSymLine[]
    for ikp in 1:nsymline
        push!( KParr, KpSymLine(nkpsymline,(kpline_raw[2*ikp-1], kpline_raw[2*ikp])) )
    end

    if avec == []
        try 
            avec    = readdlm("AVEC")
        catch
        end
    end
    if bvec == []
        try 
            bvec    = reciprocal_lattice_vectors( avec ) 
        catch
            bvec    = readdlm("BVEC")
        end
    end

    println( "lattice vectors (a_i = A[i,:]) ")
    @show avec
    println( "reciprocal lattice vectors (b_i = B[i,:]) ")
    @show bvec

    # for kp_i in KParr
    #     @show kp_i
    # end
    return KParr
end

function GetKarrAllFromKpSymLines( KPSymLines::Vector{KpSymLine} )
    karr    = Vector{Float64}[]
    kscale  = KPSymLines[1].mode=="f" ? pi*2.0 : 1.0
    karr    = vcat( [ kp_i.karrsymline for kp_i in KPSymLines  ]... ) * kscale 
    @show karr
    return karr
end

function reciprocal_lattice_vectors(lattice_vectors::Matrix)
    # Calculate the volume of the unit cell
    volume = abs(det(lattice_vectors))
    @show volume

    # Calculate the reciprocal lattice vectors
    b1 = 2*pi * cross(lattice_vectors[2, :], lattice_vectors[3, :]) / volume
    b2 = 2*pi * cross(lattice_vectors[3, :], lattice_vectors[1, :]) / volume
    b3 = 2*pi * cross(lattice_vectors[1, :], lattice_vectors[2, :]) / volume

    # println( "reciprocal lattice vectors (b_i = B[i,:]) ")
    # @show transpose([ b1 b2 b3 ])

    return transpose([b1 b2 b3])
end
# # Example usage:
# # Define the lattice vectors as columns of a 3x3 matrix
# # Each column represents a lattice vector
# lattice_vectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
# 
# reciprocal_vectors = reciprocal_lattice_vectors(lattice_vectors)
# println("Reciprocal lattice vectors:")
# println(reciprocal_vectors)

function read_lattice_vectors_from_poscar( ; file_path::String="./")
    open(joinpath(file_path,"POSCAR"), "r") do file
        lines = readlines(file)
        
        # Read the scaling factor
        scaling_factor = parse(Float64, strip(lines[2]))
        
        # Read the lattice vectors and apply the scaling factor
        a1 = parse.(Float64, split(strip(lines[3]))) * scaling_factor
        a2 = parse.(Float64, split(strip(lines[4]))) * scaling_factor
        a3 = parse.(Float64, split(strip(lines[5]))) * scaling_factor
        
        return [a1, a2, a3]
    end
end

function calc_cart_distance_kpsymlines( KPSymLines::Vector{KpSymLine} )
    bvec = get_bvec()
    dist_arr    = []
    for kpsymline in KPSymLines
        k1pos_cart  = bvec' * kpsymline.kpsym1pos
        k2pos_cart  = bvec' * kpsymline.kpsym2pos
        dist_cart   = norm(k2pos_cart .- k1pos_cart)
        println( "Distance $(kpsymline.kpsym1lab)-$(kpsymline.kpsym2lab) : ", dist_cart)
        push!( dist_arr, norm(k2pos_cart .- k1pos_cart) ) 
    end
    return dist_arr
end

function get_bvec(; avec=[])
    if avec == [] 
        try 
            avec    = readdlm("AVEC")
            println( "Reading AVEC ...")
        catch
            avec    = read_lattice_vectors_from_poscar()
            println( "Reading POSCAR ...")
        end
    end
    @show avec
    bvec    = reciprocal_lattice_vectors(avec)
    return bvec
end

function get_cum_index_cart( KPSymLines::Vector{KpSymLine} )
    dist_arr    = calc_cart_distance_kpsymlines(KPSymLines)
    nkpsymline  = length(KPSymLines)
    nkpsym      = KPSymLines[1].nkpsym
    println( "nkp for each symline = ", nkpsym )
    cum_dist_arr    = cumsum(vcat(0,dist_arr))
    cum_ind_label   = cum_dist_arr  / cum_dist_arr[end]
    symlabels   = [KPSymLines[1].kpsym1lab, KPSymLines[1].kpsym2lab]
    for kpsymline in KPSymLines[2:end]
        if symlabels[end] != kpsymline.kpsym1lab
            symlabels[end] = symlabels[end] *"|"* kpsymline.kpsym1lab
        end
        push!( symlabels, kpsymline.kpsym2lab )
    end
    # @show [ collect(LinRange(0,dist_arr[i],3)) for i in 1:nkpsymline]
    # @show vcat( [ collect(LinRange(0,dist_arr[i],nkpsym)) for i in 1:nkpsymline]... )
    # @show vcat( [ collect(LinRange(0,dist_arr[i],nkpsym)) for i in 1:nkpsymline]... )
    cum_ind_rep = vcat( [ collect(LinRange(0,dist_arr[i],2)) .+ cum_dist_arr[i] for i in 1:nkpsymline]... )
    @show cum_ind_rep
    @show cum_ind_rep / cum_ind_rep[end]
    cum_index   = vcat( [ collect(LinRange(0,dist_arr[i],nkpsym)) .+ cum_dist_arr[i] for i in 1:nkpsymline]... )
    return cum_index / cum_index[end], cum_ind_label, symlabels
end
