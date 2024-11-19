## import DFTforge
## 
## import Plots
## import FFTW
## import JSON
## import Dierckx
## using JSON
## using ImageFiltering
## using Dierckx
## using Distributed
## using DelimitedFiles
## using Printf
## using LinearAlgebra
## using TOML
## 
## include("subroutines_mft/subroutines_edmftf.jl")        #### EDMFTF related functions
## include("subroutines_mft/subroutines_etc.jl")           #### etc. function
## include("subroutines_mft/subroutines_io_gen.jl")        #### processing input/ouput
## include("subroutines_mft/subroutines_dmft.jl")          #### dmft
## include("subroutines_mft/subroutines_grid.jl")          
## include("subroutines_mft/subroutines_fourier.jl")          
## include("subroutines_mft/subroutines_jx.jl")          
## 
## @everywhere using LinearAlgebra
## using Distributed
## using DFTforge.DFTrefinery
## using DFTforge.DFTcommon;
## 
## # Julia 1.0
## using Statistics
## #
## 
## check_print_version()
## 
## @everywhere using Distributed
## @everywhere import DFTforge
## @everywhere using DFTforge.DFTrefinery
## @everywhere using DFTforge.DFTcommon
## 
## 
## arg_input = DFTcommon.Arg_Inputs();
## arg_input = parse_input(ARGS,arg_input)
## 
## arg_input = parse_TOML(arg_input.TOMLinput,arg_input)
## 
## # let argument override
## arg_input = parse_input(ARGS,arg_input)


basisTransform_rule = basisTransform_rule_type()
hamiltonian_info = [];
if (DFTcommon.OpenMX == arg_input.DFT_type || DFTcommon.EcalJ == arg_input.DFT_type)
  hamiltonian_info  = set_current_dftdataset(arg_input.result_file, arg_input.result_file_dict, arg_input.DFT_type, arg_input.spin_type, basisTransform_rule)
elseif (DFTcommon.Wannier90 == arg_input.DFT_type)
  atomnum           = arg_input.Wannier_Optional_Info.atomnum
  atompos           = arg_input.Wannier_Optional_Info.atompos
  atoms_orbitals_list   = arg_input.Wannier_Optional_Info.atoms_orbitals_list

  #hamiltonian_info = DFTforge.read_dftresult(result_file,arg_input.DFT_type,Wannier90_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
  hamiltonian_info  = set_current_dftdataset(arg_input.result_file, arg_input.result_file_dict,
                        arg_input.DFT_type, arg_input.Wannier90_type, arg_input.spin_type, atoms_orbitals_list, atomnum, atompos, basisTransform_rule)
  #hamiltonian_info = set_current_dftdataset(scf_r, arg_input.DFT_type, spin_type,basisTransform_rule)
end

DFTforge.pwork(set_current_dftdataset,(hamiltonian_info, 1));

dump( hamiltonian_info )


# ENV["PROJECT_PATH_ED"]="/Users/jun/dmft/ed/github/KaiED/envs/KED"
include("/Users/jun/github/KaiEDJ/src/mybase.jl")
@show tick()

toml_filename   = arg_input.TOMLinput
TOML_param      = TOML.parsefile(toml_filename)["KaiEDJ"]
@show TOML.parsefile(toml_filename)
@show TOML_param
# @show arg_input.Optional["KaiEDJ"]

dirnameData   = "./data/"
run( `mkdir -p $(dirnameData)` )

# hr  = read_wann90( "/Users/jun/dmft/ed/github/KaiED/LCO_2u/wannier90_hr_cut.dat" )
# 
# @show size(hr)
# @show typeof(hr)

nkarr      = TOML_param["KgridNum"]
@show nkarr
@show typeof(nkarr)


function GenGaussLegedreQuadNormalized( nk )
    xarr, warr  = gausslegendre( nk )        ## Domain is [-1,1] initially.
    xarr_norm   = [ (x + 1) for x in xarr ] ./ 2.0   ## Domain is now [0,2]
    warr_norm   = warr / 2.
    return xarr_norm, warr_norm
end

quadarr3D   = [ GenGaussLegedreQuadNormalized( nk_i ) for nk_i in nkarr ]
karr3D      = [[ [kx,ky,kz]*2*pi for (kx,ky,kz) in Iterators.product(quadarr3D[1][1],quadarr3D[2][1],quadarr3D[3][1]) ]...]
wkarr3D     = [[ wkx*wky*wkz for (wkx,wky,wkz) in Iterators.product(quadarr3D[1][2],quadarr3D[2][2],quadarr3D[3][2])]...]
nk3D     = length(karr3D)
@show nk3D

######## Correlated Orbital Info #########
Corr_atom_Ind       = TOML_param["Corr_atom_Ind"]
Corr_orbital_Ind    = TOML_param["Corr_orbital_Ind"]
Nimp    = 2 
Nimp    = length(TOML_param["Corr_atom_Ind"])
norb    = 1
norb    = length(TOML_param["Corr_orbital_Ind"][1])
nspin   = 2
nspinorb= norb * nspin
ISPIN   = 1 
NSpinOrbNimp  = [ nspin*length(orb_indarr) for orb_indarr in Corr_orbital_Ind ]
PartitionNSpinOrbNimp = cumsum([ 0, NSpinOrbNimp... ])
# @show NSpinOrbNimp
# @show PartitionNSpinOrbNimp
IndCorrOrbNimp  = [ [1,2] , [3,4] ] 
IndCorrOrbNimp  = [ collect(PartitionNSpinOrbNimp[iimp-1]+1:PartitionNSpinOrbNimp[iimp])  for iimp in 2:length(PartitionNSpinOrbNimp) ]
@show Nimp
@show norb
@show nspin
@show IndCorrOrbNimp

norblatt        = Nimp*norb
nspinorblatt    = Nimp*nspinorb

######## Get TB model from wann #########
# @time hk  = hr_to_hk( hr, norblatt, karr3D ; fnamePOSCORRORBITAL="/Users/jun/dmft/ed/github/KaiED/LCO_2u/POSCORRORBITAL" )
# @time hk  = hr_to_hk( hr, norblatt, karr3D )
@time hk  = hr_dftforge_to_hk( hamiltonian_info, norblatt, karr3D )



######## TB Model Check in Real-frequency #########
NReFreq   = 800
ReFreqGridVal   = LinRange( -10, 10, NReFreq ) 
epsilon         = 0.1
ReFreqGrid      = ReFreqGridVal .+ im * epsilon
WriteVec( ReFreqGridVal, dirnameData*"ReFreqGridVal" )

chem    = hamiltonian_info.scf_r.ChemP
@show chem

@time @show ECorrOrb    = real(get_hlocal( hk, wkarr3D ))
@show typeof(ECorrOrb)
WriteMat( ECorrOrb, "ECorrOrb")
@time G0local_w  = GetGreenLocalFromHkGrid( hk, wkarr3D, ReFreqGrid, collect(I(nspinorblatt)*chem) )



######## Imaginary Frequency Green Function Construction #########
beta    = TOML_param["beta"]       # 128
NImFreq = TOML_param["NImFreq"]    # 4*beta
@show (beta, NImFreq)
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im
WriteVec( ImFreqGridVal, dirnameData*"ImFreqGridVal" )

######## Bath Orbital Info #########
nbath   = 6
try 
    global nbath   = TOML_param["nbath"]       # 6
    println("nbath (from TOML) = $(nbath)")
catch
    global nbath   = 6
    println("nbath (default) = $(nbath)")
end
ntot    = nspinorb + nbath
nbathHalf   = div(nbath,2)

######## Bath Orbital Parameters Initialization #########
BPNimp  = [ BathParams(nspinorb,nbath) for i in 1:Nimp ]
println("")
try 
    ReadBathParam!( BPNimp )
catch
    InitBathParam!( BPNimp )
end
ShowBathParam( BPNimp )
WriteBathParam( BPNimp ; fname=dirnameData*"BathParams_init.dat")


######## Setting Up Hamiltonian Operators #########
U           = TOML_param["imp1_U"]  # 8. # 2  : Hubbard interaction strength
JHund       = TOML_param["imp1_J"]  # 0. # 0.167     : Hund coupling strength 
chem_int    = 0.5*U       # 0.5*U for single-band half-filling , 2.5*U-5.0*JHund for t2g-multiband half-filling
opccaavec   = [ GetOpUSlaterKanamori( ; U=U, JHund=JHund, norb=norb ) ]


@show (norb, nspinorb, nbath, nkarr, NImFreq, NReFreq )
@show (chem, chem_int, U, JHund)



######## Setting Up Initial Green Function #########

@time G0iw      = GetGreenLocalFromHkGrid( hk, wkarr3D, ImFreqGrid, collect(I(nspinorblatt)*chem) )
Hybiw_init      = GetHybFromGreenLocalFromHkGrid( hk, wkarr3D, ImFreqGrid, collect(I(nspinorblatt)*chem), ECorrOrb )
# Hybimpiw_init   = GetDeltaHybDiscGrid( ebathl, Vil, ImFreqGrid )

rhodiag_G0latt  = real(diag(GetDensityMatrixFromGreenImFreq( G0iw, ImFreqGrid )))
nlatt_init  = [ sum(rhodiag_G0latt), rhodiag_G0latt... ]
@show nlatt_init 

Hybiw   = deepcopy(Hybiw_init)

Hybw_init       = GetHybFromGreenLocalGrid( G0local_w, ReFreqGrid, ECorrOrb - I*chem )
# Hybimpw_init    = GetDeltaHybDiscGrid( ebathl, Vil, ReFreqGrid )

outputlevel = 0

for iorb in 1:nspinorb*Nimp
    fnameTailOrb    = "_$(iorb)_$(iorb).dat"
    WriteVecMat( real(ReFreqGrid),       GetijarrayFromVecMat( G0local_w,    iorb, iorb ), dirnameData*"G0w"*fnameTailOrb )
    WriteVecMat( real(ReFreqGrid),       GetijarrayFromVecMat( Hybw_init,    iorb, iorb ), dirnameData*"Hybw"*fnameTailOrb )
    # WriteVecMat( real(ReFreqGrid),       GetijarrayFromVecMat( Hybimpw_init, iorb, iorb ), dirnameData*"Hybimpw"*fnameTailOrb )

    # WriteVecMat( real(ImFreqGridVal),    GetijarrayFromVecMat( G0iw    ,     iorb, iorb ), dirnameData*"G0iw"*fnameTailOrb )
    # WriteVecMat( real(ImFreqGridVal),    GetijarrayFromVecMat( Hybiw_init,   iorb, iorb ), dirnameData*"Hybiw"*fnameTailOrb )
    # WriteVecMat( real(ImFreqGridVal),    GetijarrayFromVecMat( Hybimpiw_init,iorb, iorb ), dirnameData*"Hybimpiw"*fnameTailOrb )
end

Hybimpw_init = [GetDeltaHybDiscGrid( BP.ebathl, BP.Vil, ReFreqGrid ) for BP in BPNimp ]
for iimp in 1:Nimp
    for iorb in 1:nspinorb
        fnameTailOrbImp    = "_$(iorb)_$(iorb)_imp$(iimp).dat"
        WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( Hybimpw_init[iimp],  iorb, iorb ), dirnameData*"Hybimpw"*fnameTailOrbImp )
    end
end

######## Starting DMFT Calculations #########

ndmftloop   = 15
for idmft in 1:ndmftloop
    init_println( "DMFT loop ($(idmft))" )
    BathOpt = idmft<2 ? false : true

    Glattnewiw, Hyblattnewiw, SelfClusteriw, BPNimp_new, GSSectorInfoNimp, nIterOptNimp, costFinalNimp, obs  = 
        SolveEDTbDegenNimp( hk, karr3D, wkarr3D, Hybiw, BPNimp, ImFreqGrid, chem_int, chem, ECorrOrb, opccaavec ; 
                            TolEGS=1e-7, BathOpt=BathOpt, SelfDiag=false, SelfTrev=false, GimpDiag=true, SpinSep=true, IterTag=idmft,
                            IndCorrOrbNimp=IndCorrOrbNimp)
    global Hybiw    = Hyblattnewiw
    global BPNimp   = deepcopy( BPNimp_new )

    WriteBathAll( BPNimp ; iiter=idmft )
    WriteBathParam( BPNimp ; fname=dirnameData*"BathParams_i$(idmft).dat")


    Glattneww, Hyblattneww, SelfClusterw, BPNimp_new, GSSectorInfoNimp_dum, nIterOptNimp_dum, costFinalNimp_dum, obs_dum  = 
        SolveEDTbDegenNimp( hk, karr3D, wkarr3D, Hybw_init, BPNimp, ReFreqGrid, chem_int, chem, ECorrOrb, opccaavec ; 
                            GSSectorInfoNimp=GSSectorInfoNimp, BathOpt=false, IndCorrOrbNimp=IndCorrOrbNimp, show_obs=false, IterTag=idmft )


    for iorb in Iterators.flatten(IndCorrOrbNimp)
        fnameTailOrb    = "_$(iorb)_$(iorb)_i$(idmft).dat"
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Gimpiw, iorb, iorb ),  dirnameData*"Gimpiw"*fnameTailOrb )
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( G0impiw, iorb, iorb ), dirnameData*"G0impiw"*fnameTailOrb )
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Selfiw, iorb, iorb ), dirnameData*"Selfiw"*fnameTailOrb )
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Hybnewiw, iorb, iorb ), dirnameData*"Hybnewiw"*fnameTailOrb )
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Glattiw, iorb, iorb ), dirnameData*"Giw"*fnameTailOrb )
        # WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Hybimpiw, iorb, iorb ), dirnameData*"Hybimpiw"*fnameTailOrb )

        # WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( Gimpw,  iorb, iorb ), dirnameData*"Gimpw"*fnameTailOrb )
        # WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( G0impw, iorb, iorb ), dirnameData*"G0impw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( SelfClusterw,  iorb, iorb ), dirnameData*"Selfw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( Hyblattneww,  iorb, iorb ), dirnameData*"Hybneww"*fnameTailOrb )
        WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( Glattneww, iorb, iorb ), dirnameData*"Gw"*fnameTailOrb )
    end

    Hybimpw_new = [GetDeltaHybDiscGrid( BP.ebathl, BP.Vil, ReFreqGrid ) for BP in BPNimp ]
    for iimp in 1:Nimp
        for iorb in 1:nspinorb
            fnameTailOrbImp    = "_$(iorb)_$(iorb)_imp$(iimp)_i$(idmft).dat"
            WriteVecMat( real(ReFreqGrid),    GetijarrayFromVecMat( Hybimpw_new[iimp],  iorb, iorb ), dirnameData*"Hybimpw"*fnameTailOrbImp )
        end
    end
            
    fnameTail   = "_i$(idmft).dat"
    WriteMat( Hyblattnewiw,     dirnameData*"Hybnewiw"*fnameTail )

    status_occu_sz  = obs[1]
    imp_status      = obs[2]
    latt_info       = obs[3]


    # open( "status.txt", "a" ) do io
    #     writedlm( io, [ idmft sum(nIterOptNimp) sum(costFinalNimp) vcat(status_occu_sz...)... ] )
    # end
    WriteStatus( "status.txt", [ idmft sum(nIterOptNimp) sum(costFinalNimp) vcat(status_occu_sz...)... ] )
    for iimp = 1:Nimp 
        WriteStatus( "status_imp$(iimp).txt", [ idmft FlattenArr(imp_status[iimp])... ] )
        # open( "status_imp$(iimp).txt", "a" ) do io
        #     writedlm( io, [ idmft FlattenArr(imp_status[iimp])... ] )
        # end
    end
    WriteStatus( "status_latt.txt" , [ idmft latt_info... ] )
    # open( "status_latt.txt", "a" ) do io
    #     writedlm( io, [ idmft latt_info... ] )
    # end
    # WriteMat( DM,       dirnameData*"DM"*fnameTail )
    # WriteMat( DMEach,   dirnameData*"DMEach"*fnameTail )

    fnameTailHDF5   = "_i$(idmft).h5"
    WriteHDF5( dirnameData*"Selfw"*fnameTailHDF5, "Selfw", SelfClusterw )
    WriteHDF5( dirnameData*"Selfiw"*fnameTailHDF5, "Selfiw", SelfClusteriw )
    WriteHDF5( dirnameData*"Glattiw"*fnameTailHDF5, "Glattiw", Glattnewiw )

    @show laptimer()
    nIterOpt    = sum(nIterOptNimp)
    if idmft > 1 
        if nIterOpt < 0 
            error("Unknown error in optimization :: nIterOpt = $(nIterOpt)")
        elseif nIterOpt < 2*nspin*Nimp+1     # 2-factor from SpinSep in BathDisc
            println("DMFT iteration converges. ( Hybiw converges. nIterOpt = $(nIterOpt) )")
            break
        end
    end
end

@show tok()
