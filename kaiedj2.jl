import KaiEDJ
#using KaiEDJ.KaiEDJ_MFT

#KaiEDJ_MFT.check_print_version()

using Statistics

using Distributed
using DFTforge.DFTrefinery
using DFTforge.DFTcommon;




@everywhere using Distributed
@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTforge.DFTcommon

#include("src/load_dftforge.jl")


arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
arg_input = parse_TOML(arg_input.TOMLinput,arg_input)

# let argument override
arg_input = parse_input(ARGS,arg_input)



Calculation_mode    = arg_input.Optional["KaiEDJ"].Calculation_mode

if (Calculation_mode == "DMFT")
    DMFT_solver         = arg_input.Optional["KaiEDJ"].DMFT_solver
    if (DMFT_solver == "ED")
        include("src/eddmft_mag.jl")
    elseif (DMFT_solver == "EDMFTF_ctqmc" || DMFT_solver=="ComCTQMC")
        include("src/dmft_mft.jl")
    end
elseif (Calculation_mode == "DMFT+MFT" || Calculation_mode == "Jx-DMFT"|| Calculation_mode == "Jx-DMFT-private"|| Calculation_mode == "Spinwave"|| Calculation_mode == "Magnon")
    include("src/dmft_mft.jl")
else
    error("Invalid 'Calculation_mode' ($(Calculation_mode)).")
end

