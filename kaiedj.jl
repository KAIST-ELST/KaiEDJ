import KaiEDJ

include("src/load_dftforge.jl")

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

