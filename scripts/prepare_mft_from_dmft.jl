ind_QMC = isdir("imp_1")

if ind_QMC
    include("copy_dmft_qmc.jl")
else
    include("copy_dmft_ed.jl")
end
