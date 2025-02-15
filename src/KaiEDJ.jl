
println("[info] Loading KaiEDJ")




module KaiEDJ

import Plots
import JSON
import Dierckx
import FFTW
using JSON
using ImageFiltering
using Dierckx
using Distributed
#using DelimitedFiles
using Printf
using LinearAlgebra
using TOML
#using Formatting
using Format
using LinearAlgebra
using SparseArrays
using DelimitedFiles
using FastGaussQuadrature
using BenchmarkTools
using Optimization
using OptimizationOptimJL
using Optim
using TickTock
using HDF5
using TOML
using JLD2


# module MyBase
include("subroutines_eddmft/basic.jl")
include("subroutines_eddmft/printout.jl")
include("subroutines_eddmft/binary_basis.jl")
include("subroutines_eddmft/wavefunction.jl")
include("subroutines_eddmft/operators.jl")
include("subroutines_eddmft/operators_phys.jl")
include("subroutines_eddmft/hamiltonian.jl")
include("subroutines_eddmft/green.jl")
include("subroutines_eddmft/bethe.jl")
include("subroutines_eddmft/discretization.jl")
include("subroutines_eddmft/lib_continued_fraction.jl")
include("subroutines_eddmft/tb_lattice.jl")
include("subroutines_eddmft/ed_solver.jl")
include("subroutines_eddmft/observables.jl")
# end

println("[info] Loading KaiEDJ finished")

end