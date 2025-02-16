
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

import DFTforge

module KaiEDMFT

import DFTforge

# import Plots
# import FFTW
# import JSON
# import Dierckx
# using JSON
# using ImageFiltering
# using Dierckx
# using Distributed
# using DelimitedFiles
# using Printf
# using LinearAlgebra
# using TOML

include("subroutines_mft/subroutines_edmftf.jl")        #### EDMFTF related functions
include("subroutines_mft/subroutines_etc.jl")           #### etc. function
include("subroutines_mft/subroutines_io_gen.jl")        #### processing input/ouput
include("subroutines_mft/subroutines_dmft.jl")          #### dmft
include("subroutines_mft/subroutines_grid.jl")          
include("subroutines_mft/subroutines_fourier.jl")          
include("subroutines_mft/subroutines_jx.jl")          


@everywhere using LinearAlgebra
using Distributed
using DFTforge.DFTrefinery
using DFTforge.DFTcommon;

# Julia 1.0
using Statistics
#

check_print_version()

@everywhere using Distributed
@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTforge.DFTcommon


## 0.1 Read input from argument & TOML file

arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
arg_input = parse_TOML(arg_input.TOMLinput,arg_input)

# let argument override
arg_input = parse_input(ARGS,arg_input)


end

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