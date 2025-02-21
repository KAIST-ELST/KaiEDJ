using Pkg

# PWD = ENV["PWD"]
# PROJECT_PATH_KAIEDJ = PWD*"/envs/KEDJ"
# ENV["PROJECT_PATH_KAIEDJ"] = PROJECT_PATH_KAIEDJ

# Pkg.activate(PROJECT_PATH_KAIEDJ)

ModuleNames = [ 
                "Formatting",
                "LinearAlgebra",
                "SparseArrays",
                "DelimitedFiles",
                "FastGaussQuadrature",
                "BenchmarkTools",
                "Optimization",
                "OptimizationOptimJL",
                "Optim",
                "Arpack",
                # "ThreadedSparseCSR",
                "SparseMatricesCSR",
                "Plots",
                "TickTock",
                "TimerOutputs",
                "FFTW",
                "JSON",
                "Dierckx",
                "ImageFiltering",
                "HDF5",
                "DFTforge",
                "JLD2"
                ]
nmod = length(ModuleNames)
for (imod, modname) in Iterators.enumerate(ModuleNames )
    println( "[$(imod)/$(nmod)] $(modname) ..." )
    Pkg.add(modname)
end


Pkg.add(url="https://github.com/KAIST-ELST/ThreadedSparseCSR.jl")
Pkg.add(path=".")

#include("scripts/patch_hotfix.jl") # Not required for DFTforge 1.4.2 and higher
