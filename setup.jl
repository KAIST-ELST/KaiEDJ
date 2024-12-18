using Pkg

PWD = ENV["PWD"]
PROJECT_PATH_KAIEDJ = PWD*"/envs/KEDJ"
ENV["PROJECT_PATH_KAIEDJ"] = PROJECT_PATH_KAIEDJ

Pkg.activate(PROJECT_PATH_KAIEDJ)

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
                "ThreadedSparseCSR",
                "SparseMatricesCSR",
                "Plots",
                "TickTock",
                "TimerOutputs",
                "FFTW",
                "JSON",
                "Dierckx",
                "ImageFiltering",
                "HDF5",
                "DFTforge"
                ]
ModuleNames = [ ]
nmod = length(ModuleNames)
for (imod, modname) in Iterators.enumerate(ModuleNames )
    println( "[$(imod)/$(nmod)] $(modname) ..." )
    Pkg.add(modname)
end

include("scripts/patch_hotfix.jl")
