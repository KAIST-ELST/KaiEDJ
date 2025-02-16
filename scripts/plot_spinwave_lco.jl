using KaiEDJ: Plots
using DelimitedFiles

rawdat = readdlm("magnon.dat")[2:end,:]

kdat = rawdat[:,1]
edat = rawdat[:,3]

kdat = convert( Vector{Float64}, kdat )
edat = convert( Vector{Float64}, edat )

xtick = [0.0000, 0.5811, 1.4035, 1.9880, 2.5691, 3.3949]
xlabel= ["P", "M", "X", "P", "G", "X"]

Plots.plot(  kdat, edat, label="w(k)", xticks=(xtick,xlabel), xlabel="k", ylabel="Energy", title="Spin Wave Dispersion")

Plots.savefig( "output_spinwave_lco.png" )
Plots.savefig( "output_spinwave_lco.pdf" )
