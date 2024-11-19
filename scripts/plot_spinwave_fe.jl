using Plots
using DelimitedFiles

rawdat = readdlm("magnon.dat")[2:end,:]

kdat = rawdat[:,1]
edat = rawdat[:,2]

kdat = convert( Vector{Float64}, kdat )
edat = convert( Vector{Float64}, edat )

xtick = [0.0000, 1.5409, 2.6427, 4.5315, 6.7268, 8.2759]
xlabel= [ "G", "N", "P", "G", "H", "N" ] 

plot(  kdat, edat, label="w(k)", xticks=(xtick,xlabel))

savefig( "output_spinwave_fe.png" )
savefig( "output_spinwave_fe.pdf" )
