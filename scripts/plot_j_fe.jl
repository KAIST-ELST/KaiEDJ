using Plots
using DelimitedFiles

Jdat_1_1 = readdlm("Jx_1_1.dat")

cutoffnumber = 180
x1 = Jdat_1_1[1:cutoffnumber,1]
J11 = Jdat_1_1[1:cutoffnumber,2]

x1 = convert( Vector{Float64}, x1 )
J11= convert( Vector{Float64}, J11 )

plot(  x1, J11, lc=:red,  mc=:red,  markershape = :circle, label="atom1-atom1", xlabel="Distance [Angst]", ylabel="J [meV]")

savefig( "output_dmftmft_fe.png" )
savefig( "output_dmftmft_fe.pdf" )
