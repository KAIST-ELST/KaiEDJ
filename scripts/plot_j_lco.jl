using Plots
using DelimitedFiles

Jdat_1_1 = readdlm("Jx_1_1.dat")
Jdat_1_2 = readdlm("Jx_1_2.dat")

cutoffnumber = 10
x1 = Jdat_1_1[1:cutoffnumber,1]
J11 = Jdat_1_1[1:cutoffnumber,2]

x2 = Jdat_1_2[1:cutoffnumber,1]
J12 = Jdat_1_2[1:cutoffnumber,2]

x1 = convert( Vector{Float64}, x1 )
x2 = convert( Vector{Float64}, x2 )
J11= convert( Vector{Float64}, J11 )
J12= convert( Vector{Float64}, J12 )

plot(  x1, J11, lc=:red,  mc=:red,  markershape = :circle, label="atom1-atom1", xlabel="Distance [Angst]", ylabel="J [meV]")
plot!( x2, J12, lc=:blue, mc=:blue, markershape = :square, label="atom1-atom2")

savefig( "output_dmftmft_lco.png" )
savefig( "output_dmftmft_lco.pdf" )
