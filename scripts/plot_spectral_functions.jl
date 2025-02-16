using KaiEDJ: Plots
using DelimitedFiles

dat1_fullname     = "data/Gw_1_1_i7.dat"
dat2_fullname     = "data/Gw_2_2_i7.dat"
@show dat1_fullname 
@show dat2_fullname 

dat1 = transpose(readdlm( dat1_fullname ))
x1_freq = dat1[1,:]
g1_re   = dat1[2,:]
g1_im   = dat1[3,:]

dat2 = transpose(readdlm( dat2_fullname ))
x2_freq = dat2[1,:]
g2_re   = dat2[2,:]
g2_im   = dat2[3,:]

Plots.plot(  -g1_im, x1_freq, lc=:blue, label="DMFT(spin-up)")
Plots.plot!( -g2_im, x2_freq, lc=:red,  label="DMFT(spin-dn)")

Plots.savefig("output_dmft_lco.png")
Plots.savefig("output_dmft_lco.pdf")
println("The figure of dmft results is saved in output_dmft_lco.png and output_dmft_lco.pdf")
