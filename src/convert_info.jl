using DelimitedFiles
using HDF5
using TOML
using LinearAlgebra



@inline function ReadHDF5Component( hdf5_fname, attr_name ) 
    fr  = h5open( hdf5_fname, "r" )
    dat = read( fr, attr_name )
    close(fr)
    return dat
end


function GetFreqData(TOML_param ; type="glattiw")
    if type=="glattiw"
        return GetFreqGreen(TOML_param)
    elseif type=="selfiw"
        return GetFreqSelf(TOML_param)
    end
end


function GetFreqGreen(TOML_param)
    imfreq_file     = TOML_param["KaiEDJ"]["imfreq_file"]
    glatt_file      = TOML_param["KaiEDJ"]["glatt_file"]

    ImFreqVal   = readdlm(imfreq_file, Float64)
    Glatt_iw    = ReadHDF5Component( glatt_file, "Glattiw" )

    return ImFreqVal, Glatt_iw
end


function GetFreqGreenFromIterId(i ; imfreq_file="data/ImFreqGridVal")
    # imfreq_file     = TOML_param["KaiEDJ"]["imfreq_file"]
    # glatt_file      = TOML_param["KaiEDJ"]["glatt_file"]
    glatt_file  = "data/Glattiw_i$(i).h5"

    ImFreqVal   = readdlm(imfreq_file, Float64)
    Glatt_iw    = ReadHDF5Component( glatt_file, "Glattiw" )

    return ImFreqVal, Glatt_iw
end


function GetFreqSelf(TOML_param)
    imfreq_file     = TOML_param["KaiEDJ"]["imfreq_file"]
    self_file       = TOML_param["KaiEDJ"]["self_file"]

    ImFreqVal   = readdlm(imfreq_file, Float64)
    Sig_iw      = ReadHDF5Component( self_file, "Selfiw" )

    return ImFreqVal, Sig_iw
end


function GetFreqSelfConvertedFromIterId(i ; imfreq_file="data/ImFreqGridVal")
    # imfreq_file         = TOML_param["KaiEDJ"]["imfreq_file"]
    # self_file           = TOML_param["KaiEDJ"]["self_file"]
    # self_component_imps = TOML_param["KaiEDJ"]["self_component_imps"]
    self_file           = "data/Selfiw_i$(i).h5"
    self_component_imps = [
                            [ [1,1], [2,2] ],
                            [ [3,3], [4,4] ]
                          ]
    @show self_component_imps

    ImFreqVal   = readdlm(imfreq_file, Float64)
    Sig_iw      = ReadHDF5Component( self_file, "Selfiw" )

    nimp        = size(self_component_imps)[1]
    @show nimp
    Sig_iw_conv = [ ConvertRut(Sig_iw,self_component_imps[imp_i]) for imp_i in 1:nimp ]

    return ImFreqVal, Sig_iw_conv
end

eV_to_Kelvin    = 11604.52500617

function GetImFreqMax(TOML_param)
    beta        = TOML_param["KaiEDJ"]["beta"]       # beta in (eV)^-1
    NImFreq     = TOML_param["KaiEDJ"]["NImFreq"]
    ImFreqMax   = (2*NImFreq-1) * pi / beta
    T_Kelvin = 11604.52500617  / beta
    dImFreq     = 2*pi / beta

    @show beta, NImFreq, ImFreqMax, T_Kelvin, dImFreq
    @show ImFreqMax
    @show ImFreqMax+dImFreq*0.5
    @show ImFreqMax+dImFreq

    return beta, NImFreq, ImFreqMax, T_Kelvin, dImFreq
end



function ConvertRut(Sig_iw,self_component)
    @show self_component

    ## Stack the arrays along column-wise.  Such as [ v_real_1, v_imag_1, v_real_2, v_imag_2, ... ].
    return reduce( hcat, [ hcat(real(Sig_iw[:,i,j]), imag(Sig_iw[:,i,j])) for (i,j) in self_component ] ) 
end


println("")
println("## Selfiw info ##" )

# toml_fname  = "kaied.toml"
# TOML_param  = TOML.parsefile(toml_fname)
# println("")
# println("## TOML info ##" )
# println(TOML_param)
# 
# println("")
# println("## Self-energy info ##" )
# ImFreqVal, Sig_iw = GetFreqSelf(TOML_param)
# # @show Sig_iw
# @show size(Sig_iw)
# @show typeof(Sig_iw)
# @show Sig_iw[1,:,:]

IterIdFinal = Int(readdlm("status.txt", Float64)[end,1])
ImFreqVal, Sig_iw_conv = GetFreqSelfConvertedFromIterId(IterIdFinal)

@show size(Sig_iw_conv)
@show typeof(Sig_iw_conv)

# println("")
# println("## Frequency info ##" )
# beta, NImFreq, ImFreqMax, T_Kelvin, dImFreq    = GetImFreqMax(TOML_param)

nimp    = size(Sig_iw_conv)[1]
for imp_i in 1:nimp
    fnameout = "SelfE_$(imp_i).dat"
    open( fnameout, "w" ) do fr
        write( fr, "#\n" )
        writedlm( fr, hcat( ImFreqVal, Sig_iw_conv[imp_i]), " " )
    end
    println( "Saved at $(fnameout)" )
end


println("")
println("## Glattiw info ##" )

ImFreqVal, Glatt_iw = GetFreqGreenFromIterId(IterIdFinal)
@show size(Glatt_iw)
@show typeof(Glatt_iw)

println("")
epsilon     = 1e-6
beta_read   =  pi / ImFreqVal[1]
norbtot        = size(Glatt_iw)[2]
density_arr = [ sum( [ giw_i * exp(im*imfreq_i*epsilon) + conj(giw_i) * exp(-im*imfreq_i*epsilon) for (giw_i,imfreq_i) in Iterators.zip(Glatt_iw[:,iorb,iorb],ImFreqVal) ] ) for iorb in 1:norbtot ] / beta_read .+ 0.5
density     = sum(density_arr)
@show beta_read
@show size(density_arr) @show density_arr
println("orbital-resolved density (spin-summed)")
for iorb in 1:div(norbtot,2)
    @show (iorb, sum(density_arr[2*iorb-1:2*iorb]))
end
println("")

fnameout = "Nele.dat"
@show density
open( fnameout, "w" ) do fr
    write( fr, "$(real(density))" )
end
println( "Saved at $(fnameout)" )


fnameout = "DFT_CorrNele.dat"
nspinorb_per_imp = div( norbtot, nimp )
density_imp = [ sum(density_arr[1+nspinorb_per_imp*(iimp-1):nspinorb_per_imp*iimp]) for iimp in 1:nimp ]
open( fnameout, "w" ) do fr
    for iimp in 1:nimp
        @show (iimp, density_imp[iimp])
        write( fr, "latt.  $(iimp) : $(real(density_imp[iimp]))" )
        write( fr, "\n")
    end
end
println( "Saved at $(fnameout)" )

fnameout = "transform_Mat.dat"
norb_per_imp = div(nspinorb_per_imp, 2)
open( fnameout, "w" ) do fr
    for iimp in 1:nimp
        ri = I(norb_per_imp)*1.0
        ii = I(norb_per_imp)*0.0
        m_flat = zeros(norb_per_imp, norb_per_imp*2)
        # @show m_flat[1:2:norb_per_imp,1:2:2*norb_per_imp]
        # @show m_flat[1:2:norb_per_imp,2:2:2*norb_per_imp]
        # @show ri
        # @show ii
        m_flat[1:2:norb_per_imp,1:2:2*norb_per_imp] = ri
        m_flat[1:2:norb_per_imp,2:2:2*norb_per_imp] = ii
        write( fr, "atom [ $(iimp)] \n\n" )
        writedlm( fr, m_flat )
        write( fr, "\n" )
    end
end
println( "Saved at $(fnameout)" )
