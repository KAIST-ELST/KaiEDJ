DMFT_dir    = "./"
MFT_dir     = map(x->string(x), ARGS)[1]
SCF_DMFT_dir= DMFT_dir

println("")
println( "original dir : ",  DMFT_dir )
println( "destination dir : ",  MFT_dir )
println("")

try 
    readdir( MFT_dir )
catch
    run(`mkdir -p $(MFT_dir)` )
    println("Create the specified directory : ", MFT_dir)
end

cpfdir = string(pwd(),"/", SCF_DMFT_dir)

println("Copying ...")



for i = 1:10
    impdir = string(cpfdir,"/imp_",i)
    if isdir(impdir)

        pwdimpdir = string( MFT_dir, "/imp_",i)
        if !isdir(pwdimpdir)
            mkdir(pwdimpdir)
        end
 
        fname = "params.obs.json"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)

            fnlist = filter!(x-> occursin(r"params.obs",x), readdir(impdir))
            
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end         
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/params.obs*  ./imp_",i,"/params.obs*")

        end

        fname = "Sig.out"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)


            fnlist = filter!(x-> occursin(r"Sig",x), readdir(impdir))                
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/Sig*  ./imp_",i,"/Sig*")

        end

        fname = "Gf.out"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)


            fnlist = filter!(x-> occursin(r"Gf",x), readdir(impdir))                 
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/Gf*  ./imp_",i,"/Gf*")

        end

    end
end


# fname = "Jx_DMFT.jl"
# cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
# println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)




fname = filter!(x-> occursin(r"toml",x), readdir(SCF_DMFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier.win"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier_hr.dat"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "Nele.dat"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)

fname = "DFT_CorrNele.dat"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)

fname = "transform_Mat.dat"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "scf.log"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)

fname = "occ.log"
cp(string(cpfdir,"/",fname), string(MFT_dir,"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)



println("")
println("")
println("-------------")
println("## NOTICE ##")
println("-------------")
println("You MUST modify your .toml file to run the MFT calculation as follows.")
println("")
println("[KaiEDJ]")
println("Calculation_mode = \"DMFT+MFT\"    instead of \"DMFT\"" )
println("")
println("")
println("-------------------------------------------------------------------------------")
println("To calculate strength of exchange interaction, change *.toml file 'Calculation_mode = \"DMFT+MFT\" or \"MFT\" ")
println("")
println("                    DMFT+MFT   : Calculation with local self-energy")
println("                    MFT        : Calculation with local self-energy = 0.0")
println("-------------------------------------------------------------------------------")
println("")
println("")



