MFT_dir     = "./"
Magnon_dir  = map(x->string(x), ARGS)[1]

println("")
println( "original dir : ",     MFT_dir )
println( "destination dir : ",  Magnon_dir )
println("")

try 
    readdir( Magnon_dir )
catch
    run(`mkdir -p $(Magnon_dir)` )
    println("Create the specified directory : ", Magnon_dir)
end


cpfdir = string(pwd(),"/", MFT_dir)


println("Copying ...")




# fname = "Jx_DMFT.jl"
# cp(string(cpfdir,"/",fname), string(Magnon_dir,"/",fname),force=true)
# println("cp ", MFT_dir,"/",fname,"  ./",fname)


fname = filter!(x-> occursin(r"toml",x), readdir(MFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(Magnon_dir,"/",fname),force=true)
println("cp ", MFT_dir,"/",fname,"  ",Magnon_dir*"/"*fname)


fname = filter!(x-> occursin(r".win",x), readdir(MFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(Magnon_dir,"/",fname),force=true)
println("cp ", MFT_dir,"/",fname,"  ",Magnon_dir*"/"*fname)


fname = filter!(x-> occursin(r"_hr.dat",x), readdir(MFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(Magnon_dir,"/",fname),force=true)
println("cp ", MFT_dir,"/",fname,"  ",Magnon_dir*"/"*fname)


fname = filter!(x-> occursin(r"occ",x), readdir(MFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(Magnon_dir,"/",fname),force=true)
println("cp ", MFT_dir,"/",fname,"  ",Magnon_dir*"/"*fname)


fname = filter!(x-> occursin(r"raw",x), readdir(MFT_dir))
#println(fname)
for i =1:size(fname)[1]
    println("cp ", MFT_dir,"/",fname[i],"  ",Magnon_dir*"/"*fname[i])
    cp(string(cpfdir,"/",fname[i]), string(Magnon_dir,"/",fname[i]),force=true)
end

println("")
println("-------------------------------------------------------------------------------")
println("To calculate magnon dispersion, change *.toml file 'Calculation_mode = \"Magnon\"")
println("and write a file \"kpath\" in the following format :")
println("")
println("                 G 0.00  0.00  0.00    N 0.00  0.00  0.50")
println("                 N 0.00  0.00  0.50    P 0.25  0.25  0.25")
println("                 P 0.25  0.25  0.25    G 0.00  0.00  0.00")
println("                 G 0.00  0.00  0.00    H 0.50 -0.50  0.50")
println("                 H 0.50 -0.50  0.50    N 0.00  0.00  0.50")
println("-------------------------------------------------------------------------------")
println("")
println("")
