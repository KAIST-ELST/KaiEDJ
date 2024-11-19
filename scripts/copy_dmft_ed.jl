DMFT_dir    = "./"
MFT_dir     = map(x->string(x), ARGS)[1]

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


println("")
println("Converting the file format of self-energy ")
println("")
include("../src/convert_info.jl")
dmftfiles   = readdir( DMFT_dir )

println("")
println("Copying files from \"$(DMFT_dir)\" to \"$(MFT_dir)/\"")
println("")

fnames      = filter( x-> occursin("SelfE_",x) , dmftfiles )
map( x->cp( joinpath(DMFT_dir,x), joinpath(MFT_dir,x), force=true), fnames  )
println("Copied (forced by default) : ", join(fnames," "))

fnames      = filter( x-> occursin("toml",x) , dmftfiles )
map( x->cp( joinpath(DMFT_dir,x), joinpath(MFT_dir,x), force=true), fnames  )
println("Copied (forced by default) : ", join(fnames," "))

fnames      = filter( x-> occursin("wannier",x) , dmftfiles )
map( x->cp( joinpath(DMFT_dir,x), joinpath(MFT_dir,x), force=true,follow_symlinks=true), fnames  )
println("Copied (forced by default) : ", join(fnames," "))

fnames      = filter( x-> occursin("Nele",x) , dmftfiles )
map( x->cp( joinpath(DMFT_dir,x), joinpath(MFT_dir,x), force=true,follow_symlinks=true), fnames  )
println("Copied (forced by default) : ", join(fnames," "))

fnames      = filter( x-> occursin("transform_Mat",x) , dmftfiles )
map( x->cp( joinpath(DMFT_dir,x), joinpath(MFT_dir,x), force=true,follow_symlinks=true), fnames  )
println("Copied (forced by default) : ", join(fnames," "))


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
println("                    DMFT+MFT   : Calculation with local self-energy")
println("                    MFT        : Calculation with local self-energy = 0.0")
println("-------------------------------------------------------------------------------")
