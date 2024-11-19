import DFTforge

pack_dir = pathof(DFTforge)

stringind_src = findfirst("src", pack_dir)[end]
pack_src_dir = pack_dir[1:stringind_src]

@show pack_src_dir

println("\n\n")
println("## HOTFIX ## ")
println("Modified src/inputHandler.jl in DFTforge.")
target_file = pack_src_dir*"/inputHandler.jl" 
patch_file  = "hotfix/dftforge/src/inputHandler.jl"
println("Making a backup file.")
run(`chmod +w $(target_file)`)
run(`cp $(target_file) $(target_file).orig.bak`)

println("Copying ...")
println( target_file, " into " , patch_file )
run(`cp $(patch_file) $(target_file)`)

run(`chmod -w $(target_file)`)

println("Done.")
println("\n")
