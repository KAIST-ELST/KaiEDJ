using TOML

# Specify the path to your TOML file
FilePathTOML   = "ex.toml"

# Read the TOML file
TOMLHeader = TOML.parsefile(FilePathTOML)

# Access values from the parsed TOML data
value1 = TOMLHeader["section"]["key1"]
value2 = TOMLHeader["section"]["key2"]

# Print or use the values as needed
println("Value 1: $value1")
println("Value 2: $value2")


