# To run the DMFT, modify your TOML script (here, LCO_kaiedj.toml) and execute the following command.
# (Example) $ julia --project=./KaiEDJ/env/julia-1.10.3_kaiedj  -p 6 ./KaiEDJ/kaiedj.jl  -T LCO_kaiedj.toml   2>&1 | tee out
julia --project=<julia_env_path>  -p 6 <KaiEDJ_path>/KaiEDJ/kaiedj.jl  -T LCO_kaiedj.toml   2>&1 | tee out

# To run the MFT+DMFT, modify your TOML script (here, LCO_kaiedj.toml) and execute the following command in order.

# 1) Copy the DMFT output into a new MFT+DMFT directory.
julia --project=<julia_env_path>  <KaiEDJ_path>/KaiEDJ/scripts/prepare_mft_from_dmft.jl  <DMFT_path> <MFT+DMFT_path>

# 2) The MFT+DMFT will use the DMFT output as an input.
julia --project=<julia_env_path>  -p 6 <KaiEDJ_path>/KaiEDJ/kaiedj.jl  -T LCO_kaiedj.toml   2>&1 | tee out


