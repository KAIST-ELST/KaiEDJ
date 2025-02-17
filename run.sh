# Running KaiEDJ 

# 0. Installation
export JULIA="julia +1.10"
$JULIA setup.jl
    
WORKDIR=$PWD
JULIA_PROJECT_PATH="./envs/KEDJ"
# JULIA_PROJECT_PATH="~/envs_julia/dmftjx_julia-1.10.3"



# 1. La2CuO4 Example
cd examples/La2CuO4


# 1-1. La2CuO4 Example: ED-DMFT

# 1-1-a. Running ED-DMFT
echo "[info] 1-1-a Running ED-DMFT"
cd ED_DMFT_SCF/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/kaiedj.jl  -T LCO_kaiedj_eddmft.toml 2>&1 | tee log
cd ../../
echo "[info] 1-1-a Running ED-DMFT done"

# 1-1-b. Plotting and saving the ED-DMFT spectral-function result in output_dmft_lco.png
echo "[info] 1-1-b Plotting and saving the ED-DMFT spectral-function result in output_dmft_lco.png"
cd ED_DMFT_SCF/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/plot_spectral_functions.jl
open   output_dmft_lco.png || evince output_dmft_lco.pdf
cd ../../
echo "[info] 1-1-b Plotting and saving the ED-DMFT done"



# 1-2. La2CuO4 Example: DMFT+MFT

# 1-2-a. Preparation for DMFT+MFT (Files are already prepared in this file-set)
echo "[info] 1-2-a Preparation for DMFT+MFT"
cd ED_DMFT_SCF/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/prepare_mft_from_dmft.jl ../../DMFT_MFT/in
cd ../../
echo "[info] 1-2-a Preparation for DMFT+MFT done"

# 1-2-b. Edit LCO_kaiedj.toml : Calculation_mode="DMFT+MFT"
echo "[info] 1-2-b Edit LCO_kaiedj.toml : Calculation_mode=DMFT+MFT"
cp $WORKDIR/scripts/LCO_kaiedj_dmftmft.toml DMFT_MFT/in/LCO_kaiedj_dmftmft.toml
echo "[info] 1-2-b Edit LCO_kaiedj.toml : Calculation_mode=DMFT+MFT done"

# 1-2-c. Running DMFT+MFT 
    
echo "[info] 1-2-c Running DMFT+MFT"
cd DMFT_MFT/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/kaiedj.jl  -T LCO_kaiedj_dmftmft.toml 2>&1 | tee log
cd ../../
echo "[info] 1-2-c Running DMFT+MFT done"

# 1-2-d. Plotting and saving the DMFT+MFT result in output_mft_lco.png
echo "[info] 1-2-d Plotting and saving the DMFT+MFT result in output_mft_lco.png"
cd DMFT_MFT/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/plot_j_lco.jl
open   output_dmftmft_lco.png || evince output_dmftmft_lco.pdf
cd ../../
echo "[info] 1-2-d Plotting and saving the DMFT+MFT done"



# 1-3. La2CuO4 Example: Spinwave

# 1-3-a. Preparation for Spinwave (Files are already prepared in this file-set)
echo "[info] 1-3-a Preparation for Spinwave"
cd DMFT_MFT/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/prepare_magnon.jl ../../SPINWAVE/in
cd ../../
echo "[info] 1-3-a Preparation for Spinwave done"

# 1-3-b. Edit LCO_kaiedj.toml : Calculation_mode="Spinwave"
echo "[info] 1-3-b Edit LCO_kaiedj.toml : Calculation_mode=Spinwave"
cp $WORKDIR/scripts/LCO_kaiedj_spinwave.toml SPINWAVE/in/LCO_kaiedj_spinwave.toml

# 1-3-c. Edit kpath for the spinwave dispersion
echo "[info] 1-3-c Edit kpath for the spinwave dispersion"
cp $WORKDIR/scripts/kpath_spinwave_lco SPINWAVE/in/kpath

# 1-3-d. Running Spinwave calculation
echo "[info] 1-3-d Running Spinwave calculation"
cd SPINWAVE/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/kaiedj.jl  -T LCO_kaiedj_spinwave.toml 2>&1 | tee log
cd ../../
echo "[info] 1-3-d Running Spinwave calculation done"

# 1-3-e. Plotting and saving the Spinwave result in output_spinwave_lco.png
echo "[info] 1-3-e Plotting and saving the Spinwave result in output_spinwave_lco.png"
cd SPINWAVE/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/plot_spinwave_lco.jl
open   output_spinwave_lco.png || evince output_spinwave_lco.pdf
cd ../../
echo "[info] 1-3-e Plotting and saving the Spinwave done"



# 2. Fe-bcc Example
echo "[info] 2. Fe-bcc Example"
cd $WORKDIR/examples/Fe/



# 2-1. Fe-bcc Example: QMC-DMFT (Results already obtained externally from eDMFT are prepared in this file-set)
# (Omitted)



# 2-2. Fe-bcc Example: DMFT+MFT

# 2-2-a. Preparation for DMFT+MFT (Files are already prepared in this file-set)
echo "[info] 2-2-a Preparation for DMFT+MFT"
cd QMC_DMFT_SCF/out/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/prepare_mft_from_dmft.jl ../../DMFT_MFT/in
cd ../../
echo "[info] 2-2-a Preparation for DMFT+MFT done"

# 2-2-b. Edit .toml : Calculation_mode="DMFT+MFT"
cp $WORKDIR/scripts/fe_kaiedj_dmftmft.toml DMFT_MFT/in/fe_kaiedj_dmftmft.toml

# 2-2-c. Running DMFT+MFT 
echo "[info] 2-2-c Running DMFT+MFT"
cd DMFT_MFT/in
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/kaiedj.jl  -T fe_kaiedj_dmftmft.toml 2>&1 | tee log
cd ../../
echo "[info] 2-2-c Running DMFT+MFT done"

# 2-2-d. Plotting and saving the DMFT+MFT result in output_mft_fe.png
echo "[info] 2-2-d Plotting and saving the DMFT+MFT result in output_mft_fe.png"
cd DMFT_MFT/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/plot_j_fe.jl
echo "[info] 2-2-d Plotting and saving the DMFT+MFT done"
#open   output_dmftmft_fe.png || evince output_dmftmft_fe.pdf
cd ../../



# 2-3. Fe-bcc Example: Spinwave

# 2-3-a. Preparation for Spinwave (Files are already prepared in this file-set)
echo "[info] 2-3-a Preparation for Spinwave"
cd DMFT_MFT/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/prepare_magnon.jl ../../SPINWAVE/in
cd ../../
echo "[info] 2-3-a Preparation for Spinwave done"

# 2-3-b. Edit LCO_kaiedj.toml : Calculation_mode="Spinwave"
cp $WORKDIR/scripts/fe_kaiedj_spinwave.toml SPINWAVE/in/fe_kaiedj_spinwave.toml

# 2-3-c. Edit kpath for the spinwave dispersion
cp $WORKDIR/scripts/kpath_fe SPINWAVE/in/kpath

# 2-3-d. Running Spinwave calculation
echo "[info] 2-3-d Running Spinwave calculation"
cd SPINWAVE/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/kaiedj.jl  -T fe_kaiedj_spinwave.toml 2>&1 | tee log
cd ../../
echo "[info] 2-3-d Running Spinwave calculation done"

# 2-3-e. Plotting and saving the Spinwave result in output_spinwave_fe.png
echo "[info] 2-3-e Plotting and saving the Spinwave result in output_spinwave_fe.png"
cd SPINWAVE/in/
$JULIA --project=$JULIA_PROJECT_PATH $WORKDIR/scripts/plot_spinwave_fe.jl
open   output_spinwave_fe.png || evince output_spinwave_fe.pdf
cd ../../

echo "[info] 2. Fe-bcc Example done"