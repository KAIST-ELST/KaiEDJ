# KaiEDJ
*KaiEDJ: Integrated DMFT-ED and MFT Framework for Correlated Magnetic Systems*

**[Contact]**

* Hyeong Jun Lee   --  hjuntaf (at) gmail.com
* Taekjung Kim     --  tj.kim (at) kaist.ac.kr
* Hongkee yoon     --  hongkeeyoon (at) kangwon.ac.kr

The KaiEDJ is maintained by KAIST-ELST Group (https://github.com/KAIST-ELST).

**[Main Features]**
* Dynamical mean-field theory (DMFT) calculation
  * Exact diagonalization (ED) solver
  * Support an interface to external solver packages
  * DFT+DMFT interface - compatible with multiple density functional theory (DFT) codes
* Magnetic force theorem (MFT) or magnetic force response theory calculation with correlated electronic structure.
* Spin-wave dispersion analysis

**[Citation]**

(Submitted) 
*KaiEDJ: A program conducting dynamical mean-field theory and magnetic force theory calculation for correlated magnetic materials*,
Hyeong Jun Lee, Taek Jung Kim, Hongkee Yoon, Myung Joon Han


## 0. Preparation

- Download `KaiEDJ`. 
  ```
  git clone https://github.com/KAIST-ELST/KaiEDJ.git
  ```

- `Julia` should be installed. ([Download Julia here.](https://julialang.org/downloads/ "official webpage")) Julia doesn't require any configurations or compile-processes in the beginning. You will just download it and set its `bin/` path to your `PATH` to use `julia` excutable.
  ```
  $ export PATH=<your-julia-full-path>/bin/:$PATH    # a shell-command example.
  export PATH=<your-julia-full-path>/bin/:$PATH      # a shell-script example such as in .bashrc or .bash_profile
  ```
  >  Version-test successfully done for
  > * julia >= 1.8.5,  <= 1.11.3


- Install a set of required packages supported by Julia by running `setup.jl` provided in our package.
  ```
  $ julia setup.jl
  ```
  This will install the following packages.
  `Formatting` `LinearAlgebra` `SparseArrays` `DelimitedFiles` `FastGaussQuadrature` `BenchmarkTools` `Optimization` `OptimizationOptimJL` `Arpack` `ThreadedSparseCSR` `SparseMatricesCSR` `Plots`
  `TickTock` `TimerOutputs` `FFTW` `JSON` `Dierckx` `ImageFiltering` `DFTforge`

- Export your `PROJECT_PATH_KAIEDJ` to your project-environment directory. (`setup.jl` will automatically generate the project-environment directory of `envs/KEDJ`.)
  ```
  $ export PROJECT_PATH_KAIEDJ=<your-KaiEDJ-full-path>/envs/KEDJ    # shell-command
  export PROJECT_PATH_KAIEDJ=<your-KaiEDJ-full-path>/envs/KEDJ      # shell-script such as .bashrc or .bash_profile
  ```

## 1. Files

- `run.sh` is a single running script for all processes including the installation and all of the example codes.
- `setup.jl` will install the required packages, activating a Julia-project-environment TOML directory.
- `kaiedj.jl` is a main script for running KaiEDJ within Julia.
- `eddmft.jl` is a script for running ED-DMFT within Julia.
- `dmftmft.jl` is a script for running DMFT+MFT within Julia.
- `examples/` contains example codes of La2CuO4 for ED-DMFT, DMFT+MFT, Spinwave, and of Fe for (QMC-DMFT,) DMFT+MFT, Spinwave calculations.
- `src/` contains source codes.
- `envs/` is a temporary directory for Julia-project-environment TOML files. `setup.jl` will use this location and generate `envs/KED/*.toml` files.
- `scripts/` contains scripts for data-converting and patches.
- `hotfix/` contains a hotfix files.

## 2. Usage

Preparing the TOML script, e.g. `example.TOML`, you can run KaiEDJ as follow.

```
$ julia kaiedj.jl -T example.TOML
```
In order to activate the Julia-project-environment, you can add a flag as `julia --project=<julia-project-environment path> ...`.

For tutorials, you can download in [github.com/hjuntaf/KaiED_tutorials.git](https://github.com/hjuntaf/KaiED_tutorials.git) and run as follows.

```
$ julia tutorials/01_manybody_binary_basis.jl
```
Every tutorial code here contains `include("../src/mybase.jl")` with relative path to include `mybase.jl` and use `KaiED`. You will need to modify this path (`../src/`) when you want to run in a different directory/location.


### Jupyter notebook

We also provide some tutorials within `jupyter notebook` (.ipynb).
You can use the jupyter notebooks after an installation of `IJulia` as follows.

#### Installation

```
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> using Pkg

julia> Pkg.add("IJulia")
```
or enter the Pkg REPL by pressing `]` and install it :
```
(@v1.9) pkg> add IJulia
```

#### Run

```
$ jupyter notebook
```

## 3. Tutorials

We provide some tutorial codes in formats of .jl (julia source code) and .ipynb (jupyter notebook).
You can see and download them in the following repo.: [github.com/hjuntaf/KaiED_tutorials](https://github.com/hjuntaf/KaiED_tutorials)

Following basic concepts are included.
- Quantum many-body basis state manipulations
- Quantum many-body wavefunction manipulations
- Green functions
- Impurity models
- Tight-binding models
- Exact diagonalization method
- Dynamical mean-field theory


## 4. (Optional) External Packages

* Our interface supports DMFT calculations with external DMFT solver packages.
* Currently, it is compatible with two continuous-time quantum Monte Carlo (CT-QMC) solvers.
  * CT-QMC of COMSUITE package (https://github.com/comscope/comsuite)
  * CT-QMC of EDMFTF (http://hauleweb.rutgers.edu/tutorials/)

    ### CT-QMC of EDMFTF package installation
    Below is a very brief introduction to the installation of EDMFTF package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the EDMFTF package site directly (http://hauleweb.rutgers.edu/tutorials/).
    
    Download and extract compressed EDMFTF package
    ```
    $ wget http://hauleweb.rutgers.edu/downloads/EDMFTF.tgz
    $ tar -zxvf EDMFTF.tgz
    ```
    
    Modify file configure.py to suit your computer's compiler environment
    ```
    $ vim EDMFTF-*/configure.py
    ```
    
    ,and run
    ```
    $ python setup.py
    ```
    
    If at least the following files are created in ``/bini`` folder, it is possible to run Jx_DMFT.
    ```
    ../bini/ctqmc
    ../bini/atom_d.py
    ../bini/maxent_run.py
    ```

    Finally, add the following lines to ``$ vi ~/.bashrc``
    ```
    export WIEN_DMFT_ROOT=EDMFTF-installed folder/bini
    export PYTHONPATH=$PYTHONPATH:$WIEN_DMFT_ROOT
    export SCRATCH="."
    export PATH=$WIEN_DMFT_ROOT:$PATH
    ```
    

    ### CT-QMC(ComCTQMC) of COMSUITE package installation
    Below is a very brief introduction to the installation of the COMSUITE package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the COMSUITE package site directly (https://github.com/comscope/comsuite).
    
    Download zip file from https://github.com/comscope/comsuite and unzip the compressed folder
    ```
    $ unzip comsuite-master.zip 
    ```
    
    Modify file arch.mk to suit your computer's compiler environment
    ```
    $ vim comsuite-master/arch.mk
    ```
    
    and, add the following two lines to ``$ vi ~/.bashrc``
    
    ```
    export COMSUITE_BIN= comsuite-installed folder/comsuite-master/bin  
    export PATH=$COMSUITE_BIN:$PATH
    ```
   
    Finally, 
    ```
    $ make all
    ```
    
    If the following files are created in ``/bin`` folder, you are ready to run Jx_DMFT.
    
    ```
    ./comsuite-master/bin/CTQMC  
    ./comsuite-master/bin/EVALSIM 
    ```


    
