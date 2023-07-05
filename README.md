# AdaptiveHypergraphs

An interactive simulation of the adaptive voter model on hypergraphs. 

This repository contains the code and figures of the paper "Polyadic Opinion Formation: The Adaptive Voter Model on a Hypergraph". 

## Installation instructions

To run the code locally, you will need Julia (v1.8 was used during development). To install Julia, please refer to the [official documentation](https://julialang.org/downloads/).

If you are using VS Code, you will also need the official [Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia). 

After you installed Julia, start the Julia REPL, activate the virtual environment and install all packages. 

```bash
$ julia
julia> ]activate .
julia> ]instantiate 
```

If you haven't used Julia before, you might wonder about the `]` character. Pressing `]` activates Pkg, Julia's package manager. Once you are inside Pkg, the the text of the prompt should change to:

```bash
(@v1.8) pkg>
```

where the version number in the brackets can be different depending on your Julia version. `@v1.8` is the default virtual environment. The `activate .` command activates the `AdaptiveHypergraphs` environment; after you execute this command, your prompt should look like this:

```bash
(AdaptiveHypergraphs) pkg>
```

Now, execute the `instantiate` command (you don't need to press `]` again if you are still in the Pkg prompt) to install the packages. Once everything is installed, press backspace or Ctrl-C to exit Pkg and return to the Julia REPL. 

If the instantiation progress throws any errors from the `PyCall` package, try setting the `PYTHON` environment variable to an empty string in your Julia REPL:

```julia
julia> ENV["PYTHON"] = ""
```

And call `instantiate` again. 

## Reproducing the figures

The scripts to plot the figures are located in `./scripts/figure_plots`. To run a script, use

```bash
julia> include("./scripts/figure_plots/<figure_name>.jl")
```

and replace `<figure_name>` by the corresponding file name. 

## Running the simulation

### From the REPL (recommended)

To run the simulation, execute 

```bash
julia> include("./scripts/main.jl")
```

from the Julia REPL. Alternatively, if you are using VS Code, navigate to the `./scripts/main.jl` file and run it using the command "Julia: Execute active file in REPL". 

> **Note**: Long compile time at first start
>
> When you load the code into a REPL for the first time, Julia precompiles the code and all libraries which the code uses. This can take an annoyingly long amount of time, especially at the very first start. Unfortunately, most of this time comes from precompiling the dependent packages, so there is little we can do about it.
>
> To reduce this time, do not kill the Julia REPL; restarting the simulation from the same REPL is vastly faster than the first run. Do not close the window with the dashboard either; subsequent runs of the simulation will reuse it. You can also change the simulation parameters in the input file between simulation runs without restarting the REPL. 

### From the shell

This method is primarily used to run the simulation in parallel on a headless terminal without an interactive visualization. If you do not want to run large-scale parallelized simulations, it is highly recommended to use the previous method.

Run

```bash
$ julia --project="." -- ./scripts/main.jl 
```

You can optionally pass the name of an input file as a command-line argument. The path of the input file has to been given relatively to the location of the `main.jl` file. 

```bash
$ julia --project="." -- ./scripts/main.jl ../input/small_network.jl
```

To run the simulation in parallel, you will need to install an MPI implementation: [OpenMPI](https://www.open-mpi.org/software/ompi/v4.1/) or [MPICH](https://www.mpich.org/downloads/). To start the simulation on, for example, 8 ranks, do

```bash
mpiexec -n 8 julia --project="." -- ./scripts/main.jl ../input/mpi.jl
```

The input file can be different, but it has to set the option `with_mpi = true`. 

The parallelized version does not show an interactive interface since it is meant to be run on a server. Instead, the simulation results are written to a file and can be visualized later.

### Using input files

The parameters of the simulation are set using the input files in the `./input` folder. The file `./input/default_network.jl` is used per default if no input file is given. It is convenient to use this file as a sandbox to tweak parameters without creating a new file. 

To pass a different input file, either pass the filename as a command-line argument when starting Julia like this:

```bash
$ julia --project="." -- ./scripts/main.jl ../input/small_network.jl
```

or change the value of the variable `input_file` in `/scripts/main.jl`. This second approach is more convenient, because it can be done without restarting the REPL. 

The whole list of configuration settings and their descriptions can be found in `./src/simulation/InputParams.jl`.