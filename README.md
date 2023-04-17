# AdaptiveHypergraphs

🚧 Under construction 🚧

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

Everything should be installed now and you should be ready to go. 

## Reproducing the figures

> **TODO**


## Running the simulation

### From the REPL (recommended)

To run the simulation, execute 

```bash
julia> include("./scripts/main.jl")
```

from the Julia REPL. Alternatively, if you are using VS Code, navigate to the `main.jl` file and run it using the command "Julia: Execute active file in REPL". 

> **Note**
> When you load the code into a REPL for the first time, Julia precompiles the code and all libraries which the code uses. This can take an annoyingly long amount of time, especially at the very first start. Unfortunately, most of this time comes from precompiling the dependent packages, so there is little we can do about it. Here is a [video](https://youtu.be/kmGnutu7_ZY) of adorable manul kittens you can watch while the code compiles. 
>
> To reduce this time, do not kill the Julia REPL; restarting the simulation from the same REPL is vastly faster than the first run. Do not close the window with the dashboard either; subsequent runs of the simulation will reuse it. You can also change the simulation parameters in the input file between simulation runs without restarting the REPL. 

### From the shell

Run

```bash
$ julia ./scripts/main.jl 
```

You can optionally pass the name of an input file as a command-line argument, for example

```bash
$ julia ./scripts/main.jl ./input/small_network.jl
```

When starting the simulation in this way, you can also execute it in parallel. You will need to install an MPI implementation: [OpenMPI](https://www.open-mpi.org/software/ompi/v4.1/) or [MPICH](https://www.mpich.org/downloads/). To start the simulation on, for example, 8 ranks, do

```bash
mpiexec -n 8 julia ./scripts/main.jl ./input/mpi.jl
```

The input file can be different, but it has to set the option `with_mpi = true`. 

The parallelized version does not show an interactive interface since it is meant to be run on a server. Instead, the simulation results are written to a file and can be visualized later.