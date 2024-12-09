
## Setup

We solve a partitioned heat equation. For information on the non-partitioned case, please refer to [1, p.37ff]. In this tutorial the computational domain is partitioned and coupled via preCICE. The coupling roughly follows the approach described in [2].

![Case setup of partitioned-heat-conduction case](https://github.com/precice/tutorials/blob/master/partitioned-heat-conduction/images/tutorials-partitioned-heat-conduction-setup.png)

Case setup from [2]. `D` denotes the Dirichlet participant and `N` denotes the Neumann participant.

The heat equation is solved on a rectangular domain `Omega = [0,2] x [0,1]` with given Dirichlet boundary conditions. We split the domain at `x_c = 1` using a straight vertical line, the coupling interface. The left part of the domain will be referred to as the Dirichlet partition and the right part as the Neumann partition. To couple the two participants we use Dirichlet-Neumann coupling. Here, the Dirichlet participant receives Dirichlet boundary conditions (`Temperature`) at the coupling interface and solves the heat equation using these boundary conditions on the left part of the domain. Then the Dirichlet participant computes the resulting heat flux (`Flux`) from the solution and sends it to the Neumann participant. The Neumann participant uses the flux as a Neumann boundary condition to solve the heat equation on the right part of the domain. We then extract the temperature from the solution and send it back to the Dirichlet participant. This establishes the coupling between the two participants.

This simple case allows us to compare the solution for the partitioned case to a known analytical solution (method of manufactures solutions, see [1, p.37ff]). For more usage examples and details, please refer to [3, sect. 4.1].

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](https://github.com/precice/tutorials/blob/master/partitioned-heat-conduction/images/tutorials-partitioned-heat-conduction-precice-config.png)


## Running the simulation

In the G+Smo build folder, build the tutorial file.

```
make <participant file> -j <#threads>
```

Go to the gismo-executable folder and link the compiled executable to the gismo_executable.

```
cd gismo-executable
chmod +x create_symlink.sh
./create_symlink.sh
```

You can find the corresponding `run.sh` script for running the case in the folders corresponding to the participant you want to use:

```bash
cd dirichlet-gismo
./run.sh
```

and

```bash
cd neumann-gismo
./run.sh
```



## Visualization


For G+Smo, please use the file `solution.pvd` in both dirichlet-gismo and neumann-gismo directories. 

![Animation of the partitioned heat equation](https://github.com/Crazy-Rich-Meghan/tutorials/blob/partitioned-heat-conduction/partitioned-heat-conduction/images/tutorial-partitioned-heat-conduction-gismo.gif)

<!-- Visualization in paraview for `x_c = 1.5`. -->

## References

[1] Azahar Monge and Philipp Birken. "Convergence Analysis of the Dirichlet-Neumann Iteration for Finite Element Discretizations." (2016). Proceedings in Applied Mathematics and Mechanics. [doi](https://doi.org/10.1002/pamm.201610355)  
[2] Benjamin Rüth, Benjamin Uekermann, Miriam Mehl, Philipp Birken, Azahar Monge, and Hans Joachim Bungartz. "Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications." (2020). International Journal for Numerical Methods in Engineering. [doi](https://doi.org/10.1002/nme.6443)  
[3] Tobias Eppacher. "Parallel-in-Time Integration with preCICE" (2024). Bachlor's thesis at Technical University of Munich. [pdf](https://mediatum.ub.tum.de/1755012)
