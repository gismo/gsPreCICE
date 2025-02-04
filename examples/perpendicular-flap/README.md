---
title: Perpendicular Flap (with IGA Solid Participant Communicating Stress Data)
keywords: G+Smo, perpendicular flap
summary: This tutorial demonstrates the isogeometric analysis (IGA) solid solver version of the “Perpendicular Flap” tutorial. It focuses on using G+Smo to handle solid-structure interactions by exchanging stress data during simulations.
---


## Overview

This example demostrates how the Geometry + Simulation Modules (G+Smo) can be utilised through a module-type adapter to couple with other codes using preCICE. G+Smo offers a robust framework for isogeometric analysis (IGA), seamlessly integrating geometric representations with numerical solvers to enable advanced simulations and efficient code coupling.

## Setup
This tutorial uses the same setup as the **original** perpendicular flap example, with the key difference being the switch in the communicated information to stress. 

## Requirements

To run the tutorial you need to install the following components:
- [preCICE](https://precice.org/quickstart.html)
- [G+Smo and gsPreCICE](https://github.com/gismo/gismo)
- [OpenFOAM](https://openfoam.org/download/)

## Run the tutorial

In the G+Smo build folder, build the tutorial file.

```
make <participant file> -j <#threads> // In this case: make perpendicular-flap-vertex-gismo -j <#threads>
make install <participant file> 
```


You can find the corresponding `run.sh` script for running the case in the folders corresponding to the participant you want to use:

for the fluid participant with OpenFOAM, 
```bash
cd fluid-openfoam
OpenFOAM<version> // Enter the OpenFOAM environment
./run.sh
```

and

```bash
cd solid-gismo
./run.sh
```

## Post-processing
The results of this tutorial are comparable to the simulation results communicated with force under the perpendicular-flap tutorials.
![G+Smo stress](https://github.com/Crazy-Rich-Meghan/tutorials/blob/perpendicular-flap-gismo-elasticity-stress/perpendicular-flap-stress/images/tutorials-perpendicular-flap-stress-displacement-watchpoint.png)

<!-- Additionally, the mesh convergence study data is available under the `images/data` directory. A sample plot illustrating the convergence of x-displacement over time is shown below:
![G+Smo converfence](https://github.com/Crazy-Rich-Meghan/tutorials/blob/perpendicular-flap-gismo-elasticity-stress/perpendicular-flap-stress/images/x_displacement_vs_time.png) -->
