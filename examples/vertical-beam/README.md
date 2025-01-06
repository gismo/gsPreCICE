---
title: Perpendicular Flap (with IGA Solid Participant Communicating Stress Data)
keywords: G+Smo, perpendicular flap
summary: This tutorial demonstrates the isogeometric analysis (IGA) solid solver version of the “Perpendicular Flap” tutorial. It focuses on using G+Smo to handle solid-structure interactions by exchanging stress data during simulations.
---


## Overview

This example demostrates how three different coupling approaches which are suitable for IsoGeometric Analaysis based participant.
### Field to field communication
- IGA-IGA coupling method.
- IGA-Spline coupling method.

### Point to point communication
- Vertex-vertex coupling method (exchange quadrature points).

For more details of the coupling methods, the G+Smo-preCICE adapter paper will be available soon.

## Setup
This tutorial uses the similiar setup as the perpendicular flap example, however the fluid participant is communicating with fake data.

## Requirements

To run the tutorial you need to install the following components:
- [preCICE](https://precice.org/quickstart.html)
- [G+Smo and gsPreCICE](https://github.com/gismo/gismo)

## Run the tutorial

In the G+Smo build folder, build the tutorial file.

```
make <participant file> -j <#threads>

make install <participant file>
```


You can find the corresponding `run.sh` script for running the case in the folders corresponding to the participant you want to use:

```bash
cd fake-fluid
./run.sh
```

and

```bash
cd solid
./run.sh
```

## Post-processing
The performance of the three coupling methods are shown in the following figure.
![Couping Methods Performance](https://github.com/Crazy-Rich-Meghan/gsPreCICE/blob/ghpage/examples/communication_time.pdf)

