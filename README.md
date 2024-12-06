# gsPreCICE
 module enabling partitioned multiphysics simulations for G+Smo with other libraries coupled via [`preCICE`](https://precice.org).

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsPreCICE"```|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, macOS|
|Build status|[![ci](https://github.com/gismo/gsPreCICE/actions/workflows/ci.yml/badge.svg)](https://github.com/gismo/gsPreCICE/actions/workflows/ci.yml)|
|Repository|[gismo/gismo/gsPreCICE](https://github.com/gismo/gsPreCICE)|
|Developer|[Hugo Verhelst](https://github.com/hverhelst), [Jingya Li](https://github.com/Crazy-Rich-Meghan) |
|Dependency|[preCICE v.3](https://github.com/gismo/gsPreCICE)|
|Maintainer|[j.li-9@tudelft.nl](mailto:j.li-9@tudelft.nl)|

## Start here

1. Install [preCICE](https://precice.org/quickstart.html).
2. Clone the [preCICE tutorials repository](https://github.com/precice/tutorials).* 
3. Get G+Smo and the adapter. 

The `gsPreCICE` module is a submodule of the G+Smo library. Follow the steps below to download and configure it:
```
git clone https://github.com/gismo/gismo.git
cd gismo
mkdir build
cd build
cmake .. -DGISMO_OPTIONAL="<Other submodules>;gsPreCICE"
```

4. Build examples.
```
make <tutorial file name> -j <number_of_threads>
```
5. Link the compiled executable to the gismo-executable folder within the tutorial directory.*
```
cd <Your preCICE tutorial folder>/partitioned-heat-conduction/gismo-executable
ln -sf <You G+Smo build folder>/bin/<executable name> ./gismo_executable`
```
6. Open two terminals and run.
```
 cd solid-gismo
 ./run.sh
```

```
cd fluid-<other solvers>
./run.sh
```


**Note**: You need to perform steps 2 and 5 if you want to run the simulation with other libraries.

## Examples
- [Partitioned Heat Conduction](examples/partitioned-heat-conduction/README.md)

## Versions

The latest supported G+Smo version is v24.08.0, and the latest supported preCICE version is v3.1.2.

## References

[1] Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. Official preCICE Adapters for Standard Open-Source Solvers. In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017.
