![GitHub commits since latest release](https://img.shields.io/github/commits-since/gismo/gsKLShell/latest?color=008A00)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/gismo/gsKLShell?color=008A00)

# gsPreCICE

Module enabling partitioned multiphysics simulations for G+Smo with other libraries coupled via [`preCICE`](https://precice.org).

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsPreCICE"```|
|--:|---|
|License|![GitHub License](https://img.shields.io/github/license/gismo/gismo?color=008A00)|
|OS support|Linux, Windows, macOS|
|Build status|[![ci](https://github.com/gismo/gsPreCICE/actions/workflows/ci.yml/badge.svg)](https://github.com/gismo/gsPreCICE/actions/workflows/ci.yml)|
|Developers/maintainers| [![Static Badge](https://img.shields.io/badge/@Crazy--Rich--Meghan-008A00)](https://github.com/Crazy-Rich-Meghan) [![Static Badge](https://img.shields.io/badge/@hverhelst-008A00)](https://github.com/hverhelst)|
|Dependency|[preCICE v.3](https://github.com/gismo/gsPreCICE)|


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

4. Build examples and make the examples discoverable in the system.
```
make <tutorial file name> -j <number_of_threads>
make install <tutorial file name>

```
5. Open two terminals and run.
```
 cd solid-gismo
 ./run.sh
```

```
cd fluid-<other solvers>
./run.sh
```


## Examples
- [Partitioned Heat Conduction](examples/partitioned-heat-conduction/README.md)
- [Perpendicular Flap](https://github.com/gismo/gsPreCICE/tree/main/examples/perpendicular-flap/README.md)
  
  **Note:** To run these two examples, `gsElasticity`, `gsKLShell` and `gsStructuralAnalysis` are also needed.
  ```
  cmake .. -DGISMO_OPTIONAL="gsKLShell;gsPreCICE;gsElasticity;gsStructuralAnalysis"
  ```

## Versions

The submodule is up to date with the recent G+Smo release , and the latest supported preCICE version is v3.1.2.

## References

[1] Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. Official preCICE Adapters for Standard Open-Source Solvers. In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017.
