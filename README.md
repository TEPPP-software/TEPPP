<h1 align="center">TEPPP</h1>
<p align="center">
  <a href="./LICENSE">
    <img alt="GitHub" src="https://img.shields.io/github/license/tomhers/TEPPP">
  </a>
  <img alt="GitHub commit activity" src="https://img.shields.io/github/commit-activity/m/TEPPP-software/TEPPP">
  <img alt="GitHub issues" src="https://img.shields.io/github/issues/TEPPP-software/TEPPP">
</p>

<h3 align="center">The Topological Entanglement of Polymers, Proteins, and Periodic systems (TEPPP) Software</h3>
<hr />

>© Copyright 2021 Tom Herschberg, Kyle Pifer and Eleni Panagiotou \
> \
>If you use this code you must reference the following paper: \
> \
>Herschberg, T., Pifer, K. and Panagiotou, E., \
>A computational package for measuring Topological Entanglement in Polymers, Proteins and Periodic systems (TEPPP), 2021,(submitted)
> \
> \
>visit www.elenipanagiotou.com for updated information

## Table of Contents
- [What is TEPPP?](#what-is-teppp)
- [Getting Started](#getting-started)
   - [Building](#building)
   - [Usage](#usage)
- [Examples](#examples)
   - [Gauss linking integral](#gauss-linking-integral)
   - [Periodic Writhe](#periodic-writhe)
   - [Jones Polynomial](#jones-polynomial)
   - [Scan Jones Polynomial](#scan-jones-polynomial)
- [License](#license)
- [Contributors](#contributors)

<hr />

## What is TEPPP?
<p>
  TEPPP is a software package designed to aid in the calculation of several topological entanglement values in molecular systems. It is designed to work both in serial and in parallel when coupled with a functioning MPI installation.
</p>

## Getting Started

### Building

TEPPP requires a compiler which fully supports C++17 (GCC 7 and up). If parallel execution is desired, a functioning MPI installation must be present in the PATH variable. 

To build only the serial version of the software, run the following command: 
```bash
make serial
``` 
To build only the parallel version of the software, run the following command: 
```bash
make mpi
``` 
To make both versions of the software, type the following command: 
```bash
make all
```

### Usage

The current version of TEPPP only supports running individual commands through the command line to obtain the desired results. If the coordinates of the system to be analyzed are in a file with an extension other than .teppp, they must be converted into a file readable by TEPPP. To do this, run the following command in the top-level installation directory: 
```bash
./convertor "/path/to/filename.ext" CHAIN_LENGTH NUM_CHAINS BOX_DIM
``` 
where `CHAIN_LENGTH` is the number of atoms in each chain, `NUM_CHAINS` is the number of chains in the system, and `BOX_DIM` is the length of one side of the periodic box if the system uses periodic boundary conditions. If the system does not use periodic boundary conditions, enter 0 for `BOX_DIM`. 

> :bulb: Note that only .teppp files are supported. All other data files must be converted to .teppp either by the user or by using the convertor utility. Only .read_data files (with all coordinates in unwrapped form) are supported by the convertor utility at this time. Once the convertor command has been run, the file with the converted data will be located in the `TEPPP/converted` directory for further use.

Once a .teppp file with the desired coordinates has been generated, any of the software commands can be used in conjunction with the file to generate results. The `base` commands that are currently available are:

* jones | Calculates the Jones polynomial of each chain in the system
* lk | Calculates the linking number between each pair of chains in the systemLinking Numbers
* wr | Calculates the Writhe of each chain in the system

All `base` commands are called using the same syntax:

The filename (including path) of the data file containing the coordinates of the system to analyze followed by CHAIN_LENGTH NUM_CHAINS BOX_DIM

In addition to these `base` commands, there are several types of variant commands also included in TEPPP. `periodic` commands analyze the topological entanglement of a given system while accounting for periodic boundary conditions. The `periodic` commands that are currently available are:

* periodic_wr | Calculates the periodic Writhe of each chain in the system
* periodic_lk | Calculates the periodic linking number between each pair of chains in the system

The syntax for calling `periodic` commands is the same as the syntax for calling `base` commands.

`scan` commands are used to analyze the topological entanglement of certain parts of chains rather than the entire chain. For example, if a user wants to the part of a single chain that contributes the most to the overall Writhe of that chain, they would use a `scan` command. The `scan` commands that are currently available are:

* jones_scan | Calculates the Jones polynomial along each chain at given intervals
* lk_scan | Calculates the linking number along each pair of chains at given intervals
* wr_scan | Calculates the Writhe along each chain at given intervals

Calling `scan` commands requires 4 parameters, which must be provided in the command line in the order shown below:

1. The filename (including path) of the data file containing the coordinates of the system to analyze.
2. The length of the initial interval at which to scan.
3. The length of the final interval at which to scan.
4. The amount to increase the interval after a scan completes.

`mpi` commands are parallel versions of the `base`, `periodic`, and `scan` commands discussed above. They leverage MPI to split the workload between a given number of processors rather than performing the work serially. The `mpi` commands that are currently available are:

* jones_mpi | Calculates the Jones polynomial of each chain in the system in parallel
* lk_mpi | Calculates the linking number between each pair of chains in the system in parallel
* wr_mpi | Calculates the Writhe of each chain in the system in parallel
* periodic_wr_mpi | Calculates the periodic Writhe of each chain in the system in parallel
* periodic_lk_mpi | Calculates the periodic linking number between each pair of chains in the system in parallel
* jones_scan_mpi | Calculates the Jones polynomial along each chain at given intervals in parallel
* lk_scan_mpi | Calculates the linking number along each pair of chains at given intervals in parallel
* wr_scan_mpi | Calculates the Writhe along each chain at given intervals in parallel

`mpi` commands have the same syntax as their `base`, `periodic`, and `scan` counterparts but must be called using `mpirun` rather than running the command itself.

## Examples

### Gauss linking integral:

To calculate the Gauss linking integral between each pair of chains in a system found in "../data/systemA.teppp" with 100 chains each of length 20 in a cubic periodic box of length 13.35315:

```bash
./lk "../data/systemA.teppp" 20 100 13.35315
```
using MPI to split the work between 4 different processes:

```bash
mpirun -np 4 ./lk_mpi "../data/systemA.teppp" 20 100 13.35315
```

### Periodic Writhe:

To calculate the Periodic Writhe of each chain in a system found in "../data/systemA.teppp" with 100 chains each of length 20 in a cubic periodic box of length 13.35315:

```bash
./periodic_wr "../data/systemA.teppp" 20 100 13.35315
```
### Jones Polynomial:

To compute the Jones polynomial of each chain in a system found in "../data/systemA.teppp" with 100 chains each of length 20 in a cubic periodic box of length 13.35315:

```bash
./jones "../data/systemA.teppp" 20 100 13.35315
```

### Scan Jones Polynomial:

To scan along each chain and calculate the Jones polynomial of each subchain in a system found in "../data/systemA.teppp" with 100 chains each of length 20 starting with scanning length 5 up to scanning length 10 with a step of 5:

```bash
./jones_scan "../data/systemC.teppp" 20 100 5 10 5
```


## License
<h3><a href="./LICENSE">BSD 3-Clause "New" or "Revised" License</a></h3>

>Copyright (c) 2021, Tom Herschberg, Kyle Pifer and Eleni Panagiotou \
>All rights reserved.

## Contributors
<a href="https://github.com/TEPPP-software/TEPPP/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=tomhers/TEPPP" />
  <img src="https://contrib.rocks/image?repo=TEPPP-software/TEPPP" />
</a>


