<h1 align="center">TEPPP</h1>
<p align="center">
  <a href="./LICENSE">
    <img alt="GitHub" src="https://img.shields.io/github/license/tomhers/TEPPP">
  </a>
  <img alt="GitHub commit activity" src="https://img.shields.io/github/commit-activity/m/Zer0F8th/TEPPP">
  <img alt="GitHub issues" src="https://img.shields.io/github/issues/Zer0F8th/TEPPP">
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
   - [Linking Numbers](#linking-numbers)
   - [Periodic Writhe](#periodic-writhe)
   - [Jones Polynomial](#jones-polynomial)
   - [Linking Number Between Pairs of Chains](#linking-number-between-pairs-of-chains)
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

> :bulb: Note that only .read_data, .dump, and .dcd files are supported at this time and that all coordinates must be provided to the convertor in unwrapped form. Once the convertor command has been run, the file with the converted data will be located in the `TEPPP/converted` directory for further use.

Once a .teppp file with the desired coordinates has been generated, any of the software commands can be used in conjunction with the file to generate results. The `base` commands that are currently available are:

* jones | Calculates the Jones polynomial of each chain in the system
* lk | Calculates the linking number between each pair of chains in the systemLinking Numbers
* wr | Calculates the Writhe of each chain in the system

All `base` commands are called using the same syntax. There is one required parameter:

1. The filename (including path) of the data file containing the coordinates of the system to analyze.

In addition to these `base` commands, there are several types of variant commands also included in TEPPP. `periodic` commands analyze the topological entanglement of a given system while accounting for periodic boundary conditions. The `periodic` commands that are currently available are:

* periodic_wr | Calculates the periodic Writhe of each chain in the system
* periodic_lk | Calculates the periodic linking number between each pair of chains in the system

The syntax for calling `periodic` commands is the same as the syntax for calling `base` commands; there is one required parameter:

1. The filename (including path) of the data file containing the coordinates of the system to analyze.

`scan` commands are used to analyze the topological entanglement of certain parts of chains rather than the entire chain. For example, if a user wants to the part of a single chain that contributes the most to the overall Writhe of that chain, they would use a `scan` command. The `scan` commands that are currently available are:

* jones_scan | Calculates the Jones polynomial along each chain at given intervals
* lk_scan | Calculates the linking number along each pair of chains at given intervals
* wr_scan | Calculates the Writhe along each chain at given intervals

Calling `scan` commands requires 4 parameters, which must be provided in the command line in the order shown below:

1. The filename (including path) of the data file containing the coordinates of the system to analyze.
2. The length of the initial interval at which to scan.
3. The length of the final interval at which to scan.
4. The amount to increase the interval after a scan completes.

Additionally, there is an optional parameter that can be passed only to the `jones_scan` command. The user may scan for a specific knot type in order to find the location of a knot within a chain by passing the name of the knot to `jones_scan` in the command line. Currently, the only supported knot types are trefoil, figure-8, pentafoil, and stevedore.

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

### Linking Numbers:

To calculate the linking numbers between each pair of chains in a system found in "../data/systemA.teppp" with 100 chains each of length 20:

```bash
./lk "../data/systemA.teppp"
```
### Periodic Writhe:

To calculate the periodic Writhe of each chain in a system found in "../data/systemB.teppp" with 50 chains each of length 25 and a periodic box with length 9.275:

```bash
./periodic_wr "../data/systemB.teppp"
```
### Jones Polynomial:

To scan along each chain and calculate the Jones polynomial of each subset of each chain from size 10 to size 20, skipping by 2, and searching for a trefoil knot in a system found in "../data/systemC.teppp" with 10 chains each of length 200:

```bash
./jones_scan "../data/systemC.teppp" 10 20 2 "trefoil"
```
### Linking Number Between Pairs of Chains:

To calculate the linking number between each pair of chains in a system found in "../data/systemD.teppp" with 20 chains each of length 20, using MPI to split the work between 4 different processes:

```bash
mpirun -np 4 ./lk_mpi "../data/systemD.teppp"
```
## License
<h3><a href="./LICENSE">BSD 3-Clause "New" or "Revised" License</a></h3>

>Copyright (c) 2021, Eleni Panagiotou and Tom Herschberg \
>All rights reserved.

## Contributors
<a href="https://github.com/tomhers/TEPPP/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=tomhers/TEPPP" />
</a>


