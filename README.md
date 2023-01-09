# hydro_tests

## Purpose

This is a collection of all test cases used to explore the capabilities of
MFM.

## Usage

A run can be executed using run_test_problem with the options described within
the file (run run_hydro_test -h to see the help).

Parameter files are converted between different codes and brough to the same
format using standardize_param.py.
A more detailed description on the usage can be found within that file.

## Structure

The main directory contains a few files that are useful to setup and run a problem.
The individual problems are located in sub-folders.

```
├── blob
|
├── grav_freefall
|
├── hydrostatic_sphere
|
├── hydrostatic_square
|
├── kepler_disk
|
├── kh
|
├── nifty_cluster
|
├── rt
|
├── sedov_blast
|
├── shock_tube
|
├── soundwave
|
├── turbulence
|
└── zeldovich_pancake

Besides parameter and Config files, these directories contain routines to create
ICs (if applicable) and plotting scripts used for the plots in the paper.
