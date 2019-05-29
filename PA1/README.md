CX4220/CSE 6220 Programming Assignment 1
=================================

## Code structure

Only one file that contains all the code at the root of this project

- `prog1.cpp`: Implement the dart board algorithm to estimate the value of PI

## Compiling

In order to compile everything, simply run
```sh
mpicc prog1.cpp -o prog1
```


## Running
For running on local system, do:
```sh
mpirun -np P ./prog1 N R
```

Where P, N, and R are program arguments that must be set by user.
