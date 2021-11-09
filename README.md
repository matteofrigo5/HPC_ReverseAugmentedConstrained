## REVERSE AUGMENTED OVERVIEW

NOTE: Before using this driver, the user must install the Chronos package and activate its license.
In order to download and install Chronos, please contact the developer team by using the
form provided at [DEMO_Chronos](https://www.m3eweb.it/chronos/ "DEMO_Chronos")

This package is composed by the following directories:

`Drivers` contains the source codes (c++) of the Driver_Lagrange, the Reverse Augmented Constraint preconditioner and the Saddle Point matrix classes.

`Chronos` is the symbolic link to the INSTALL directory of the the linear solver package (version 1.0).

`ExternalLib` is the symbolic link to the the external libraries necessary for Chronos:

 - `lapack` lapack package (version 3.8.0).

 - `parmetis` metis package (version 4.0.3).

 - `pugixml` pugixml package (version 11.1).

 - `jwt-cpp` A header only library for creating and validating json web tokens in c++ (version 0.5.0). 

 - `curlpp` C++ wrapper around libcURL (version 0.8.1).

 - `LexActivator` interface for Chronos with LexActivator package (https://docs.cryptlex.com/changelog/lexactivator).

Both the INSTALL and the ExternalLib directories are provided in the Chronos installation package.

`Binary` contains scripts (bash) to compile the driver.

`Benchmarks` contains the examples used for the validation of the code.

## TO COMPILE 

### - REQUIRED PACKAGES
cmake 3.10.2, GNU 7.5.0, OMP 4.5, MPI 3.1

### - COMPILE DRIVER
Create directories to store binary file in `Binary` directory:

```
mkdir Driver_Lagrange
```

Run scripts in `Binary` directory:

```
./run_cmake
./run_make
```

### - COMPILE DRIVER AFTER MINOR CHANGES
Run script in the `Binary` directory:

```
./run_make
```

### - COMPILE DRIVER AFTER MAJOR CHANGES
Run scripts in the `Binary` directory:

```
./rm_binary
./run_cmake
./run_make
```

## TO RUN THE BENCHMARK
The executable `driver_Lagrange` is located at the path `./Binary/Driver_Lagrange/src/Core/driver_Lagrange`.

Move to `./Benchmarks/TestLagrange/` directory and extract the matrices from the archive with the following bash command:

```
tar -xvf mat.tar.gz
```

Finally, run the script with `./RUN` in order to test the program.
