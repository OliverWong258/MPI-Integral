# Project Description

This project requires several libraries for compilation and execution, including LAPACK, ScaLAPACK, OpenMPI, and BLAS. The program is designed to run in parallel using MPI and OpenMP.

## Environment Setup

Before compiling and running the project, you need to install the following dependencies:
- **LAPACK**
- **ScaLAPACK**
- **OpenMPI**
- **BLAS**

Please make sure these libraries are installed and properly configured on your system.

## Compilation and Execution

### 1. Clean the Build
To remove any previously compiled executable files, run:

```bash
make clean
```

### 2. Compile the Project
Once you've cleaned the previous build, run the following command to compile the project:

```bash
make
```

If the compilation is successful, the executable file `final_openmp` will be generated.

### 3. Run the Program
To run the program, use the following command:

```bash
mpirun -np x ./final_openmp
```

Where `x` is the number of processes to use. Typically, using 4 processes is sufficient.

### 4. Input File
When running the program, you will be prompted to provide the path to the `INPUT.txt` file. You can either:
- Input the full path to the file, or
- Type `default` to use the default path `../input/INPUT.txt`.

### 5. Parallel Threads
After specifying the input file, you will be prompted to enter the number of parallel threads. Enter the desired number based on the resources available.

## Documentation

For more detailed instructions and explanations, please refer to the `Report.docx` file.