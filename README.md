# satutils

A set of utilities common for both the simulation and processing of signals from various satellite systems.

## Prerequisites
First a valid C++20 compiler should be installed. I recommend clang++-18 which is easy on ubuntu 24.04.
```sh
sudo apt install clang-18
```

Make sure Eigen and Pybind11 are installed on your computer.
```sh
sudo apt install libeigen3-dev
sudo apt install python3-pybind11
```

You must also make sure that `navtools` is available for the `satutils` package to find! This most likely implies building from a higher level cmake lists that merely include `navtools` and `satutils` as subdirectories.
```cmake
cmake_minimum_required(VERSION 3.15...3.27)
add_subdirectory(src/navtools)
add_subdirectory(src/satutils)
```

## Building
To build the the standalone C++ project use the build script:
```sh
./build.sh
```
You can turn off the python installation by setting the `build_python` option to `False` in the `build.sh` file.

To build the python project, first create a virtual environment, and then pip install the project (replace the "." with the path to the satutils folder):
```sh
python3 -m venv .venv
. .venv/bin/activate
pip install numpy
pip install .
```

## Python Linting
```sh
pip install pybind11-stubgen
pybind11-stubgen satutils -o src
```
