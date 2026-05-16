This simulator currently allows uncontrolled, unpowered kinematic simulations of miscellaneous objects on both a flat and spherical Earth.
Flat Earth will orient using the standard NED-cartesian coordinate system while spherical Earth uses standard Earth-Centered Earth-Fixed coordinate system.
Models currently available for use are the NASA standard spheroid, brick and ISS used for check cases.
Things not currently modeled:
- Earth's Rotation
- Differentiation between rocket shapes and aircraft shapes
- Advanced simulation of aerodynamic surfaces
- Powered Flight
- Controlled Flight

## Dependencies

- CMake 3.16+
- A C++26-compatible compiler (GCC 15+ or Clang 21+)
- clangd 21+ (optional, for IDE support)
- clang-format / clang-tidy (optional, for formatting and static analysis)

### Linux (Ubuntu)
```bash
sudo apt install cmake g++ clangd clang-format clang-tidy
```

### Windows

Install [Visual Studio](https://visualstudio.microsoft.com/) with the **Desktop development with C++** workload — CMake support is included.

## Building

### Linux
```bash
cmake -B build
cmake --build build
./build/QuaternionSim
```

### Windows
Open the project folder directly in Visual Studio (`File > Open > Folder`). It will detect `CMakeLists.txt` and configure automatically. Use the build button or `Ctrl+Shift+B` to build.
