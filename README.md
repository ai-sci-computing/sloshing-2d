# Sloshing-2D

Real-time 2D sloshing tank simulator solving the incompressible Navier-Stokes equations with VOF (Volume of Fluid) free-surface tracking.

![C++20](https://img.shields.io/badge/C%2B%2B-20-blue) ![OpenGL](https://img.shields.io/badge/OpenGL-2.1-green)

## Features

- **MAC staggered grid** eliminates checkerboard pressure oscillations
- **Fractional-step projection** method with semi-Lagrangian advection
- **VOF transport** with Superbee TVD limiter and Weymouth-Yue divergence correction
- **Red-Black SOR** pressure solver with warm-starting
- **CSF surface tension** model (Brackbill et al.)
- **Texture-mapped rendering** via OpenGL for efficient visualization
- **Interactive controls** — move the tank with arrow keys and watch the fluid slosh
- **Physics-validated** against analytical hydrostatic and sloshing frequency solutions

## Building

Requires CMake 3.20+, a C++20 compiler, GLFW 3.3+, and doctest.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Running

```bash
./build/sloshing_tank           # Default 256x128 grid
./build/sloshing_tank 512 256   # Custom resolution
```

### Controls

| Key | Action |
|-----|--------|
| Arrow keys / WASD | Accelerate tank |
| Space | Pause / Resume |
| R | Reset |
| +/- | Speed up / down |
| Q / Esc | Quit |

## Tests

```bash
cd build && ctest
```

Four test suites validate grid operations, solver correctness, VOF transport, and physics:

- **Hydrostatic pressure**: RMS error ~5% vs analytical p = rho g (H - y)
- **Sloshing frequency**: 3.7% error vs linear wave theory at 256x128

## Project Structure

```
src/
  grid.{hpp,cpp}        MAC staggered grid and field storage
  solver.{hpp,cpp}      Navier-Stokes fractional-step solver
  vof.{hpp,cpp}         VOF advection with TVD flux limiting
  renderer.{hpp,cpp}    OpenGL texture-based renderer
  simulation.{hpp,cpp}  Main loop orchestrator
tests/
  test_grid.cpp         Grid data structure tests
  test_solver.cpp       Solver unit tests
  test_vof.cpp          VOF transport tests
  test_physics.cpp      Physics validation tests
```

## License

Free for personal, academic, and research use. No warranty. Commercial use requires written permission from the author.
