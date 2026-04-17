# Sloshing-2D

Real-time 2D sloshing tank simulator solving the incompressible Navier-Stokes equations with PLIC (Piecewise Linear Interface Calculation) free-surface tracking.

![C++20](https://img.shields.io/badge/C%2B%2B-20-blue) ![OpenGL](https://img.shields.io/badge/OpenGL-2.1-green)

## Features

- **MAC staggered grid** eliminates checkerboard pressure oscillations
- **Fractional-step projection** method with semi-Lagrangian advection
- **PLIC geometric VOF** — sharp interface by construction, no numerical diffusion
- **Velocity extrapolation** into the air band (normal-direction weighted) for accurate free-surface dynamics
- **Non-homogeneous Neumann wall BCs** folded into the pressure Poisson equation
- **Red-Black SOR** pressure solver with warm-starting
- **CSF surface tension** model (Brackbill et al.) with SOLID-aware curvature
- **Internal obstacles** (baffles) via persistent obstacle mask
- **9 initial conditions** selectable at runtime, including OpenFOAM damBreak geometry
- **Texture-mapped rendering** via OpenGL for efficient visualization
- **Interactive controls** — move the tank with arrow keys and watch the fluid slosh

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
| 1-9 | Select initial condition |
| T | Toggle surface tension |
| Space | Pause / Resume |
| R | Reset |
| +/- | Speed up / down |
| Q / Esc | Quit |

### Initial Conditions

| Key | Setup |
|-----|-------|
| 1 | Still water |
| 2 | Dambreak |
| 3 | Dambreak with baffle (centre) |
| 4 | Tilted surface (sloshing mode) |
| 5 | Central column |
| 6 | Two-column collision |
| 7 | Shallow fill |
| 8 | Dambreak with baffle (3/4) |
| 9 | OpenFOAM damBreak tutorial geometry |

## Tests

```bash
cd build && ctest
```

Eight physics tests validate the solver against analytical solutions:

| Test | Metric | Value |
|------|--------|-------|
| Hydrostatic pressure | RMS error | 5.1% |
| Sloshing frequency | period error | 0.8% |
| Ritter dambreak wavefront | sim/theory ratio | 0.71 |
| Thin-film spreading | sim/theory ratio | 0.64 |
| Interactive dambreak (256x128) | Ritter ratio at t=0.5s | 0.71 |
| Baffle impermeability | leakage | 0% |
| Volume under violent shaking | drift | 0% |
| VOF boundedness | all cells in [0,1] | pass |

## Project Structure

```
src/
  grid.{hpp,cpp}        MAC staggered grid, field storage, initial conditions
  solver.{hpp,cpp}      Navier-Stokes fractional-step solver with velocity extrapolation
  vof.{hpp,cpp}         PLIC geometric VOF transport
  renderer.{hpp,cpp}    OpenGL texture-based renderer
  simulation.{hpp,cpp}  Main loop orchestrator
tests/
  test_grid.cpp         Grid data structure tests
  test_solver.cpp       Solver unit tests
  test_vof.cpp          VOF transport tests
  test_physics.cpp      Physics validation tests (8 cases)
docs/
  update-plic.md        Development notes on the PLIC migration
```

## License

Free for personal, academic, and research use. No warranty. Commercial use requires written permission from the author.
