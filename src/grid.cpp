/**
 * @file grid.cpp
 * @brief Implementation of the MAC staggered grid.
 */

#include "grid.hpp"
#include <algorithm>
#include <cmath>

Grid::Grid(int nx, int ny, double lx, double ly)
    : NX(nx), NY(ny), dx(lx / nx), dy(ly / ny), Lx(lx), Ly(ly),
      u((nx + 1) * ny, 0.0),
      v(nx * (ny + 1), 0.0),
      p(nx * ny, 0.0),
      vof(nx * ny, 0.0),
      cell_type(nx * ny, CellType::EMPTY)
{
}

void Grid::classify_cells(double threshold) {
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            // Boundary cells are SOLID walls
            if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1) {
                Type(i, j) = CellType::SOLID;
            } else if (VOF(i, j) > threshold) {
                Type(i, j) = CellType::FLUID;
            } else {
                Type(i, j) = CellType::EMPTY;
            }
        }
    }
}

void Grid::enforce_boundary_conditions() {
    // Left and right walls: u = 0 at wall faces
    for (int j = 0; j < NY; ++j) {
        U(0, j) = 0.0;  // Left wall
        U(1, j) = 0.0;  // Inner face of left SOLID layer
        U(NX - 1, j) = 0.0; // Inner face of right SOLID layer
        U(NX, j) = 0.0; // Right wall
    }

    // Bottom and top walls: v = 0 at wall faces
    for (int i = 0; i < NX; ++i) {
        V(i, 0) = 0.0;  // Bottom wall
        V(i, 1) = 0.0;  // Inner face of bottom SOLID layer
        V(i, NY - 1) = 0.0; // Inner face of top SOLID layer
        V(i, NY) = 0.0; // Top wall
    }

    // Free-slip tangential BC: reflect velocity into solid cells
    // Bottom/top: u in solid row mirrors the adjacent fluid row
    for (int i = 1; i < NX; ++i) {
        // j=0 is solid: tangential u mirrors from j=1
        U(i, 0) = U(i, 1);
        // j=NY-1 is solid: tangential u mirrors from j=NY-2
        if (NY >= 2) U(i, NY - 1) = U(i, NY - 2);
    }

    // Left/right: v in solid column mirrors the adjacent fluid column
    for (int j = 1; j < NY; ++j) {
        // i=0 is solid: tangential v mirrors from i=1
        V(0, j) = V(1, j);
        // i=NX-1 is solid: tangential v mirrors from i=NX-2
        if (NX >= 2) V(NX - 1, j) = V(NX - 2, j);
    }
}

void Grid::initialize_water(double fill_fraction) {
    double water_height = fill_fraction * Ly;

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            // Boundary cells are SOLID walls — no fluid
            if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1) {
                VOF(i, j) = 0.0;
                continue;
            }

            double y_bottom = j * dy;
            double y_top = (j + 1) * dy;

            if (y_top <= water_height) {
                VOF(i, j) = 1.0;
            } else if (y_bottom >= water_height) {
                VOF(i, j) = 0.0;
            } else {
                VOF(i, j) = (water_height - y_bottom) / dy;
            }
        }
    }

    // Zero velocity everywhere
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);

    classify_cells();
}

std::pair<double, double> Grid::interpolate_velocity(double x, double y) const {
    // Bilinear interpolation of u-velocity (defined at vertical faces)
    // u(i,j) lives at (i*dx, (j+0.5)*dy)
    double u_x = x / dx;
    double u_y = y / dy - 0.5;

    int i0 = std::max(0, std::min(NX, static_cast<int>(std::floor(u_x))));
    int j0 = std::max(0, std::min(NY - 2, static_cast<int>(std::floor(u_y))));
    int i1 = std::min(NX, i0 + 1);
    int j1 = std::min(NY - 1, j0 + 1);

    double sx = u_x - i0;
    double sy = u_y - j0;
    sx = std::clamp(sx, 0.0, 1.0);
    sy = std::clamp(sy, 0.0, 1.0);

    double u_val = (1 - sx) * (1 - sy) * U(i0, j0)
                 + sx       * (1 - sy) * U(i1, j0)
                 + (1 - sx) * sy       * U(i0, j1)
                 + sx       * sy       * U(i1, j1);

    // Bilinear interpolation of v-velocity (defined at horizontal faces)
    // v(i,j) lives at ((i+0.5)*dx, j*dy)
    double v_x = x / dx - 0.5;
    double v_y = y / dy;

    i0 = std::max(0, std::min(NX - 2, static_cast<int>(std::floor(v_x))));
    j0 = std::max(0, std::min(NY, static_cast<int>(std::floor(v_y))));
    i1 = std::min(NX - 1, i0 + 1);
    j1 = std::min(NY, j0 + 1);

    sx = v_x - i0;
    sy = v_y - j0;
    sx = std::clamp(sx, 0.0, 1.0);
    sy = std::clamp(sy, 0.0, 1.0);

    double v_val = (1 - sx) * (1 - sy) * V(i0, j0)
                 + sx       * (1 - sy) * V(i1, j0)
                 + (1 - sx) * sy       * V(i0, j1)
                 + sx       * sy       * V(i1, j1);

    return {u_val, v_val};
}
