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
      cell_type(nx * ny, CellType::EMPTY),
      is_obstacle(nx * ny, 0)
{
}

void Grid::classify_cells(double threshold) {
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            if (is_obstacle[idx(i, j)]) {
                Type(i, j) = CellType::SOLID;
            } else if (VOF(i, j) > threshold) {
                Type(i, j) = CellType::FLUID;
            } else {
                Type(i, j) = CellType::EMPTY;
            }
        }
    }
}

void Grid::initialize_water(double fill_fraction) {
    std::fill(is_obstacle.begin(), is_obstacle.end(), 0);

    double water_height = fill_fraction * Ly;

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double y_bottom = j * dy;
            double y_top    = (j + 1) * dy;

            if (y_top <= water_height) {
                VOF(i, j) = 1.0;
            } else if (y_bottom >= water_height) {
                VOF(i, j) = 0.0;
            } else {
                VOF(i, j) = (water_height - y_bottom) / dy;
            }
        }
    }

    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);

    classify_cells();
}

void Grid::initialize_dambreak(double column_width_frac, double column_height_frac) {
    std::fill(is_obstacle.begin(), is_obstacle.end(), 0);
    std::fill(vof.begin(), vof.end(), 0.0);

    double col_w = std::clamp(column_width_frac,  0.05, 0.95) * Lx;
    double col_h = std::clamp(column_height_frac, 0.05, 0.95) * Ly;

    for (int i = 0; i < NX; ++i) {
        double x_left  = i * dx;
        double x_right = (i + 1) * dx;
        for (int j = 0; j < NY; ++j) {
            double y_bottom = j * dy;
            double y_top    = (j + 1) * dy;

            if (x_left >= col_w || y_bottom >= col_h) { VOF(i, j) = 0.0; continue; }

            double frac_x = (std::min(x_right, col_w) - x_left) / dx;
            double frac_y = (std::min(y_top,    col_h) - y_bottom) / dy;
            VOF(i, j) = std::clamp(frac_x * frac_y, 0.0, 1.0);
        }
    }

    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);

    classify_cells();
}

void Grid::initialize_dambreak_baffle(double column_width_frac,
                                      double column_height_frac,
                                      double baffle_x_frac,
                                      double baffle_height_frac,
                                      int baffle_thickness_cells) {
    initialize_dambreak(column_width_frac, column_height_frac);

    int i_center = std::clamp(static_cast<int>(baffle_x_frac * NX), 1, NX - 2);
    int half = std::max(1, baffle_thickness_cells) / 2;
    int i_lo = std::max(0, i_center - half);
    int i_hi = std::min(NX - 1, i_center + (baffle_thickness_cells - half - 1));
    int j_top = std::clamp(static_cast<int>(baffle_height_frac * NY), 1, NY - 1);

    for (int i = i_lo; i <= i_hi; ++i) {
        for (int j = 0; j < j_top; ++j) {
            is_obstacle[idx(i, j)] = 1;
            VOF(i, j) = 0.0;
        }
    }

    classify_cells();
}

void Grid::initialize_tilted_surface(double fill_fraction, double amplitude) {
    std::fill(is_obstacle.begin(), is_obstacle.end(), 0);

    double mean_h = fill_fraction * Ly;

    for (int i = 0; i < NX; ++i) {
        double x_center = (i + 0.5) * dx;
        double local_h = mean_h + amplitude * (x_center - Lx / 2.0) / (Lx / 2.0);
        local_h = std::max(local_h, 0.0);

        for (int j = 0; j < NY; ++j) {
            double y_bottom = j * dy;
            double y_top = (j + 1) * dy;

            if (y_top <= local_h) VOF(i, j) = 1.0;
            else if (y_bottom >= local_h) VOF(i, j) = 0.0;
            else VOF(i, j) = (local_h - y_bottom) / dy;
        }
    }

    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);
    classify_cells();
}

void Grid::initialize_central_column(double column_width_frac, double column_height_frac) {
    std::fill(is_obstacle.begin(), is_obstacle.end(), 0);
    std::fill(vof.begin(), vof.end(), 0.0);

    double col_w = std::clamp(column_width_frac, 0.05, 0.95) * Lx;
    double col_h = std::clamp(column_height_frac, 0.05, 0.95) * Ly;
    double x_lo = (Lx - col_w) / 2.0;
    double x_hi = (Lx + col_w) / 2.0;

    for (int i = 0; i < NX; ++i) {
        double xl = i * dx, xr = (i + 1) * dx;
        if (xl >= x_hi || xr <= x_lo) continue;
        double fx = (std::min(xr, x_hi) - std::max(xl, x_lo)) / dx;

        for (int j = 0; j < NY; ++j) {
            double yb = j * dy, yt = (j + 1) * dy;
            if (yb >= col_h) continue;
            double fy = (std::min(yt, col_h) - yb) / dy;
            VOF(i, j) = std::clamp(fx * fy, 0.0, 1.0);
        }
    }

    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);
    classify_cells();
}

void Grid::initialize_two_columns(double column_width_frac, double column_height_frac,
                                  double gap_frac) {
    std::fill(is_obstacle.begin(), is_obstacle.end(), 0);
    std::fill(vof.begin(), vof.end(), 0.0);

    double col_w = std::clamp(column_width_frac, 0.05, 0.4) * Lx;
    double col_h = std::clamp(column_height_frac, 0.05, 0.95) * Ly;
    double gap = std::clamp(gap_frac, 0.1, 0.8) * Lx;
    double x_left_hi = (Lx - gap) / 2.0;
    double x_left_lo = x_left_hi - col_w;
    double x_right_lo = (Lx + gap) / 2.0;
    double x_right_hi = x_right_lo + col_w;

    for (int i = 0; i < NX; ++i) {
        double xl = i * dx, xr = (i + 1) * dx;

        double fx = 0.0;
        if (xr > x_left_lo && xl < x_left_hi)
            fx = (std::min(xr, x_left_hi) - std::max(xl, x_left_lo)) / dx;
        else if (xr > x_right_lo && xl < x_right_hi)
            fx = (std::min(xr, x_right_hi) - std::max(xl, x_right_lo)) / dx;
        if (fx <= 0.0) continue;

        for (int j = 0; j < NY; ++j) {
            double yb = j * dy, yt = (j + 1) * dy;
            if (yb >= col_h) continue;
            double fy = (std::min(yt, col_h) - yb) / dy;
            VOF(i, j) = std::clamp(fx * fy, 0.0, 1.0);
        }
    }

    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);
    classify_cells();
}

void Grid::initialize_shallow(double fill_fraction) {
    initialize_water(std::clamp(fill_fraction, 0.02, 0.25));
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
