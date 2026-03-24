/**
 * @file vof.cpp
 * @brief Implementation of VOF transport with Weymouth-Yue operator-split advection,
 *        Superbee TVD flux limiter, and interface compression.
 */

#include "vof.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

/// @brief Superbee flux limiter — the most compressive TVD limiter.
static double superbee(double r) {
    if (r <= 0.0) return 0.0;
    return std::max({0.0, std::min(2.0 * r, 1.0), std::min(r, 2.0)});
}

/// @brief Compute TVD flux through a face using upwind + Superbee correction.
static double tvd_flux(double f_upup, double f_up, double f_down,
                       double vel, double dt, double dh) {
    double courant = std::abs(vel) * dt / dh;
    courant = std::min(courant, 1.0);

    double flux = f_up * courant;

    double df = f_down - f_up;
    if (std::abs(df) > 1e-12) {
        double r = (f_up - f_upup) / df;
        double phi = superbee(r);
        flux += 0.5 * phi * courant * (1.0 - courant) * df;
    }

    return flux;
}

void VOFTransport::advect_x(Grid& grid, double dt) {
    int NX = grid.NX;
    int NY = grid.NY;

    std::vector<double> face_flux((NX + 1) * NY, 0.0);
    auto flux_idx = [&](int i, int j) { return i * NY + j; };

    for (int i = 2; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            double vel = grid.U(i, j);
            if (std::abs(vel) < 1e-14) continue;

            double f_up, f_down, f_upup;
            if (vel > 0.0) {
                f_up   = grid.VOF(i - 1, j);
                f_down = grid.VOF(i, j);
                f_upup = (i >= 3) ? grid.VOF(i - 2, j) : f_up;
            } else {
                f_up   = grid.VOF(i, j);
                f_down = grid.VOF(i - 1, j);
                f_upup = (i < NX - 2) ? grid.VOF(i + 1, j) : f_up;
            }

            face_flux[flux_idx(i, j)] = tvd_flux(f_upup, f_up, f_down,
                                                   vel, dt, grid.dx);
        }
    }

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            if (grid.Type(i, j) == CellType::SOLID) continue;

            double F_L = face_flux[flux_idx(i, j)];
            double F_R = face_flux[flux_idx(i + 1, j)];

            double u_L = grid.U(i, j);
            double u_R = grid.U(i + 1, j);

            double gain_left  = (u_L > 0.0) ?  F_L : -F_L;
            double lose_right = (u_R > 0.0) ?  F_R : -F_R;

            double stretch = 1.0 - dt / grid.dx * (u_R - u_L);
            stretch = std::max(stretch, 0.1);

            double f_new = (grid.VOF(i, j) + gain_left - lose_right) / stretch;
            grid.VOF(i, j) = std::clamp(f_new, 0.0, 1.0);
        }
    }
}

void VOFTransport::advect_y(Grid& grid, double dt) {
    int NX = grid.NX;
    int NY = grid.NY;

    std::vector<double> face_flux(NX * (NY + 1), 0.0);
    auto flux_idx = [&](int i, int j) { return i * (NY + 1) + j; };

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 2; j < NY - 1; ++j) {
            double vel = grid.V(i, j);
            if (std::abs(vel) < 1e-14) continue;

            double f_up, f_down, f_upup;
            if (vel > 0.0) {
                f_up   = grid.VOF(i, j - 1);
                f_down = grid.VOF(i, j);
                f_upup = (j >= 3) ? grid.VOF(i, j - 2) : f_up;
            } else {
                f_up   = grid.VOF(i, j);
                f_down = grid.VOF(i, j - 1);
                f_upup = (j < NY - 2) ? grid.VOF(i, j + 1) : f_up;
            }

            face_flux[flux_idx(i, j)] = tvd_flux(f_upup, f_up, f_down,
                                                   vel, dt, grid.dy);
        }
    }

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            if (grid.Type(i, j) == CellType::SOLID) continue;

            double F_B = face_flux[flux_idx(i, j)];
            double F_T = face_flux[flux_idx(i, j + 1)];

            double v_B = grid.V(i, j);
            double v_T = grid.V(i, j + 1);

            double gain_bottom = (v_B > 0.0) ?  F_B : -F_B;
            double lose_top    = (v_T > 0.0) ?  F_T : -F_T;

            double stretch = 1.0 - dt / grid.dy * (v_T - v_B);
            stretch = std::max(stretch, 0.1);

            double f_new = (grid.VOF(i, j) + gain_bottom - lose_top) / stretch;
            grid.VOF(i, j) = std::clamp(f_new, 0.0, 1.0);
        }
    }
}

void VOFTransport::sharpen_interface(Grid& grid) {
    // Gentle interface sharpening: snap very-nearly-full and very-nearly-empty
    // cells to their binary value. This prevents slow accumulation of diffuse
    // VOF "haze" without creating staircase artifacts.
    // Cells in the active transition zone (0.02 < f < 0.98) are left untouched.
    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            double f = grid.VOF(i, j);
            if (f > 0.0 && f < 0.02) grid.VOF(i, j) = 0.0;
            if (f < 1.0 && f > 0.98) grid.VOF(i, j) = 1.0;
        }
    }
}

void VOFTransport::enforce_vof_boundaries(Grid& grid) {
    // Zero VOF in all boundary (SOLID) cells
    for (int i = 0; i < grid.NX; ++i) {
        grid.VOF(i, 0) = 0.0;
        grid.VOF(i, grid.NY - 1) = 0.0;
    }
    for (int j = 0; j < grid.NY; ++j) {
        grid.VOF(0, j) = 0.0;
        grid.VOF(grid.NX - 1, j) = 0.0;
    }
}

void VOFTransport::set_reference_volume(const Grid& grid) {
    reference_volume_ = total_volume(grid);
}

void VOFTransport::advect(Grid& grid, double dt) {
    if (reference_volume_ < 0.0) {
        set_reference_volume(grid);
    }

    // Alternate sweep order for symmetry
    if (step_count_ % 2 == 0) {
        advect_x(grid, dt);
        advect_y(grid, dt);
    } else {
        advect_y(grid, dt);
        advect_x(grid, dt);
    }
    ++step_count_;

    sharpen_interface(grid);
    enforce_vof_boundaries(grid);
    correct_mass(grid);
}

void VOFTransport::correct_mass(Grid& grid) {
    double current = total_volume(grid);
    if (current < 1e-14) return;

    double ratio = reference_volume_ / current;

    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            grid.VOF(i, j) = std::clamp(grid.VOF(i, j) * ratio, 0.0, 1.0);
        }
    }
}

double VOFTransport::total_volume(const Grid& grid) const {
    double vol = 0.0;
    double cell_area = grid.dx * grid.dy;
    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            vol += grid.VOF(i, j) * cell_area;
        }
    }
    return vol;
}
