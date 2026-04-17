/**
 * @file vof.cpp
 * @brief VOF transport with PLIC (Piecewise Linear Interface Calculation).
 */

#include "vof.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

// Fraction of the unit square [0,1]^2 satisfying  mx*X + my*Y <= d.
// Requires mx >= 0, my >= 0.
static double plic_area_unit(double mx, double my, double d) {
    if (d <= 0.0) return 0.0;
    if (d >= mx + my) return 1.0;
    if (my < 1e-14) return std::clamp(d / mx, 0.0, 1.0);
    if (mx < 1e-14) return std::clamp(d / my, 0.0, 1.0);

    double xi = d / mx, yi = d / my;
    if (xi <= 1.0 && yi <= 1.0) return 0.5 * xi * yi;
    if (xi >  1.0 && yi <= 1.0) return 0.5 * (yi + (d - mx) / my);
    if (xi <= 1.0 && yi >  1.0) return 0.5 * (xi + (d - my) / mx);
    double gap = mx + my - d;
    return 1.0 - 0.5 * gap * gap / (mx * my);
}

// Fluid area in the rectangle [x0,x1] x [y0,y1] for line  nx*x + ny*y <= d
// with arbitrary sign of nx, ny.
static double fluid_area_in_rect(double nx, double ny, double d,
                                  double x0, double x1, double y0, double y1) {
    double sx = x1 - x0, sy = y1 - y0;
    if (sx <= 0.0 || sy <= 0.0) return 0.0;
    double d2 = d - nx * x0 - ny * y0;
    double mx = nx * sx, my = ny * sy;
    if (mx < 0.0) { d2 -= mx; mx = -mx; }
    if (my < 0.0) { d2 -= my; my = -my; }
    return sx * sy * plic_area_unit(mx, my, d2);
}

void VOFTransport::advect_plic(Grid& grid, double dt) {
    int NX = grid.NX, NY = grid.NY;
    double dx = grid.dx, dy = grid.dy;
    double cell_area = dx * dy;

    auto vof_s = [&](int i, int j) -> double {
        i = std::clamp(i, 0, NX - 1);
        j = std::clamp(j, 0, NY - 1);
        if (grid.Type(i, j) == CellType::SOLID) return 0.5;
        return grid.VOF(i, j);
    };

    // --- 1. Reconstruct interface in each partial cell ---
    struct Rec { double nx, ny, d; };
    std::vector<Rec>  rec(NX * NY, {0, 0, 0});
    std::vector<char> is_if(NX * NY, 0);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double f = grid.VOF(i, j);
            if (f <= 0.0 || f >= 1.0) continue;

            // Youngs' gradient (3x3 weighted stencil)
            double gx = -(vof_s(i+1,j-1) + 2.0*vof_s(i+1,j) + vof_s(i+1,j+1)
                         -vof_s(i-1,j-1) - 2.0*vof_s(i-1,j) - vof_s(i-1,j+1)) / (8.0*dx);
            double gy = -(vof_s(i-1,j+1) + 2.0*vof_s(i,j+1) + vof_s(i+1,j+1)
                         -vof_s(i-1,j-1) - 2.0*vof_s(i,j-1) - vof_s(i+1,j-1)) / (8.0*dy);

            double mag = std::sqrt(gx*gx + gy*gy);
            if (mag < 1e-14) { gx = 0.0; gy = -1.0; }
            else { gx /= mag; gy /= mag; }

            // Bisect for line constant d such that the fluid area in the
            // cell-local rectangle [0,dx]x[0,dy] equals f * cell_area.
            double d_lo = std::min({0.0, gx*dx, gy*dy, gx*dx+gy*dy}) - 1e-10;
            double d_hi = std::max({0.0, gx*dx, gy*dy, gx*dx+gy*dy}) + 1e-10;
            for (int it = 0; it < 40; ++it) {
                double dm = 0.5 * (d_lo + d_hi);
                if (fluid_area_in_rect(gx, gy, dm, 0, dx, 0, dy) < f * cell_area)
                    d_lo = dm;
                else
                    d_hi = dm;
            }
            int idx = grid.idx(i, j);
            rec[idx] = {gx, gy, 0.5*(d_lo + d_hi)};
            is_if[idx] = 1;
        }
    }

    // --- 2. Geometric fluxes (unsigned fluid area swept through each face) ---
    std::vector<double> flux_x((NX+1)*NY, 0.0);
    std::vector<double> flux_y(NX*(NY+1), 0.0);
    auto fxi = [&](int i, int j) { return i*NY + j; };
    auto fyi = [&](int i, int j) { return i*(NY+1) + j; };

    for (int i = 1; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double vel = grid.U(i, j);
            if (std::abs(vel) < 1e-14) continue;
            int iD = (vel > 0) ? i-1 : i;
            double w = std::min(std::abs(vel) * dt, dx);
            double fD = grid.VOF(iD, j);

            if      (fD <= 0.0) {}
            else if (fD >= 1.0) { flux_x[fxi(i,j)] = w * dy; }
            else if (is_if[grid.idx(iD,j)]) {
                auto& r = rec[grid.idx(iD,j)];
                double sx0 = (vel > 0) ? dx - w : 0.0;
                double sx1 = (vel > 0) ? dx     : w;
                flux_x[fxi(i,j)] = fluid_area_in_rect(r.nx, r.ny, r.d, sx0, sx1, 0, dy);
            }
            else { flux_x[fxi(i,j)] = fD * w * dy; }
        }
    }

    for (int i = 0; i < NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            double vel = grid.V(i, j);
            if (std::abs(vel) < 1e-14) continue;
            int jD = (vel > 0) ? j-1 : j;
            double w = std::min(std::abs(vel) * dt, dy);
            double fD = grid.VOF(i, jD);

            if      (fD <= 0.0) {}
            else if (fD >= 1.0) { flux_y[fyi(i,j)] = dx * w; }
            else if (is_if[grid.idx(i,jD)]) {
                auto& r = rec[grid.idx(i,jD)];
                double sy0 = (vel > 0) ? dy - w : 0.0;
                double sy1 = (vel > 0) ? dy     : w;
                flux_y[fyi(i,j)] = fluid_area_in_rect(r.nx, r.ny, r.d, 0, dx, sy0, sy1);
            }
            else { flux_y[fyi(i,j)] = fD * dx * w; }
        }
    }

    // --- 3. Conservative update ---
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            if (grid.Type(i,j) == CellType::SOLID) continue;

            double uL = grid.U(i,j),   uR = grid.U(i+1,j);
            double vB = grid.V(i,j),   vT = grid.V(i,j+1);

            double gain = 0.0, lose = 0.0;
            if (uL > 0) gain += flux_x[fxi(i,j)];   else lose += flux_x[fxi(i,j)];
            if (uR < 0) gain += flux_x[fxi(i+1,j)]; else lose += flux_x[fxi(i+1,j)];
            if (vB > 0) gain += flux_y[fyi(i,j)];   else lose += flux_y[fyi(i,j)];
            if (vT < 0) gain += flux_y[fyi(i,j+1)]; else lose += flux_y[fyi(i,j+1)];

            grid.VOF(i,j) = std::clamp(grid.VOF(i,j) + (gain - lose) / cell_area, 0.0, 1.0);
        }
    }
}


// =========================================================================
//  Post-processing and main entry point
// =========================================================================

void VOFTransport::sharpen_interface(Grid& grid) {
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j) {
            double f = grid.VOF(i, j);
            if (f > 0.0 && f < 0.02) grid.VOF(i, j) = 0.0;
            if (f < 1.0 && f > 0.98) grid.VOF(i, j) = 1.0;
        }
}

void VOFTransport::enforce_vof_boundaries(Grid& grid) {
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            if (grid.Type(i, j) == CellType::SOLID) grid.VOF(i, j) = 0.0;
}

void VOFTransport::set_reference_volume(const Grid& grid) {
    reference_volume_ = total_volume(grid);
}

void VOFTransport::advect(Grid& grid, double dt) {
    if (reference_volume_ < 0.0) set_reference_volume(grid);

    advect_plic(grid, dt);

    sharpen_interface(grid);
    enforce_vof_boundaries(grid);
    correct_mass(grid);
}

void VOFTransport::correct_mass(Grid& grid) {
    double current = total_volume(grid);
    if (current < 1e-14) return;

    double cell_area = grid.dx * grid.dy;
    double volume_error = current - reference_volume_;

    double interface_volume = 0.0;
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j) {
            if (grid.Type(i, j) == CellType::SOLID) continue;
            double f = grid.VOF(i, j);
            if (f > 0.0 && f < 1.0) interface_volume += f * cell_area;
        }

    if (interface_volume < 1e-14) {
        double ratio = reference_volume_ / current;
        for (int i = 0; i < grid.NX; ++i)
            for (int j = 0; j < grid.NY; ++j) {
                if (grid.Type(i, j) == CellType::SOLID) continue;
                grid.VOF(i, j) = std::clamp(grid.VOF(i, j) * ratio, 0.0, 1.0);
            }
        return;
    }

    double correction_ratio = std::max((interface_volume - volume_error) / interface_volume, 0.0);
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j) {
            if (grid.Type(i, j) == CellType::SOLID) continue;
            double f = grid.VOF(i, j);
            if (f > 0.0 && f < 1.0)
                grid.VOF(i, j) = std::clamp(f * correction_ratio, 0.0, 1.0);
        }
}

double VOFTransport::total_volume(const Grid& grid) const {
    double vol = 0.0;
    double cell_area = grid.dx * grid.dy;
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            vol += grid.VOF(i, j) * cell_area;
    return vol;
}
