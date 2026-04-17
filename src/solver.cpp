/**
 * @file solver.cpp
 * @brief Implementation of the fractional-step Navier-Stokes solver.
 *
 * The outer-domain walls are face boundaries (no SOLID border cells). The
 * no-penetration condition u·n = 0 at every wall — outer faces and internal
 * SOLID (baffle) faces — is folded into the pressure Poisson equation as a
 * non-homogeneous Neumann boundary condition:
 *
 *     ∂p/∂n |wall = (ρ/Δt) u*·n
 *
 * Discretely, this is implemented by:
 *   (a) treating the wall-side velocity as zero in the divergence sum used
 *       as the Poisson RHS (the BC contribution exactly cancels it), and
 *   (b) dropping the wall-side neighbor from the pressure stencil
 *       (homogeneous-Neumann ghost-cell elimination).
 * After the pressure solve, projection at a wall face yields u·n = 0 by
 * construction, so no separate BC enforcement step is needed.
 */

#include "solver.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace {
/// True if cell (i, j) acts as a wall for a face: out of bounds or SOLID.
inline bool is_wall_cell(const Grid& g, int i, int j) {
    if (i < 0 || i >= g.NX || j < 0 || j >= g.NY) return true;
    return g.Type(i, j) == CellType::SOLID;
}

/// Effective u at the LEFT face of cell (i, j) for divergence/projection.
/// Treated as 0 if the cell on the other side is a wall.
inline double u_left(const Grid& g, int i, int j) {
    return is_wall_cell(g, i - 1, j) ? 0.0 : g.U(i, j);
}
inline double u_right(const Grid& g, int i, int j) {
    return is_wall_cell(g, i + 1, j) ? 0.0 : g.U(i + 1, j);
}
inline double v_below(const Grid& g, int i, int j) {
    return is_wall_cell(g, i, j - 1) ? 0.0 : g.V(i, j);
}
inline double v_above(const Grid& g, int i, int j) {
    return is_wall_cell(g, i, j + 1) ? 0.0 : g.V(i, j + 1);
}
} // namespace

double Solver::compute_dt(const Grid& grid) const {
    double max_vel = 1e-10;

    for (int i = 0; i <= grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            max_vel = std::max(max_vel, std::abs(grid.U(i, j)));
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j <= grid.NY; ++j)
            max_vel = std::max(max_vel, std::abs(grid.V(i, j)));

    double dt = params.cfl * std::min(grid.dx, grid.dy) / max_vel;
    return std::clamp(dt, params.min_dt, params.max_dt);
}

void Solver::extrapolate_velocity_to_air(Grid& grid, int n_layers) {
    // Normal-direction velocity extrapolation into the air band.
    //
    // For each undefined face, defined neighbours are weighted by how
    // well aligned they are with the direction toward the fluid (the
    // VOF gradient direction). This copies the surface velocity outward
    // along the interface normal rather than blurring it isotropically.
    // Where no gradient information is available, falls back to equal
    // weighting.

    int NX = grid.NX, NY = grid.NY;

    // --- Precompute VOF gradient at cell centres ---
    std::vector<double> gx(NX * NY, 0.0);
    std::vector<double> gy(NX * NY, 0.0);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double fxp = (i+1 < NX) ? grid.VOF(i+1,j) : grid.VOF(i,j);
            double fxm = (i-1 >= 0) ? grid.VOF(i-1,j) : grid.VOF(i,j);
            double fyp = (j+1 < NY) ? grid.VOF(i,j+1) : grid.VOF(i,j);
            double fym = (j-1 >= 0) ? grid.VOF(i,j-1) : grid.VOF(i,j);
            if (i+1 < NX && grid.Type(i+1,j) == CellType::SOLID) fxp = grid.VOF(i,j);
            if (i-1 >= 0 && grid.Type(i-1,j) == CellType::SOLID) fxm = grid.VOF(i,j);
            if (j+1 < NY && grid.Type(i,j+1) == CellType::SOLID) fyp = grid.VOF(i,j);
            if (j-1 >= 0 && grid.Type(i,j-1) == CellType::SOLID) fym = grid.VOF(i,j);
            gx[grid.idx(i,j)] = (fxp - fxm) / (2.0 * grid.dx);
            gy[grid.idx(i,j)] = (fyp - fym) / (2.0 * grid.dy);
        }
    }

    // --- Mark defined faces ---
    std::vector<char> def_u(grid.u.size(), 0);
    std::vector<char> def_v(grid.v.size(), 0);

    for (int i = 0; i <= NX; ++i)
        for (int j = 0; j < NY; ++j) {
            int idx = grid.u_idx(i, j);
            bool fl = (i > 0)  && grid.Type(i-1,j) == CellType::FLUID;
            bool fr = (i < NX) && grid.Type(i,j)   == CellType::FLUID;
            if (fl || fr) { def_u[idx] = 1; continue; }
            bool wl = (i == 0)  || grid.Type(i-1,j) == CellType::SOLID;
            bool wr = (i == NX) || grid.Type(i,j)   == CellType::SOLID;
            if (wl || wr) def_u[idx] = 1;
        }
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j <= NY; ++j) {
            int idx = grid.v_idx(i, j);
            bool fb = (j > 0)  && grid.Type(i,j-1) == CellType::FLUID;
            bool ft = (j < NY) && grid.Type(i,j)   == CellType::FLUID;
            if (fb || ft) { def_v[idx] = 1; continue; }
            bool wb = (j == 0)  || grid.Type(i,j-1) == CellType::SOLID;
            bool wt = (j == NY) || grid.Type(i,j)   == CellType::SOLID;
            if (wb || wt) def_v[idx] = 1;
        }

    // --- Iterative normal-weighted extrapolation ---
    // Neighbour direction vectors on the face grid.
    static constexpr int udirs[][2] = {{-1,0},{1,0},{0,-1},{0,1}};

    for (int pass = 0; pass < n_layers; ++pass) {
        auto def_u_prev = def_u;
        auto def_v_prev = def_v;

        // u-faces
        for (int i = 0; i <= NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                int idx = grid.u_idx(i, j);
                if (def_u_prev[idx]) continue;

                // Face-level gradient: average from neighbouring cells
                double fgx = 0, fgy = 0;
                int gc = 0;
                if (i > 0 && i-1 < NX) { fgx += gx[grid.idx(i-1,j)]; fgy += gy[grid.idx(i-1,j)]; gc++; }
                if (i < NX)            { fgx += gx[grid.idx(i,j)];   fgy += gy[grid.idx(i,j)];   gc++; }
                if (gc > 0) { fgx /= gc; fgy /= gc; }
                double gmag = std::sqrt(fgx*fgx + fgy*fgy);

                double wsum = 0, vsum = 0;
                for (auto& d : udirs) {
                    int ni = i + d[0], nj = j + d[1];
                    if (ni < 0 || ni > NX || nj < 0 || nj >= NY) continue;
                    int nidx = grid.u_idx(ni, nj);
                    if (!def_u_prev[nidx]) continue;

                    double w;
                    if (gmag > 1e-12) {
                        // Weight by alignment with ∇VOF (toward fluid)
                        w = (fgx * d[0] + fgy * d[1]) / gmag;
                        w = std::max(w, 0.0);
                        w = w * w + 0.01; // cos² + small fallback
                    } else {
                        w = 1.0; // no gradient: isotropic
                    }
                    wsum += w;
                    vsum += w * grid.u[nidx];
                }
                if (wsum > 0) { grid.u[idx] = vsum / wsum; def_u[idx] = 1; }
            }
        }

        // v-faces
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j <= NY; ++j) {
                int idx = grid.v_idx(i, j);
                if (def_v_prev[idx]) continue;

                double fgx = 0, fgy = 0;
                int gc = 0;
                if (j > 0 && j-1 < NY) { fgx += gx[grid.idx(i,j-1)]; fgy += gy[grid.idx(i,j-1)]; gc++; }
                if (j < NY)            { fgx += gx[grid.idx(i,j)];   fgy += gy[grid.idx(i,j)];   gc++; }
                if (gc > 0) { fgx /= gc; fgy /= gc; }
                double gmag = std::sqrt(fgx*fgx + fgy*fgy);

                double wsum = 0, vsum = 0;
                for (auto& d : udirs) {
                    int ni = i + d[0], nj = j + d[1];
                    if (ni < 0 || ni >= NX || nj < 0 || nj > NY) continue;
                    int nidx = grid.v_idx(ni, nj);
                    if (!def_v_prev[nidx]) continue;

                    double w;
                    if (gmag > 1e-12) {
                        w = (fgx * d[0] + fgy * d[1]) / gmag;
                        w = std::max(w, 0.0);
                        w = w * w + 0.01;
                    } else {
                        w = 1.0;
                    }
                    wsum += w;
                    vsum += w * grid.v[nidx];
                }
                if (wsum > 0) { grid.v[idx] = vsum / wsum; def_v[idx] = 1; }
            }
        }
    }
}

void Solver::advect_velocity(Grid& grid, double dt) {
    // Semi-Lagrangian back-trace at every interior face. Outer wall faces
    // (i=0, i=NX, j=0, j=NY) are not advected — they stay at zero, which
    // is what the projection produces from the Neumann BC.
    std::vector<double> u_new(grid.u.size(), 0.0);
    std::vector<double> v_new(grid.v.size(), 0.0);

    for (int i = 1; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            double x = i * grid.dx;
            double y = (j + 0.5) * grid.dy;

            auto [vel_u, vel_v] = grid.interpolate_velocity(x, y);

            double x_dep = std::clamp(x - vel_u * dt, 0.0, grid.Lx);
            double y_dep = std::clamp(y - vel_v * dt, 0.0, grid.Ly);

            auto [u_dep, _] = grid.interpolate_velocity(x_dep, y_dep);
            u_new[grid.u_idx(i, j)] = u_dep;
        }
    }

    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 1; j < grid.NY; ++j) {
            double x = (i + 0.5) * grid.dx;
            double y = j * grid.dy;

            auto [vel_u, vel_v] = grid.interpolate_velocity(x, y);

            double x_dep = std::clamp(x - vel_u * dt, 0.0, grid.Lx);
            double y_dep = std::clamp(y - vel_v * dt, 0.0, grid.Ly);

            auto [_, v_dep] = grid.interpolate_velocity(x_dep, y_dep);
            v_new[grid.v_idx(i, j)] = v_dep;
        }
    }

    grid.u = std::move(u_new);
    grid.v = std::move(v_new);
}

void Solver::apply_body_forces(Grid& grid, double dt, double tank_ax, double tank_ay) {
    // Effective acceleration in the tank frame: g - a_tank.
    double ax = -tank_ax;
    double ay = -params.gravity - tank_ay;

    // u-faces: apply at every interior face touching at least one FLUID cell.
    // (Wall-adjacent contributions get absorbed by the non-homogeneous Neumann
    //  BC in the pressure solve, so no skipping at SOLID-adjacent faces.)
    for (int i = 1; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            CellType tl = grid.Type(i - 1, j);
            CellType tr = grid.Type(i, j);
            if (tl == CellType::FLUID || tr == CellType::FLUID) {
                grid.U(i, j) += ax * dt;
            }
        }
    }

    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 1; j < grid.NY; ++j) {
            CellType tb = grid.Type(i, j - 1);
            CellType tt = grid.Type(i, j);
            if (tb == CellType::FLUID || tt == CellType::FLUID) {
                grid.V(i, j) += ay * dt;
            }
        }
    }
}

int Solver::solve_pressure(Grid& grid, double dt) {
    // Warm-start: keep previous pressure as initial guess; zero in non-fluid.
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            if (grid.Type(i, j) != CellType::FLUID)
                grid.P(i, j) = 0.0;

    double scale = params.density / dt;
    double dx2 = grid.dx * grid.dx;
    double dy2 = grid.dy * grid.dy;

    auto sor_sweep = [&](int color, double& max_residual) {
        for (int i = 0; i < grid.NX; ++i) {
            for (int j = 0; j < grid.NY; ++j) {
                if ((i + j) % 2 != color) continue;
                if (grid.Type(i, j) != CellType::FLUID) continue;

                // BC-absorbed divergence: wall-adjacent face velocities are 0.
                // The (ρ/Δt) u*·n contribution from the non-homogeneous
                // Neumann BC exactly cancels the actual stored value.
                double div = (u_right(grid, i, j) - u_left(grid, i, j)) / grid.dx
                           + (v_above(grid, i, j) - v_below(grid, i, j)) / grid.dy;

                double p_sum = 0.0;
                double coeff = 0.0;

                auto add_neighbor = [&](int ni, int nj, double h2) {
                    // Out-of-bounds (outer wall) or SOLID neighbor: drop from
                    // the stencil (Neumann ghost-cell elimination).
                    if (ni < 0 || ni >= grid.NX || nj < 0 || nj >= grid.NY) return;
                    CellType t = grid.Type(ni, nj);
                    if (t == CellType::FLUID) {
                        p_sum += grid.P(ni, nj) / h2;
                        coeff += 1.0 / h2;
                    } else if (t == CellType::EMPTY) {
                        coeff += 1.0 / h2; // Dirichlet p=0 at free surface
                    }
                };

                add_neighbor(i - 1, j, dx2);
                add_neighbor(i + 1, j, dx2);
                add_neighbor(i, j - 1, dy2);
                add_neighbor(i, j + 1, dy2);

                if (coeff < 1e-10) continue;

                double p_new = (p_sum - scale * div) / coeff;
                double residual = std::abs(p_new - grid.P(i, j));
                max_residual = std::max(max_residual, residual);

                grid.P(i, j) += params.sor_omega * (p_new - grid.P(i, j));
            }
        }
    };

    int iter = 0;
    for (; iter < params.pressure_iters; ++iter) {
        double max_residual = 0.0;
        sor_sweep(0, max_residual);
        sor_sweep(1, max_residual);
        if (max_residual < params.pressure_tol) break;
    }

    return iter + 1;
}

void Solver::project_velocity(Grid& grid, double dt) {
    double inv_rho_dt = dt / params.density;

    // u-faces. For interior faces between non-wall cells, project with the
    // stored pressure gradient. For wall-adjacent faces (outer boundary or
    // FLUID/SOLID interface), the Neumann BC outcome is u·n = 0 — assign
    // directly. There's no enforce_boundary_conditions safety net.
    for (int i = 0; i <= grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            bool wall_l = is_wall_cell(grid, i - 1, j);
            bool wall_r = is_wall_cell(grid, i, j);
            if (wall_l || wall_r) {
                grid.U(i, j) = 0.0;
                continue;
            }
            CellType tl = grid.Type(i - 1, j);
            CellType tr = grid.Type(i, j);
            if (tl == CellType::FLUID || tr == CellType::FLUID) {
                grid.U(i, j) -= inv_rho_dt * (grid.P(i, j) - grid.P(i - 1, j)) / grid.dx;
            }
        }
    }

    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j <= grid.NY; ++j) {
            bool wall_b = is_wall_cell(grid, i, j - 1);
            bool wall_t = is_wall_cell(grid, i, j);
            if (wall_b || wall_t) {
                grid.V(i, j) = 0.0;
                continue;
            }
            CellType tb = grid.Type(i, j - 1);
            CellType tt = grid.Type(i, j);
            if (tb == CellType::FLUID || tt == CellType::FLUID) {
                grid.V(i, j) -= inv_rho_dt * (grid.P(i, j) - grid.P(i, j - 1)) / grid.dy;
            }
        }
    }
}

double Solver::max_divergence(const Grid& grid) const {
    double max_div = 0.0;

    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            if (grid.Type(i, j) != CellType::FLUID) continue;

            double div = (u_right(grid, i, j) - u_left(grid, i, j)) / grid.dx
                       + (v_above(grid, i, j) - v_below(grid, i, j)) / grid.dy;
            max_div = std::max(max_div, std::abs(div));
        }
    }

    return max_div;
}

void Solver::apply_surface_tension(Grid& grid, double dt) {
    if (params.surface_tension <= 0.0) return;

    int NX = grid.NX;
    int NY = grid.NY;
    double dx = grid.dx;
    double dy = grid.dy;

    // Curvature κ = -∇·n̂ where n̂ = ∇α / |∇α|.
    // When reading VOF from a SOLID (baffle) neighbour, mirror the center
    // cell's own value so the artificial VOF=0 inside the obstacle is not
    // mistaken for a water-air interface.
    std::vector<double> kappa(NX * NY, 0.0);

    auto fv = [&](int ci, int cj, int ni, int nj) -> double {
        if (ni < 0 || ni >= NX || nj < 0 || nj >= NY) return grid.VOF(ci, cj);
        if (grid.Type(ni, nj) == CellType::SOLID)       return grid.VOF(ci, cj);
        return grid.VOF(ni, nj);
    };

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            double f_ip = fv(i, j, i+1, j), f_im = fv(i, j, i-1, j);
            double f_jp = fv(i, j, i, j+1), f_jm = fv(i, j, i, j-1);

            double dfdx = (f_ip - f_im) / (2.0 * dx);
            double dfdy = (f_jp - f_jm) / (2.0 * dy);
            double mag = std::sqrt(dfdx * dfdx + dfdy * dfdy);
            if (mag < 1e-8) continue;

            double f_c = grid.VOF(i, j);
            double d2fdx2 = (f_ip - 2.0*f_c + f_im) / (dx * dx);
            double d2fdy2 = (f_jp - 2.0*f_c + f_jm) / (dy * dy);

            double f_pp = fv(i, j, i+1, j+1), f_mp = fv(i, j, i-1, j+1);
            double f_pm = fv(i, j, i+1, j-1), f_mm = fv(i, j, i-1, j-1);
            double d2fdxdy = (f_pp - f_mp - f_pm + f_mm) / (4.0 * dx * dy);

            double mag3 = mag * mag * mag;
            kappa[grid.idx(i, j)] = -(dfdx * dfdx * d2fdy2
                                     - 2.0 * dfdx * dfdy * d2fdxdy
                                     + dfdy * dfdy * d2fdx2) / mag3;
        }
    }

    double sigma_over_rho = params.surface_tension / params.density;

    for (int i = 1; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            if (is_wall_cell(grid, i - 1, j) || is_wall_cell(grid, i, j)) continue;
            double kappa_face = 0.5 * (kappa[grid.idx(i - 1, j)] + kappa[grid.idx(i, j)]);
            double dfdx = (grid.VOF(i, j) - grid.VOF(i - 1, j)) / dx;
            grid.U(i, j) += dt * sigma_over_rho * kappa_face * dfdx;
        }
    }

    for (int i = 0; i < NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            if (is_wall_cell(grid, i, j - 1) || is_wall_cell(grid, i, j)) continue;
            double kappa_face = 0.5 * (kappa[grid.idx(i, j - 1)] + kappa[grid.idx(i, j)]);
            double dfdy = (grid.VOF(i, j) - grid.VOF(i, j - 1)) / dy;
            grid.V(i, j) += dt * sigma_over_rho * kappa_face * dfdy;
        }
    }
}

double Solver::step(Grid& grid, double tank_ax, double tank_ay) {
    double dt = compute_dt(grid);

    grid.classify_cells();
    extrapolate_velocity_to_air(grid, 1);
    advect_velocity(grid, dt);
    apply_body_forces(grid, dt, tank_ax, tank_ay);
    apply_surface_tension(grid, dt);
    solve_pressure(grid, dt);
    project_velocity(grid, dt);

    last_max_divergence = max_divergence(grid);

    return dt;
}
