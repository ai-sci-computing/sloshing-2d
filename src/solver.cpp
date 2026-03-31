/**
 * @file solver.cpp
 * @brief Implementation of the fractional-step Navier-Stokes solver.
 */

#include "solver.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

double Solver::compute_dt(const Grid& grid) const {
    double max_vel = 1e-10; // Avoid division by zero

    for (int i = 1; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            max_vel = std::max(max_vel, std::abs(grid.U(i, j)));
        }
    }
    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 1; j < grid.NY; ++j) {
            max_vel = std::max(max_vel, std::abs(grid.V(i, j)));
        }
    }

    double dt = params.cfl * std::min(grid.dx, grid.dy) / max_vel;
    return std::clamp(dt, params.min_dt, params.max_dt);
}

void Solver::advect_velocity(Grid& grid, double dt) {
    // Semi-Lagrangian advection: trace back from each velocity sample point

    std::vector<double> u_new(grid.u.size(), 0.0);
    std::vector<double> v_new(grid.v.size(), 0.0);

    // Advect u-velocity (stored at vertical faces: position (i*dx, (j+0.5)*dy))
    for (int i = 1; i < grid.NX; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            double x = i * grid.dx;
            double y = (j + 0.5) * grid.dy;

            // Velocity at this u-sample point (need to interpolate v here too)
            auto [vel_u, vel_v] = grid.interpolate_velocity(x, y);

            // Trace backward
            double x_dep = x - vel_u * dt;
            double y_dep = y - vel_v * dt;

            // Clamp to domain
            x_dep = std::clamp(x_dep, 0.0, grid.Lx);
            y_dep = std::clamp(y_dep, 0.0, grid.Ly);

            // Interpolate u at departure point
            auto [u_dep, v_dep_unused] = grid.interpolate_velocity(x_dep, y_dep);
            u_new[grid.u_idx(i, j)] = u_dep;
        }
    }

    // Advect v-velocity (stored at horizontal faces: position ((i+0.5)*dx, j*dy))
    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY; ++j) {
            double x = (i + 0.5) * grid.dx;
            double y = j * grid.dy;

            auto [vel_u, vel_v] = grid.interpolate_velocity(x, y);

            double x_dep = x - vel_u * dt;
            double y_dep = y - vel_v * dt;

            x_dep = std::clamp(x_dep, 0.0, grid.Lx);
            y_dep = std::clamp(y_dep, 0.0, grid.Ly);

            auto [u_dep_unused, v_dep] = grid.interpolate_velocity(x_dep, y_dep);
            v_new[grid.v_idx(i, j)] = v_dep;
        }
    }

    grid.u = std::move(u_new);
    grid.v = std::move(v_new);
}

void Solver::apply_body_forces(Grid& grid, double dt, double tank_ax, double tank_ay) {
    // In the tank's reference frame, effective acceleration is (g - a_tank).
    // Gravity acts downward: g_x = 0, g_y = -gravity
    double ax = 0.0 - tank_ax;
    double ay = -params.gravity - tank_ay;

    // Apply body forces to velocity faces adjacent to fluid cells.
    // Full-strength gravity is applied at any face where at least one
    // neighbor is a FLUID cell. The pressure solver (with warm-start and
    // free-surface BC p=0) handles the balance at the surface.
    for (int i = 2; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            if (grid.Type(i - 1, j) == CellType::FLUID ||
                grid.Type(i, j) == CellType::FLUID) {
                grid.U(i, j) += ax * dt;
            }
        }
    }

    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 2; j < grid.NY - 1; ++j) {
            if (grid.Type(i, j - 1) == CellType::FLUID ||
                grid.Type(i, j) == CellType::FLUID) {
                grid.V(i, j) += ay * dt;
            }
        }
    }
}

int Solver::solve_pressure(Grid& grid, double dt) {
    // Warm-start: keep previous pressure field as initial guess.
    // Only zero out pressure in cells that are no longer fluid.
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            if (grid.Type(i, j) != CellType::FLUID)
                grid.P(i, j) = 0.0;

    double scale = params.density / dt;
    double dx2 = grid.dx * grid.dx;
    double dy2 = grid.dy * grid.dy;

    // Red-Black SOR: sweep red cells (i+j even) then black cells (i+j odd).
    // Within each color, updates are independent and can be parallelized.
    auto sor_sweep = [&](int color, double& max_residual) {
        for (int i = 1; i < grid.NX - 1; ++i) {
            for (int j = 1; j < grid.NY - 1; ++j) {
                if ((i + j) % 2 != color) continue;
                if (grid.Type(i, j) != CellType::FLUID) continue;

                double div = (grid.U(i + 1, j) - grid.U(i, j)) / grid.dx
                           + (grid.V(i, j + 1) - grid.V(i, j)) / grid.dy;

                double p_sum = 0.0;
                double coeff = 0.0;

                auto add_neighbor = [&](int ni, int nj, double h2) {
                    CellType t = grid.Type(ni, nj);
                    if (t == CellType::FLUID) {
                        p_sum += grid.P(ni, nj) / h2;
                        coeff += 1.0 / h2;
                    } else if (t == CellType::EMPTY) {
                        coeff += 1.0 / h2;
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
        sor_sweep(0, max_residual); // Red
        sor_sweep(1, max_residual); // Black
        if (max_residual < params.pressure_tol) break;
    }

    return iter + 1;
}

void Solver::project_velocity(Grid& grid, double dt) {
    double inv_rho_dt = dt / params.density;

    // Update u-velocity: u = u* - (dt/rho) * dp/dx
    for (int i = 2; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            // Only update if at least one side is fluid
            if (grid.Type(i - 1, j) == CellType::FLUID ||
                grid.Type(i, j) == CellType::FLUID) {
                grid.U(i, j) -= inv_rho_dt * (grid.P(i, j) - grid.P(i - 1, j)) / grid.dx;
            }
        }
    }

    // Update v-velocity: v = v* - (dt/rho) * dp/dy
    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 2; j < grid.NY - 1; ++j) {
            if (grid.Type(i, j - 1) == CellType::FLUID ||
                grid.Type(i, j) == CellType::FLUID) {
                grid.V(i, j) -= inv_rho_dt * (grid.P(i, j) - grid.P(i, j - 1)) / grid.dy;
            }
        }
    }
}

double Solver::max_divergence(const Grid& grid) const {
    double max_div = 0.0;

    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY - 1; ++j) {
            if (grid.Type(i, j) != CellType::FLUID) continue;

            double div = (grid.U(i + 1, j) - grid.U(i, j)) / grid.dx
                       + (grid.V(i, j + 1) - grid.V(i, j)) / grid.dy;
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

    // Compute smoothed VOF gradient components and curvature at cell centers.
    // Curvature: κ = -∇·n̂, where n̂ = ∇α / |∇α|
    std::vector<double> kappa(NX * NY, 0.0);

    // Start at i=2,j=2 to avoid stencil reaching into SOLID boundary cells
    // where VOF=0 creates artificial gradients and spurious curvature.
    for (int i = 2; i < NX - 2; ++i) {
        for (int j = 2; j < NY - 2; ++j) {
            double dfdx = (grid.VOF(i + 1, j) - grid.VOF(i - 1, j)) / (2.0 * dx);
            double dfdy = (grid.VOF(i, j + 1) - grid.VOF(i, j - 1)) / (2.0 * dy);
            double mag = std::sqrt(dfdx * dfdx + dfdy * dfdy);

            if (mag < 1e-8) continue; // Not at interface

            // Second derivatives for curvature via divergence of the normal
            double d2fdx2 = (grid.VOF(i + 1, j) - 2.0 * grid.VOF(i, j) + grid.VOF(i - 1, j)) / (dx * dx);
            double d2fdy2 = (grid.VOF(i, j + 1) - 2.0 * grid.VOF(i, j) + grid.VOF(i, j - 1)) / (dy * dy);

            // Mixed derivative
            double d2fdxdy = (grid.VOF(i + 1, j + 1) - grid.VOF(i - 1, j + 1)
                            - grid.VOF(i + 1, j - 1) + grid.VOF(i - 1, j - 1))
                           / (4.0 * dx * dy);

            // κ = -(dfdx²·d2fdy2 - 2·dfdx·dfdy·d2fdxdy + dfdy²·d2fdx2) / |∇f|³
            double mag3 = mag * mag * mag;
            kappa[grid.idx(i, j)] = -(dfdx * dfdx * d2fdy2
                                     - 2.0 * dfdx * dfdy * d2fdxdy
                                     + dfdy * dfdy * d2fdx2) / mag3;
        }
    }

    // Apply CSF force to u-velocity faces: F_x = σ·κ·(∂α/∂x) / ρ
    double sigma_over_rho = params.surface_tension / params.density;

    for (int i = 2; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            // Average curvature to face
            double kappa_face = 0.5 * (kappa[grid.idx(i - 1, j)] + kappa[grid.idx(i, j)]);
            // VOF gradient at face
            double dfdx = (grid.VOF(i, j) - grid.VOF(i - 1, j)) / dx;

            grid.U(i, j) += dt * sigma_over_rho * kappa_face * dfdx;
        }
    }

    // Apply CSF force to v-velocity faces: F_y = σ·κ·(∂α/∂y) / ρ
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 2; j < NY - 1; ++j) {
            double kappa_face = 0.5 * (kappa[grid.idx(i, j - 1)] + kappa[grid.idx(i, j)]);
            double dfdy = (grid.VOF(i, j) - grid.VOF(i, j - 1)) / dy;

            grid.V(i, j) += dt * sigma_over_rho * kappa_face * dfdy;
        }
    }
}

double Solver::step(Grid& grid, double tank_ax, double tank_ay) {
    double dt = compute_dt(grid);

    // 1. Classify cells based on current VOF
    grid.classify_cells();

    // 2. Advect velocity (semi-Lagrangian)
    advect_velocity(grid, dt);

    // 3. Apply body forces (gravity + tank acceleration)
    apply_body_forces(grid, dt, tank_ax, tank_ay);

    // 3b. Apply surface tension (CSF model)
    apply_surface_tension(grid, dt);

    // 4. Enforce BCs before pressure solve
    grid.enforce_boundary_conditions();

    // 5. Solve pressure Poisson equation
    solve_pressure(grid, dt);

    // 6. Project velocity to be divergence-free
    project_velocity(grid, dt);

    // 7. Enforce BCs after projection
    grid.enforce_boundary_conditions();

    return dt;
}
