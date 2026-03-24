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

    double dt_cfl = params.cfl * std::min(grid.dx, grid.dy) / max_vel;

    // Viscous stability limit (explicit diffusion)
    double dt_visc = 0.25 * std::min(grid.dx * grid.dx, grid.dy * grid.dy)
                     * params.density / std::max(params.viscosity, 1e-10);

    double dt = std::min(dt_cfl, dt_visc);
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

    int iter = 0;
    for (; iter < params.pressure_iters; ++iter) {
        double max_residual = 0.0;

        for (int i = 1; i < grid.NX - 1; ++i) {
            for (int j = 1; j < grid.NY - 1; ++j) {
                if (grid.Type(i, j) != CellType::FLUID) continue;

                // Divergence of intermediate velocity at cell (i,j)
                double div = (grid.U(i + 1, j) - grid.U(i, j)) / grid.dx
                           + (grid.V(i, j + 1) - grid.V(i, j)) / grid.dy;

                // Count valid neighbors and sum their pressures
                // SOLID neighbor: ∂p/∂n = 0 → use p_center
                // EMPTY neighbor: p = 0 (free surface)
                // FLUID neighbor: use p_neighbor
                double p_sum = 0.0;
                int n_neighbors = 0;
                double coeff = 0.0;

                auto add_neighbor = [&](int ni, int nj, double h2) {
                    CellType t = grid.Type(ni, nj);
                    if (t == CellType::FLUID) {
                        p_sum += grid.P(ni, nj) / h2;
                        coeff += 1.0 / h2;
                    } else if (t == CellType::EMPTY) {
                        // p = 0 at free surface → contributes 0 to sum
                        coeff += 1.0 / h2;
                    }
                    // SOLID: Neumann BC → skip (effectively p_neighbor = p_center)
                };

                add_neighbor(i - 1, j, dx2);
                add_neighbor(i + 1, j, dx2);
                add_neighbor(i, j - 1, dy2);
                add_neighbor(i, j + 1, dy2);

                if (coeff < 1e-10) continue; // Isolated cell

                double p_new = (p_sum - scale * div) / coeff;
                double residual = std::abs(p_new - grid.P(i, j));
                max_residual = std::max(max_residual, residual);

                // SOR update
                grid.P(i, j) += params.sor_omega * (p_new - grid.P(i, j));
            }
        }

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

double Solver::step(Grid& grid, double tank_ax, double tank_ay) {
    double dt = compute_dt(grid);

    // 1. Classify cells based on current VOF
    grid.classify_cells();

    // 2. Advect velocity (semi-Lagrangian)
    advect_velocity(grid, dt);

    // 3. Apply body forces (gravity + tank acceleration)
    apply_body_forces(grid, dt, tank_ax, tank_ay);

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
