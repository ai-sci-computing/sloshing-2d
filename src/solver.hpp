/**
 * @file solver.hpp
 * @brief Incompressible Navier-Stokes solver using the fractional-step method.
 *
 * Implements a projection method on the MAC staggered grid:
 * 1. Semi-Lagrangian advection of velocity
 * 2. Application of body forces (gravity + tank acceleration)
 * 3. Pressure Poisson solve via SOR (only in fluid/surface cells)
 * 4. Velocity projection to enforce divergence-free constraint
 *
 * The solver operates only on the interior of the grid; boundary conditions
 * are enforced separately by Grid::enforce_boundary_conditions().
 */

#pragma once

#include "grid.hpp"

/**
 * @struct SolverParams
 * @brief Physical and numerical parameters for the solver.
 */
struct SolverParams {
    double gravity = 9.81;       ///< Gravitational acceleration [m/s²]
    double density = 1000.0;     ///< Fluid density [kg/m³]
    double viscosity = 0.001;    ///< Dynamic viscosity [Pa·s] (may be increased for stability)
    double cfl = 0.5;            ///< CFL number for adaptive timestep
    double min_dt = 1e-6;        ///< Minimum timestep [s]
    double max_dt = 0.01;        ///< Maximum timestep [s]
    int pressure_iters = 100;    ///< Maximum SOR iterations for pressure solve
    double pressure_tol = 1e-4;  ///< Pressure solver convergence tolerance
    double sor_omega = 1.7;      ///< SOR relaxation factor (1 < ω < 2)
};

/**
 * @class Solver
 * @brief Fractional-step Navier-Stokes solver for incompressible free-surface flow.
 */
class Solver {
public:
    SolverParams params;

    Solver() = default;
    explicit Solver(const SolverParams& p) : params(p) {}

    /**
     * @brief Perform one complete simulation timestep.
     *
     * Executes: advection → body forces → pressure solve → projection → BCs.
     *
     * @param grid The simulation grid (modified in place)
     * @param tank_ax Tank acceleration in x [m/s²] (in world frame)
     * @param tank_ay Tank acceleration in y [m/s²] (in world frame)
     * @return The timestep Δt that was used [s]
     */
    double step(Grid& grid, double tank_ax, double tank_ay);

    /**
     * @brief Compute CFL-limited timestep from current velocity field.
     * @param grid The simulation grid
     * @return Adaptive timestep [s]
     */
    double compute_dt(const Grid& grid) const;

    /**
     * @brief Semi-Lagrangian advection of the velocity field.
     *
     * Traces particles backward in time from each velocity sample point and
     * interpolates the velocity at the departure point. Unconditionally stable
     * regardless of CFL number, though accuracy degrades for large timesteps.
     *
     * @param grid The simulation grid
     * @param dt Timestep [s]
     */
    void advect_velocity(Grid& grid, double dt);

    /**
     * @brief Apply body forces (gravity and tank pseudo-forces) to velocity.
     *
     * In the tank's non-inertial reference frame, the effective body force is:
     * F = ρ(g - a_tank), so the acceleration added to each velocity component is
     * (g - a_tank) * dt.
     *
     * @param grid The simulation grid
     * @param dt Timestep [s]
     * @param tank_ax Tank acceleration in x [m/s²]
     * @param tank_ay Tank acceleration in y [m/s²]
     */
    void apply_body_forces(Grid& grid, double dt, double tank_ax, double tank_ay);

    /**
     * @brief Solve the pressure Poisson equation using SOR.
     *
     * Solves ∇²p = (ρ/Δt)∇·u* in fluid cells, with:
     * - p = 0 at free-surface (EMPTY neighbor) boundaries
     * - ∂p/∂n = 0 at SOLID wall boundaries
     *
     * @param grid The simulation grid
     * @param dt Timestep [s]
     * @return Number of iterations used
     */
    int solve_pressure(Grid& grid, double dt);

    /**
     * @brief Project velocity to enforce the divergence-free constraint.
     *
     * Subtracts the pressure gradient from the intermediate velocity:
     * u = u* - (Δt/ρ) ∇p
     *
     * @param grid The simulation grid
     * @param dt Timestep [s]
     */
    void project_velocity(Grid& grid, double dt);

    /**
     * @brief Compute the maximum divergence of the velocity field.
     *
     * Used for verifying that the projection step succeeded.
     *
     * @param grid The simulation grid
     * @return Maximum absolute divergence across all fluid cells
     */
    double max_divergence(const Grid& grid) const;
};
