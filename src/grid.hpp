/**
 * @file grid.hpp
 * @brief MAC (Marker-And-Cell) staggered grid for 2D incompressible flow.
 *
 * The grid uses a staggered arrangement where:
 * - Pressure and VOF are stored at cell centers: p[i][j], vof[i][j]
 * - Horizontal velocity is stored at vertical cell faces: u[i][j] (NX+1 x NY)
 * - Vertical velocity is stored at horizontal cell faces: v[i][j] (NX x NY+1)
 *
 * This staggering prevents checkerboard pressure oscillations and provides
 * natural locations for applying the divergence and gradient operators.
 */

#pragma once

#include <vector>
#include <cstddef>

/// @brief Cell classification for boundary condition handling.
enum class CellType : uint8_t {
    EMPTY, ///< Air cell (no fluid)
    FLUID, ///< Fully or partially filled with fluid
    SOLID  ///< Solid wall boundary
};

/**
 * @class Grid
 * @brief 2D MAC staggered grid holding all simulation field data.
 *
 * Coordinate system: i is the horizontal index (0..NX-1), j is the vertical
 * index (0..NY-1). j=0 is the bottom of the tank.
 */
class Grid {
public:
    const int NX; ///< Number of cells in horizontal direction
    const int NY; ///< Number of cells in vertical direction
    const double dx; ///< Cell width  [m]
    const double dy; ///< Cell height [m]
    const double Lx; ///< Tank physical width  [m]
    const double Ly; ///< Tank physical height [m]

    // --- Field arrays (row-major: index = i * NY + j for cell-centered) ---

    std::vector<double> u;   ///< Horizontal velocity at vertical faces, size (NX+1)*NY
    std::vector<double> v;   ///< Vertical velocity at horizontal faces, size NX*(NY+1)
    std::vector<double> p;   ///< Pressure at cell centers, size NX*NY
    std::vector<double> vof; ///< Volume fraction at cell centers, size NX*NY [0=empty, 1=full]
    std::vector<CellType> cell_type; ///< Cell classification, size NX*NY

    /**
     * @brief Construct a grid with given resolution and physical dimensions.
     * @param nx Number of horizontal cells
     * @param ny Number of vertical cells
     * @param lx Physical width [m]
     * @param ly Physical height [m]
     */
    Grid(int nx, int ny, double lx, double ly);

    // --- Indexing helpers ---

    /// @brief Cell-centered field index (pressure, VOF, cell_type).
    inline int idx(int i, int j) const { return i * NY + j; }

    /// @brief Horizontal velocity index. u is (NX+1) wide, NY tall.
    inline int u_idx(int i, int j) const { return i * NY + j; }

    /// @brief Vertical velocity index. v is NX wide, (NY+1) tall.
    inline int v_idx(int i, int j) const { return i * (NY + 1) + j; }

    // --- Field accessors with bounds awareness ---

    /// @brief Access u-velocity at face (i, j). i ∈ [0, NX], j ∈ [0, NY-1].
    inline double& U(int i, int j) { return u[u_idx(i, j)]; }
    inline double  U(int i, int j) const { return u[u_idx(i, j)]; }

    /// @brief Access v-velocity at face (i, j). i ∈ [0, NX-1], j ∈ [0, NY].
    inline double& V(int i, int j) { return v[v_idx(i, j)]; }
    inline double  V(int i, int j) const { return v[v_idx(i, j)]; }

    /// @brief Access pressure at cell center (i, j).
    inline double& P(int i, int j) { return p[idx(i, j)]; }
    inline double  P(int i, int j) const { return p[idx(i, j)]; }

    /// @brief Access VOF at cell center (i, j).
    inline double& VOF(int i, int j) { return vof[idx(i, j)]; }
    inline double  VOF(int i, int j) const { return vof[idx(i, j)]; }

    /// @brief Access cell type at (i, j).
    inline CellType& Type(int i, int j) { return cell_type[idx(i, j)]; }
    inline CellType  Type(int i, int j) const { return cell_type[idx(i, j)]; }

    /**
     * @brief Classify cells as FLUID, EMPTY, or SOLID based on VOF values.
     *
     * Wall cells (i=0, i=NX-1, j=0, j=NY-1) are always SOLID.
     * Interior cells with VOF > threshold are FLUID, otherwise EMPTY.
     *
     * @param threshold VOF threshold for fluid classification (default 0.01)
     */
    void classify_cells(double threshold = 0.01);

    /**
     * @brief Enforce no-slip / no-penetration boundary conditions at walls.
     *
     * Sets normal velocity to zero at wall faces and reflects tangential
     * velocity for free-slip (or negates for no-slip).
     */
    void enforce_boundary_conditions();

    /**
     * @brief Initialize a still-water configuration.
     * @param fill_fraction Fraction of tank height filled with water (0 to 1)
     */
    void initialize_water(double fill_fraction);

    /**
     * @brief Compute the velocity at an arbitrary point using bilinear interpolation.
     * @param x Physical x-coordinate
     * @param y Physical y-coordinate
     * @return Pair of (u, v) interpolated velocities
     */
    std::pair<double, double> interpolate_velocity(double x, double y) const;
};
