/**
 * @file main.cpp
 * @brief Entry point for the 2D sloshing tank simulation.
 *
 * Parses optional command-line arguments and launches the interactive
 * simulation with OpenGL rendering.
 *
 * Usage: sloshing_tank [nx ny]
 *   nx, ny: Grid resolution (default 256 128)
 */

#include "simulation.hpp"
#include <cstdlib>
#include <cstdio>

int main(int argc, char* argv[]) {
    SimConfig config;

    // Optional: override grid resolution from command line
    if (argc >= 3) {
        config.nx = std::atoi(argv[1]);
        config.ny = std::atoi(argv[2]);
        if (config.nx < 16 || config.ny < 16) {
            std::fprintf(stderr, "Grid resolution must be at least 16x16\n");
            return 1;
        }
        // Scale tank dimensions to maintain aspect ratio
        config.lx = 2.0 * config.nx / 256.0;
        config.ly = 1.0 * config.ny / 128.0;
        // Scale window to match
        config.win_width  = std::min(1920, config.nx * 4);
        config.win_height = std::min(1080, config.ny * 4);
    }

    Simulation sim;
    if (!sim.init(config)) {
        std::fprintf(stderr, "Failed to initialize simulation\n");
        return 1;
    }

    sim.run();
    sim.shutdown();

    return 0;
}
