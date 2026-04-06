/**
 * @file renderer.hpp
 * @brief OpenGL renderer for visualizing the 2D sloshing tank simulation.
 *
 * Uses GLFW for window management and legacy OpenGL for rendering the VOF
 * field as colored quads. Each grid cell is drawn with a color interpolated
 * between air (white/light blue) and water (deep blue) based on its VOF value.
 */

#pragma once

#include "grid.hpp"
#include <vector>
#include <cstdint>

// Forward-declare GLFWwindow to avoid including GLFW in the header
struct GLFWwindow;

/**
 * @struct InputState
 * @brief Captures keyboard input for tank motion control.
 */
struct InputState {
    bool left  = false; ///< Tank accelerating left
    bool right = false; ///< Tank accelerating right
    bool up    = false; ///< Tank accelerating up
    bool down  = false; ///< Tank accelerating down
    bool reset = false; ///< Reset simulation
    bool pause = false; ///< Toggle pause
    bool quit  = false; ///< Exit application
    bool speed_up   = false; ///< Increase simulation speed
    bool speed_down = false; ///< Decrease simulation speed
};

/**
 * @class Renderer
 * @brief Manages the GLFW window and OpenGL rendering of the simulation state.
 */
class Renderer {
public:
    /**
     * @brief Initialize GLFW and create the rendering window.
     * @param width Window width in pixels
     * @param height Window height in pixels
     * @param title Window title
     * @return true on success, false on failure
     */
    bool init(int width, int height, const char* title);

    /**
     * @brief Shut down GLFW and release resources.
     */
    void shutdown();

    /**
     * @brief Render the current simulation state.
     *
     * Draws the VOF field as colored quads and optionally overlays
     * additional information (velocity arrows, simulation stats).
     *
     * @param grid The simulation grid to render
     * @param sim_time Current simulation time [s]
     * @param fps Current frames per second
     */
    void render(const Grid& grid, double sim_time, double fps);

    /**
     * @brief Poll input events and return the current input state.
     * @return InputState with currently pressed keys
     */
    InputState poll_input();

    /**
     * @brief Check if the window should close.
     * @return true if the user requested window close
     */
    bool should_close() const;

private:
    GLFWwindow* window_ = nullptr;
    int win_width_  = 0;
    int win_height_ = 0;
    unsigned int texture_id_ = 0;
    std::vector<uint8_t> pixel_buf_; ///< RGB pixel buffer for texture upload

    /**
     * @brief Map a VOF value to an RGBA color.
     *
     * VOF=0 → light sky blue (air), VOF=1 → deep ocean blue (water).
     * Partial fill creates a smooth gradient.
     *
     * @param vof_value Volume fraction [0,1]
     * @param r Output red component [0,1]
     * @param g Output green component [0,1]
     * @param b Output blue component [0,1]
     */
    void vof_to_color(double vof_value, float& r, float& g, float& b) const;
};
