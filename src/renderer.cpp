/**
 * @file renderer.cpp
 * @brief OpenGL rendering implementation for the sloshing tank simulation.
 */

#include "renderer.hpp"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLFW/glfw3.h>
#include <cmath>

bool Renderer::init(int width, int height, const char* title) {
    if (!glfwInit()) {
        return false;
    }

    // Request OpenGL 2.1 compatibility profile (legacy GL for simple 2D)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    window_ = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!window_) {
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // VSync

    win_width_ = width;
    win_height_ = height;

    // Set up orthographic projection matching the window
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0); // Normalized [0,1] x [0,1]
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    return true;
}

void Renderer::shutdown() {
    if (window_) {
        glfwDestroyWindow(window_);
        window_ = nullptr;
    }
    glfwTerminate();
}

void Renderer::vof_to_color(double vof_value, float& r, float& g, float& b) const {
    double f = std::clamp(vof_value, 0.0, 1.0);

    // Smoothstep to sharpen the visual transition without altering physics.
    // Maps the diffuse 0..1 range into a crisper curve centered at 0.5.
    f = f * f * (3.0 - 2.0 * f); // Hermite smoothstep
    f = f * f * (3.0 - 2.0 * f); // Apply twice for a steeper curve

    // Air (f=0):   light sky blue  (0.85, 0.92, 1.0)
    // Water (f=1): deep ocean blue (0.05, 0.15, 0.6)
    r = static_cast<float>(0.85 * (1.0 - f) + 0.05 * f);
    g = static_cast<float>(0.92 * (1.0 - f) + 0.15 * f);
    b = static_cast<float>(1.0  * (1.0 - f) + 0.6  * f);
}

void Renderer::render(const Grid& grid, double sim_time, double fps) {
    glClear(GL_COLOR_BUFFER_BIT);

    int nx = grid.NX;
    int ny = grid.NY;

    glBegin(GL_QUADS);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            float r, g, b;

            if (grid.Type(i, j) == CellType::SOLID) {
                // Walls: dark gray
                r = 0.3f; g = 0.3f; b = 0.3f;
            } else {
                vof_to_color(grid.VOF(i, j), r, g, b);
            }

            glColor3f(r, g, b);

            // Normalized coordinates [0,1]
            float x0 = static_cast<float>(i) / nx;
            float x1 = static_cast<float>(i + 1) / nx;
            float y0 = static_cast<float>(j) / ny;
            float y1 = static_cast<float>(j + 1) / ny;

            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    glfwSwapBuffers(window_);
}

InputState Renderer::poll_input() {
    glfwPollEvents();

    InputState state;
    state.left       = glfwGetKey(window_, GLFW_KEY_LEFT)  == GLFW_PRESS ||
                       glfwGetKey(window_, GLFW_KEY_A)     == GLFW_PRESS;
    state.right      = glfwGetKey(window_, GLFW_KEY_RIGHT) == GLFW_PRESS ||
                       glfwGetKey(window_, GLFW_KEY_D)     == GLFW_PRESS;
    state.up         = glfwGetKey(window_, GLFW_KEY_UP)    == GLFW_PRESS ||
                       glfwGetKey(window_, GLFW_KEY_W)     == GLFW_PRESS;
    state.down       = glfwGetKey(window_, GLFW_KEY_DOWN)  == GLFW_PRESS ||
                       glfwGetKey(window_, GLFW_KEY_S)     == GLFW_PRESS;
    state.reset      = glfwGetKey(window_, GLFW_KEY_R)     == GLFW_PRESS;
    state.pause      = glfwGetKey(window_, GLFW_KEY_SPACE) == GLFW_PRESS;
    state.quit       = glfwGetKey(window_, GLFW_KEY_ESCAPE)== GLFW_PRESS ||
                       glfwGetKey(window_, GLFW_KEY_Q)     == GLFW_PRESS;
    state.speed_up   = glfwGetKey(window_, GLFW_KEY_EQUAL) == GLFW_PRESS; // +/=
    state.speed_down = glfwGetKey(window_, GLFW_KEY_MINUS) == GLFW_PRESS;

    return state;
}

bool Renderer::should_close() const {
    return window_ && glfwWindowShouldClose(window_);
}
