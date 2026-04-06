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
#include <cstdio>

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
    glfwSwapInterval(0); // No VSync — let simulation run as fast as possible

    win_width_ = width;
    win_height_ = height;

    // Store this pointer for the resize callback
    glfwSetWindowUserPointer(window_, this);

    // Framebuffer resize callback — updates viewport and projection
    glfwSetFramebufferSizeCallback(window_, [](GLFWwindow* w, int width, int height) {
        auto* self = static_cast<Renderer*>(glfwGetWindowUserPointer(w));
        self->win_width_ = width;
        self->win_height_ = height;
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
    });

    // Set up initial orthographic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Create texture for VOF field rendering
    glGenTextures(1, &texture_id_);
    glBindTexture(GL_TEXTURE_2D, texture_id_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glEnable(GL_TEXTURE_2D);

    return true;
}

void Renderer::shutdown() {
    if (texture_id_) {
        glDeleteTextures(1, &texture_id_);
        texture_id_ = 0;
    }
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

    // Build RGB pixel buffer (column-major grid → row-major texture)
    // Texture pixel (x, y) maps to grid cell (i=x, j=y) with y=0 at bottom.
    pixel_buf_.resize(nx * ny * 3);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float r, g, b;
            if (grid.Type(i, j) == CellType::SOLID) {
                r = 0.3f; g = 0.3f; b = 0.3f;
            } else {
                vof_to_color(grid.VOF(i, j), r, g, b);
            }
            int offset = (j * nx + i) * 3;
            pixel_buf_[offset + 0] = static_cast<uint8_t>(r * 255.0f);
            pixel_buf_[offset + 1] = static_cast<uint8_t>(g * 255.0f);
            pixel_buf_[offset + 2] = static_cast<uint8_t>(b * 255.0f);
        }
    }

    // Upload to texture
    glBindTexture(GL_TEXTURE_2D, texture_id_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, pixel_buf_.data());

    // Draw a single textured quad covering [0,1] x [0,1]
    glColor3f(1.0f, 1.0f, 1.0f); // Texture color is not modulated
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(0.0f, 0.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0f, 0.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0f, 1.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(0.0f, 1.0f);
    glEnd();

    glfwSwapBuffers(window_);

    // Update window title with FPS and simulation time
    char title[128];
    std::snprintf(title, sizeof(title), "Sloshing Tank  |  %.1f FPS  |  t = %.2f s", fps, sim_time);
    glfwSetWindowTitle(window_, title);
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
