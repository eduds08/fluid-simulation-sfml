#pragma once

#include <SFML/Graphics.hpp>
#include <memory>
#include <vector>
#include "Constants.h"

using namespace constants;

class FluidField
{
public:

    FluidField(int size, float dt, float diffusion, float viscosity);
    ~FluidField() {}

    void addDensity(int x, int y, float amount);

    void addVelocity(int x, int y, float amountX, float amountY);

    void setBoundaries(int b, std::vector<float>& x);

    void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c, int iter);

    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, int iter);

    void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div, int iter);

    void advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY);

    void update();

    void render(sf::RenderWindow& window);

    void fade();

    const int xy(int x, int y) const 
    { 
        if (x < 0) x = 0;
        if (y < 0) y = 0;
        if (x >= m_size) x = m_size - 1;
        if (y >= m_size) y = m_size - 1;
        
        return x + y * m_size;
    }

private:
    sf::RectangleShape m_shape{ sf::Vector2f{ SCALE, SCALE } };

    int m_size{};
    float m_dt{};
    float m_diff{};
    float m_visc{};

    std::vector<float> m_s{};
    std::vector<float> m_density{};

    std::vector<float> m_Vx{};
    std::vector<float> m_Vy{};

    std::vector<float> m_Vx0{};
    std::vector<float> m_Vy0{};
};
