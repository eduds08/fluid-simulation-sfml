#include "FluidField.h"

FluidField::FluidField(int size, float dt, int diffusion, int viscosity)
{
    m_size = size;
    m_dt = dt;
    m_diff = static_cast<float>(diffusion);
    m_visc = static_cast<float>(viscosity);

    for (int i = 0; i < size * size; ++i)
    {
        m_s.emplace_back(0.f);
        m_density.emplace_back(0.f);
        m_Vx.emplace_back(0.f);
        m_Vy.emplace_back(0.f);
        m_Vx0.emplace_back(0.f);
        m_Vy0.emplace_back(0.f);
    }
}

void FluidField::addDensity(int x, int y, float amount)
{
    m_density[xy(x, y)] += amount;
}

void FluidField::addVelocity(int x, int y, float amountX, float amountY)
{
    m_Vx[xy(x, y)] += amountX;
    m_Vy[xy(x, y)] += amountY;
}

void FluidField::setBoundaries(int b, std::vector<float>& x)
{
    for (int i = 1; i < m_size - 1; ++i)
    {
        x[xy(i, 0)] = b == 2 ? -x[xy(i, 1)] : x[xy(i, 1)];
        x[xy(i, m_size - 1)] = b == 2 ? -x[xy(i, m_size - 2)] : x[xy(i, m_size - 2)];

        x[xy(0, i)] = b == 1 ? -x[xy(1, i)] : x[xy(1, i)];
        x[xy(m_size - 1, i)] = b == 1 ? -x[xy(m_size - 2, i)] : x[xy(m_size - 2, i)];
    }

    x[xy(0, 0)] = 0.33f * (x[xy(1, 0)] + x[xy(0, 1)]);

    x[xy(0, m_size - 1)] = 0.33f * (x[xy(1, m_size - 1)] + x[xy(0, m_size - 2)]);

    x[xy(m_size - 1, 0)] = 0.33f * (x[xy(m_size - 2, 0)] + x[xy(m_size - 1, 1)]);

    x[xy(m_size - 1, m_size - 1)] = 0.33f * (x[xy(m_size - 2, m_size - 1)] + x[xy(m_size - 1, m_size - 2)]);
}

void FluidField::lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c, int iter)
{
    float cRecip = 1.f / c;

    for (int k = 0; k < iter; ++k) 
    {
        for (int j = 1; j < m_size - 1; ++j)
        {
            for (int i = 1; i < m_size - 1; ++i)
            {
                x[xy(i, j)] = (x0[xy(i, j)] + a * (x[xy(i + 1, j)] + x[xy(i - 1, j)] + x[xy(i, j + 1)] + x[xy(i, j - 1)])) * cRecip;
            }
        }
        setBoundaries(b, x);
    }
}

void FluidField::diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, int iter)
{
    float a = m_dt * diff * (m_size - 2) * (m_size - 2);

    lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}

void FluidField::project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div, int iter)
{
    for (int j = 1; j < m_size - 1; ++j)
    {
        for (int i = 1; i < m_size - 1; ++i)
        {
            div[xy(i, j)] = -0.5f * (velocX[xy(i + 1, j)] - velocX[xy(i - 1, j)] + velocY[xy(i, j + 1)] - velocY[xy(i, j - 1)]) / m_size;
            p[xy(i, j)] = 0;
        }
    }

    setBoundaries(0, div);
    setBoundaries(0, p);
    lin_solve(0, p, div, 1, 6, iter);

    for (int j = 1; j < m_size - 1; ++j)
    {
        for (int i = 1; i < m_size - 1; ++i)
        {
            velocX[xy(i, j)] -= 0.5f * (p[xy(i + 1, j)] - p[xy(i - 1, j)]) * m_size;
            velocY[xy(i, j)] -= 0.5f * (p[xy(i, j + 1)] - p[xy(i, j - 1)]) * m_size;
        }
    }

    setBoundaries(1, velocX);
    setBoundaries(2, velocY);
}

void FluidField::advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY)
{
    float i0, i1, j0, j1;

    float dtx = m_dt * (m_size - 2);
    float dty = m_dt * (m_size - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    int i, j;

    for (j = 1; j < m_size - 1; ++j)
    {
        for (i = 1; i < m_size - 1; ++i)
        {
            tmp1 = dtx * velocX[xy(i, j)];
            tmp2 = dty * velocY[xy(i, j)];

            x = static_cast<float>(i) - tmp1;
            y = static_cast<float>(j) - tmp2;

            if (x < 0.5f) 
            {
                x = 0.5f;
            }

            if (x > static_cast<float>(m_size) + 0.5f) 
            {
                x = static_cast<float>(m_size) + 0.5f;
            }

            i0 = floorf(x);
            i1 = i0 + 1.0f;

            if (y < 0.5f)
            {
                y = 0.5f;
            }

            if (y > static_cast<float>(m_size) + 0.5f)
            {
                y = static_cast<float>(m_size) + 0.5f;
            }

            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = static_cast<int>(i0);
            int i1i = static_cast<int>(i1);
            int j0i = static_cast<int>(j0);
            int j1i = static_cast<int>(j1);

            d[xy(i, j)] =
                s0 * (t0 * d0[xy(i0i, j0i)] + (t1 * d0[xy(i0i, j1i)])) +
                s1 * (t0 * d0[xy(i1i, j0i)] + (t1 * d0[xy(i1i, j1i)]));
        }
    }

    setBoundaries(b, d);
}

void FluidField::update()
{
    diffuse(1, m_Vx0, m_Vx, m_visc, 4);
    diffuse(2, m_Vy0, m_Vy, m_visc, 4);

    project(m_Vx0, m_Vy0, m_Vx, m_Vy, 4);

    advect(1, m_Vx, m_Vx0, m_Vx0, m_Vy0);
    advect(2, m_Vy, m_Vy0, m_Vx0, m_Vy0);

    project(m_Vx, m_Vy, m_Vx0, m_Vy0, 4);

    diffuse(0, m_s, m_density, m_diff, 4);
    advect(0, m_density, m_s, m_Vx, m_Vy);
}

void FluidField::render(sf::RenderWindow& window)
{
    for (int i = 0; i < m_size; ++i)
    {
        for (int j = 0; j < m_size; ++j)
        {
            float x = static_cast<float>(i * SCALE);
            float y = static_cast<float>(j * SCALE);

            float d = m_density[xy(i, j)];

            m_shape.setPosition(x, y);

            m_shape.setFillColor(sf::Color{ 255, 255, 255, static_cast<sf::Uint8>(d) });

            window.draw(m_shape);
        }
    }
}

void FluidField::fade()
{
    for (int i = 0; i < m_density.size(); ++i)
    {
        float d = m_density[i];
        if (d < 0.f) d = 0.f;
        if (d > 255.f) d = 255.f;
        m_density[i] = d - 0.01f;
    }
}
