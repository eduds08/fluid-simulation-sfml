#include "FluidField.h"

int main()
{
    sf::RenderWindow window(sf::VideoMode(SIZE * SCALE, SIZE * SCALE), "Fluid Simulation");

    FluidField fluidField{ SIZE, 0.01f, 0.f, 0.f};

    sf::Vector2i currentPosition = sf::Mouse::getPosition(window);
    sf::Vector2i previousPosition = currentPosition;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        if (currentPosition != sf::Mouse::getPosition(window))
        {
            previousPosition = currentPosition;
            currentPosition = sf::Mouse::getPosition(window);
        }

        fluidField.addDensity(currentPosition.x / SCALE, currentPosition.y / SCALE, 100);

        float amtX = static_cast<float>(currentPosition.x - previousPosition.x);
        float amtY = static_cast<float>(currentPosition.y - previousPosition.y);

        fluidField.addVelocity(currentPosition.x / SCALE, currentPosition.y / SCALE, amtX, amtY);

        fluidField.update();

        window.clear(sf::Color::Cyan);

        fluidField.fade();
        fluidField.render(window);

        window.display();
    }

    return 0;
}
