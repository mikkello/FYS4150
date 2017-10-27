#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#define _USE_MATH_DEFINES

#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <math.h>


class SolarSystem
{
public:
    SolarSystem();
	double G = 4 * M_PI*M_PI; // Gravitational constant [4*pi^2*AU^3 / M*yr^2]
	double c = 63198; // Speed of light [AU/yr]

    void createCelestialBody(vec3 position, vec3 velocity, double mass);
    void calculateForcesAndEnergy();
	void SolarSystem::calculateRelativisticForcesAndEnergy();
	vec3 SolarSystem::angularMomentum() const;
	int numberOfBodies() const;
    double totalE() const;
    double potentialE() const;
    double kineticE() const;
    void writeToFile(std::string filename);
    std::vector<CelestialBody> &bodies();

private:
    std::vector<CelestialBody> m_bodies;
    std::ofstream m_file;
	vec3 m_angularMomentum;
    double m_kineticEnergy;
    double m_potentialEnergy;
};

#endif // SOLARSYSTEM_H
