#include "solarsystem.h"
#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
}

void SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
}

// Newtonian gravitational force calculation
void SolarSystem::calculateForcesAndEnergy()
{
	m_kineticEnergy = 0; m_potentialEnergy = 0; m_angularMomentum.zeros();

	for (CelestialBody &body : m_bodies) {
		body.force.zeros(); // Reset forces
	}

	for (int i = 0; i<numberOfBodies(); i++) {
		CelestialBody &body1 = m_bodies[i];
		for (int j = i + 1; j<numberOfBodies(); j++) {
			CelestialBody &body2 = m_bodies[j];
			vec3 deltaRVector = body1.position - body2.position;
			double dr = deltaRVector.length();
			vec3 force = -G*body1.mass*body2.mass* deltaRVector / pow(dr, 4);
			body1.force += force; body2.force -= force;
			m_potentialEnergy -= G*body1.mass*body2.mass / dr;
		}
		m_angularMomentum += body1.mass*body1.position.cross(body1.velocity);
		m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
	}
    // Remove comment to make the Sun stationary
	// vec3 force_sun = vec3(0, 0, 0);
    // m_bodies[0].force = force_sun;
}



// Relativistic correction to gravitational force
void SolarSystem::calculateRelativisticForcesAndEnergy()
{
	m_kineticEnergy = 0;
	m_potentialEnergy = 0;
	m_angularMomentum.zeros();

	CelestialBody &body1 = m_bodies[0]; CelestialBody &body2 = m_bodies[1];

	vec3 deltaRVector = body1.position - body2.position;
	double dr = deltaRVector.length();

	// Angular momentum
	double l1 = body1.position.cross(body1.velocity - body2.velocity).length();
	double l2 = body2.position.cross(body2.velocity - body1.velocity).length();

	// Relativistiv force calculation
	body1.force = -G*body1.mass*body2.mass*deltaRVector / (dr*dr*dr) * (1 + 3 * l1*l1 / (dr*dr*c*c));
	body2.force = G*body2.mass*body1.mass*deltaRVector / (dr*dr*dr) * (1 + 3 * l2*l2 / (dr*dr*c*c));

	m_angularMomentum = body1.position.cross(body1.momentum) + body2.position.cross(body2.momentum);
	m_potentialEnergy = -G*body1.mass*body2.mass / dr;
	m_kineticEnergy = 0.5*body1.mass*body1.velocity.lengthSquared() + 0.5*body2.mass*body2.velocity.lengthSquared();

	// Stationary Sun
	vec3 force_sun = vec3(0, 0, 0);
	m_bodies[0].force = force_sun;
}



int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalE() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialE() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticE() const
{
    return m_kineticEnergy;
}

vec3 SolarSystem::angularMomentum() const
{
	return m_angularMomentum;
}

void SolarSystem::writeToFile(string filename)
{
	if (!m_file.good()) {
		m_file.open(filename.c_str(), ofstream::out);
		if (!m_file.good()) {
			cout << "Error opening file " << filename << ". Aborting!" << endl;
			terminate();
		}
	}
	
	// Writes position of bodies, energies and angular momentum to file positionsx.xyz
	for (CelestialBody &body : m_bodies) {
		m_file << setprecision(15) << body.position.x() << " " << body.position.y() << " " << body.position.z() << " " << totalE() << " " << potentialE() << " " << kineticE() << " " << angularMomentum()[2] << "\n";
	}
}


std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}
