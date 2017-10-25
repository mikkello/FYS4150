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

void SolarSystem::calculateForcesAndEnergy()
{
	m_kineticEnergy = 0;
	m_potentialEnergy = 0;
	m_angularMomentum.zeros();

	for (CelestialBody &body : m_bodies) {
		// Reset forces on all bodies
		body.force.zeros();
	}

	for (int i = 0; i<numberOfBodies(); i++) {
		CelestialBody &body1 = m_bodies[i];
		for (int j = i + 1; j<numberOfBodies(); j++) {
			CelestialBody &body2 = m_bodies[j];
			vec3 deltaRVector = body1.position - body2.position;
			double dr = deltaRVector.length();
			vec3 force = -G*body1.mass*body2.mass* deltaRVector / pow(dr, 3);
			body1.force += force;
			body2.force -= force;
			m_potentialEnergy += G*body1.mass*body2.mass / dr;
		}
		m_angularMomentum += body1.mass*body1.position.cross(body1.velocity);
		m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
	}
 // Stationary Sun
	//vec3 force_sun = vec3(0, 0, 0);
	//m_bodies[0].force = force_sun;
}

void SolarSystem::calculateForces_GR() {

	for (CelestialBody &body : m_bodies) {
		body.force.zeros();
	}

	for (int i = 0; i<numberOfBodies(); i++) {
		CelestialBody &body1 = m_bodies[i];
		for (int j = i + 1; j<numberOfBodies(); j++) {
			CelestialBody &body2 = m_bodies[j];
			vec3 deltaRVector = body1.position - body2.position;
			double dr = deltaRVector.length();
			vec3 l_vec = body1.position.cross(body1.velocity);
			double l = l_vec.length();
			double c = 63239.726;
			vec3 force = -G*body1.mass / (dr*dr*dr)*deltaRVector*(1 + 3 * l * l / (dr*dr*c*c));
			body1.force += force;
			body2.force -= force;
		}
	}
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
		m_file << setprecision(15) << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n"
//<< " " << totalE() << " " << potentialE() << " " << kineticE() << " " << angularMomentum()[2] << "\n";
	}
}




std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}
