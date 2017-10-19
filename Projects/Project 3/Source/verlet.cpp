#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet(double dt) :
	m_dt(dt)
{

}

void Verlet::VelocityVerlet(SolarSystem &system)
{
	system.calculateForcesAndEnergy();

	for (CelestialBody &body : system.bodies()) {
		body.velocity += m_dt / 2 * body.force / body.mass;
		body.position += m_dt*body.velocity;
	}
	system.calculateForcesAndEnergy();

	for (CelestialBody &body : system.bodies()) {
		body.velocity += m_dt / 2 * body.force / body.mass;
	}
}