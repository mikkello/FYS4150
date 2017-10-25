#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet(double dt) :
	m_dt(dt)
{

}

void Verlet::VelocityVerlet(SolarSystem &system)
{
	system.calculateForcesAndEnergy();

	for (CelestialBody &body : system.bodies())
	{
		body.position += m_dt*body.velocity + m_dt*m_dt / 2 * body.force / body.mass;
		vec3 acc = body.force / body.mass;
		system.calculateForcesAndEnergy();
		body.velocity += m_dt / 2 * (body.force / body.mass + acc);
	}
}

void Verlet::Mercury_VelocityVerlet(SolarSystem &system)
{
	system.calculateForces_GR();
	for (CelestialBody &body : system.bodies()) {
		body.velocity += 0.5*m_dt*body.force / body.mass;
		body.position += m_dt*body.velocity;
	}
	system.calculateForces_GR();

	for (CelestialBody &body : system.bodies()) {
		body.velocity += 0.5*m_dt*body.force / body.mass;
	}
}