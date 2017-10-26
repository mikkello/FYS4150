#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet(double dt) :
	m_dt(dt)
{

}

// Velocity Verlet algorithm
void Verlet::VelocityVerlet(SolarSystem &system)
{
	system.calculateForcesAndEnergy();
	for (CelestialBody &body : system.bodies())
	{
		body.velocity += 0.5*m_dt*(body.force / body.mass);   
		body.position += m_dt*body.velocity;                
		system.calculateForcesAndEnergy();                  
		body.velocity += 0.5*m_dt*(body.force / body.mass);  
	}
}

// Relativistic Velocity Verlet Algorithm
void Verlet::relativistic_VelocityVerlet(SolarSystem &system)
{
	system.calculateRelativisticForcesAndEnergy();

	for (CelestialBody &body : system.bodies()) {
		body.velocity += 0.5*(body.force / body.mass)*m_dt;
		body.position += body.velocity*m_dt;               
		system.calculateRelativisticForcesAndEnergy();
		body.velocity += 0.5*(body.force / body.mass)*m_dt;
	}
}