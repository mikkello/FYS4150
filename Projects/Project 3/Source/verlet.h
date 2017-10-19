#ifndef VERLET_H
#define VERLET_H


class Verlet
{
public:
	double m_dt;
	Verlet(double dt);
	void VelocityVerlet(class SolarSystem &system);
};

#endif // EULER_H
