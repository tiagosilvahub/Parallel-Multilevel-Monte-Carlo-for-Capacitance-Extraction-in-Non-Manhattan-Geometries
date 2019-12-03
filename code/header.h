#include "distance.cpp"
#include "integrator.cpp"
#include "voltage.cpp"
#include "points.cpp"

double distance(Point *p, Circuit circuit, Vector *normal);
void integrator(Point *p, double d, double normalx, double normaly, double delta, double rndAngle);
double voltage(Point *p, Circuit circuit);
void points(Circuit circuit, Point * integrationPoints, double guassDelta);