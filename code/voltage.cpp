
double voltage(Point *p, Circuit circuit){

	return circuit.rectangles[p->electrode-1].voltage;

};