
#include <algorithm>
#include <cmath>
#include <assert.h>

#include "circuit.h"

#define OUTERLAYER 0

using namespace std;

double rectangle(Point *p, Rectangle r){

	double x = p->x;
	double y = p->y;

	if(r.p2.y-r.corners[0] < y){
		if(x < r.p1.x+r.corners[0]){ // quadrant 1
			
			return sqrt( pow(x-(r.p1.x+r.corners[0]),2) + pow(y-(r.p2.y-r.corners[0]),2)) - r.corners[0];

		} else if(r.p1.x +r.corners[0] < x && x < r.p2.x-r.corners[1]){ // quadrant 2
			return abs(y-r.p2.y);

		} else if(r.p2.x-r.corners[1] < x){ // quadrant 3
			return sqrt( pow(x-(r.p2.x-r.corners[1]),2) + pow(y-(r.p2.y-r.corners[1]),2)) - r.corners[1];

		}
	} else if(r.p1.y+r.corners[3] < y && y < r.p2.y -r.corners[0]){ 
		if(x < r.p1.x+r.corners[0]){ // quadrant 4
			return abs(x-r.p1.x);

		} else if(r.p1.x +r.corners[0] < x && x < r.p2.x-r.corners[1]){ // quadrant 5  --  INSIDE circuit\t don't need distance? 
			return -1 * min( abs(x-r.p1.x), min(abs(x-r.p2.x), min(abs(y-r.p1.y), abs(y-r.p2.y))));

		} else if(r.p2.x-r.corners[1] < x){ // quadrant 6
			return abs(x-r.p2.x);

		}
	} else if(y < r.p1.y+r.corners[3]){
		if(x < r.p1.x+r.corners[0]){ // quadrant 7			
			return sqrt( pow(x-(r.p1.x+r.corners[3]),2) + pow(y-(r.p1.y+r.corners[3]),2)) - r.corners[3];

		} else if(r.p1.x +r.corners[0] < x && x < r.p2.x-r.corners[1]){ // quadrant 8
			return abs(y-r.p1.y);

		} else if(r.p2.x-r.corners[1] < x){ // quadrant 9
			return sqrt( pow(x-(r.p2.x-r.corners[2]),2) + pow(y-(r.p1.y+r.corners[2]),2)) - r.corners[2];

		}
	}
	return -1;
}


double distance(Point *p, double delta, double delta2, Circuit c, Vector * normal){	

	double distancePtoCentre;
	double dl,dr,dr2;
	int closestRect=0;

	// distance to outer layer is equal to radius - distance between the point and the centre of the outer layer
	distancePtoCentre = sqrt( pow( p->x - c.centre.x, 2.0) + pow( p->y - c.centre.y, 2.0));
	dl = c.radius - distancePtoCentre;

	// distance to all rectangles in circuit
	dr = rectangle(p, c.rectangles[0]);

	for(int i=1; i< c.nrectangles; i++){
		dr2 = rectangle(p, c.rectangles[i]);
		if(dr > dr2){
			dr = dr2;
			closestRect=i;
		}
	}

	if( dl < dr){
		p->electrode = OUTERLAYER;
		normal->x = (p->x - c.centre.x) / distancePtoCentre; 	// normal vector on the outerlayer
		normal->y = (p->y - c.centre.y) / distancePtoCentre;
		
		if(dl < delta)
			p->firstL = JUSTCOLLIDED; // with outerlayer

		return dl;	
	} else {
		p->electrode = c.rectangles[closestRect].electrode;		
		normal->x = 0;
		normal->y = 0;
		
		if(dr < delta && p->firstL != COLLIDED)
			p->firstL = JUSTCOLLIDED;

		if(dr < delta2)
			p->secondL=COLLIDED;

		return dr;
	} 
}