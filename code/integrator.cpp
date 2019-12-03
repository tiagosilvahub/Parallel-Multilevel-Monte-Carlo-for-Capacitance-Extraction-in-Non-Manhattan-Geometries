
#define PI 3.1415926535897
#define OUTERLAYER 0
using namespace std;

// p original point, returns point after step of size d. delta is hitting distance
// electrode is the identifing number of the surface nearest to the point, to which the distance is calculated
// rndAngle = from -0.5 to 0.5 * 2pi

void integrator(Point *p, double d, double normalx, double normaly, double rndAngle, double delta){

	if(p->firstL == JUSTCOLLIDED){
		if(p->electrode != OUTERLAYER )
			return;			// performance increase return
		else { 				// particle is reflected by outer layer				
			p->firstL = NOCOLLISION;			
			if(d==0){		// it was possible for the particle to get struck on the outerlayer with double precision... d=0 was undefined behaviour
							// delta is needed in this function so that a reflecting step can be performed					
				p->x -= 2 * delta * normalx;
				p->y -= 2 * delta * normaly;
			} else {
				p->x -= 2 * abs(d) * normalx;	// p = p - 2|d|N 
				p->y -= 2 * abs(d) * normaly;			
			}
		}
	} else {
		p->x += cos(rndAngle)*d;
		p->y += sin(rndAngle)*d;
	}
}
