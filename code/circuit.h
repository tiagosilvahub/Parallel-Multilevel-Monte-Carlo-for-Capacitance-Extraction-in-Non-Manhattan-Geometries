#include <vector>

#define NOELECTRODE -1			// particle state
#define NOCOLLISION 0
#define JUSTCOLLIDED 1
#define COLLIDED 2

struct Point {		// each point should have DIM coordinates to support 3D, only 2D for now
	double x;
	double y;
	int firstL;	// NOCOLLISION JUSTCOLLIDED COLLIDED
	int secondL; 	// NOCOLLISION COLLIDED
	int electrode;	// NOELECTRODE:-1 OUTERLAYER:0 1 2 3 4 ...
};

struct Rectangle {		
	double electrode;		// electrode number the rectangle belongs to 
	double voltage;
	Point p1;				// rectangle defined by two opposing points, or in another perpective, 4 coordinates
	Point p2;
	std::vector<double> corners;	// radius of the rounded corners, can be 0 for right angle corners
									// vector indexes: 
									// 0 _ 1
									//  | |
									//  |_|
									// 3   2
};

struct Circuit {	// rectangles* and circle shapped outer layer
	int nrectangles;
	std::vector<Rectangle> rectangles;
	Point centre;
	double radius;
};

struct Vector {		
	double x;
	double y;
};
