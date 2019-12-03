#include <iostream>
#include <stdio.h>
#include <limits>
#include <fstream>

#define POINTSTAKEN 9

void points(Circuit circuit, Point * integrationPoints, double guassDelta){

		int nElectrodes = 9;
		int * nPointsElectrode = new int[nElectrodes];

		nPointsElectrode[0] = 4;		// TODO: from input
		nPointsElectrode[1] = 6;
		nPointsElectrode[2] = 8;
		nPointsElectrode[3] = 4;
		nPointsElectrode[4] = 4;
		nPointsElectrode[5] = 6;
		nPointsElectrode[6] = 6;
		nPointsElectrode[7] = 4;
		nPointsElectrode[8] = 6;

		Point * pointsElectrode = new Point[48];
		
		pointsElectrode[0].x = circuit.rectangles[0].p1.x;
		pointsElectrode[0].y = circuit.rectangles[0].p2.y;
		pointsElectrode[1].x = circuit.rectangles[0].p2.x;
		pointsElectrode[1].y = circuit.rectangles[0].p2.y;
		pointsElectrode[2].x = circuit.rectangles[0].p2.x;
		pointsElectrode[2].y = circuit.rectangles[0].p1.y;
		pointsElectrode[3].x = circuit.rectangles[0].p1.x;
		pointsElectrode[3].y = circuit.rectangles[0].p1.y;
		
		pointsElectrode[4].x = circuit.rectangles[2].p1.x;
		pointsElectrode[4].y = circuit.rectangles[2].p2.y;
		pointsElectrode[5].x = circuit.rectangles[1].p1.x;
		pointsElectrode[5].y = circuit.rectangles[1].p1.y;
		pointsElectrode[6].x = circuit.rectangles[1].p1.x;
		pointsElectrode[6].y = circuit.rectangles[1].p2.y;
		pointsElectrode[7].x = circuit.rectangles[1].p2.x;
		pointsElectrode[7].y = circuit.rectangles[1].p2.y;
		pointsElectrode[8].x = circuit.rectangles[2].p2.x;
		pointsElectrode[8].y = circuit.rectangles[2].p1.y;
		pointsElectrode[9].x = circuit.rectangles[2].p1.x;
		pointsElectrode[9].y = circuit.rectangles[2].p1.y;
		
		pointsElectrode[10].x = circuit.rectangles[3].p1.x;
		pointsElectrode[10].y = circuit.rectangles[3].p2.y;
		pointsElectrode[11].x = circuit.rectangles[3].p2.x;
		pointsElectrode[11].y = circuit.rectangles[3].p2.y;
		pointsElectrode[12].x = circuit.rectangles[3].p2.x;
		pointsElectrode[12].y = circuit.rectangles[3].p1.y;	
		pointsElectrode[13].x = circuit.rectangles[4].p2.x;
		pointsElectrode[13].y = circuit.rectangles[4].p2.y;
		pointsElectrode[14].x = circuit.rectangles[4].p2.x;
		pointsElectrode[14].y = circuit.rectangles[4].p1.y;
		pointsElectrode[15].x = circuit.rectangles[4].p1.x;
		pointsElectrode[15].y = circuit.rectangles[4].p1.y;	
		pointsElectrode[16].x = circuit.rectangles[4].p1.x;
		pointsElectrode[16].y = circuit.rectangles[4].p2.y;
		pointsElectrode[17].x = circuit.rectangles[3].p1.x;
		pointsElectrode[17].y = circuit.rectangles[3].p1.y;

		pointsElectrode[18].x = circuit.rectangles[5].p1.x;
		pointsElectrode[18].y = circuit.rectangles[5].p2.y;
		pointsElectrode[19].x = circuit.rectangles[5].p2.x;
		pointsElectrode[19].y = circuit.rectangles[5].p2.y;
		pointsElectrode[20].x = circuit.rectangles[5].p2.x;
		pointsElectrode[20].y = circuit.rectangles[5].p1.y;
		pointsElectrode[21].x = circuit.rectangles[5].p1.x;
		pointsElectrode[21].y = circuit.rectangles[5].p1.y;

		pointsElectrode[22].x = circuit.rectangles[6].p1.x;
		pointsElectrode[22].y = circuit.rectangles[6].p2.y;
		pointsElectrode[23].x = circuit.rectangles[6].p2.x;
		pointsElectrode[23].y = circuit.rectangles[6].p2.y;
		pointsElectrode[24].x = circuit.rectangles[6].p2.x;
		pointsElectrode[24].y = circuit.rectangles[6].p1.y;
		pointsElectrode[25].x = circuit.rectangles[6].p1.x;
		pointsElectrode[25].y = circuit.rectangles[6].p1.y;

		pointsElectrode[26].x = circuit.rectangles[7].p1.x;
		pointsElectrode[26].y = circuit.rectangles[7].p2.y;
		pointsElectrode[27].x = circuit.rectangles[7].p2.x;
		pointsElectrode[27].y = circuit.rectangles[7].p2.y;
		pointsElectrode[28].x = circuit.rectangles[7].p2.x;
		pointsElectrode[28].y = circuit.rectangles[7].p1.y;
		pointsElectrode[29].x = circuit.rectangles[8].p2.x;
		pointsElectrode[29].y = circuit.rectangles[8].p2.y;
		pointsElectrode[30].x = circuit.rectangles[8].p2.x;
		pointsElectrode[30].y = circuit.rectangles[8].p1.y;
		pointsElectrode[31].x = circuit.rectangles[8].p1.x;
		pointsElectrode[31].y = circuit.rectangles[8].p1.y;

		pointsElectrode[32].x = circuit.rectangles[10].p1.x;
		pointsElectrode[32].y = circuit.rectangles[10].p2.y;
		pointsElectrode[33].x = circuit.rectangles[10].p2.x;
		pointsElectrode[33].y = circuit.rectangles[10].p2.y;
		pointsElectrode[34].x = circuit.rectangles[9].p2.x;
		pointsElectrode[34].y = circuit.rectangles[9].p1.y;
		pointsElectrode[35].x = circuit.rectangles[9].p1.x;
		pointsElectrode[35].y = circuit.rectangles[9].p1.y;	
		pointsElectrode[36].x = circuit.rectangles[9].p1.x;
		pointsElectrode[36].y = circuit.rectangles[9].p2.y;	
		pointsElectrode[37].x = circuit.rectangles[10].p1.x;
		pointsElectrode[37].y = circuit.rectangles[10].p1.y;

		pointsElectrode[38].x = circuit.rectangles[11].p1.x;
		pointsElectrode[38].y = circuit.rectangles[11].p2.y;
		pointsElectrode[39].x = circuit.rectangles[11].p2.x;
		pointsElectrode[39].y = circuit.rectangles[11].p2.y;
		pointsElectrode[40].x = circuit.rectangles[11].p2.x;
		pointsElectrode[40].y = circuit.rectangles[11].p1.y;
		pointsElectrode[41].x = circuit.rectangles[11].p1.x;
		pointsElectrode[41].y = circuit.rectangles[11].p1.y;

		pointsElectrode[42].x = circuit.rectangles[12].p1.x;
		pointsElectrode[42].y = circuit.rectangles[12].p2.y;
		pointsElectrode[43].x = circuit.rectangles[12].p2.x;
		pointsElectrode[43].y = circuit.rectangles[12].p2.y;
		pointsElectrode[44].x = circuit.rectangles[12].p2.x;
		pointsElectrode[44].y = circuit.rectangles[12].p1.y;
		pointsElectrode[45].x = circuit.rectangles[13].p2.x;
		pointsElectrode[45].y = circuit.rectangles[13].p2.y;
		pointsElectrode[46].x = circuit.rectangles[13].p2.x;
		pointsElectrode[46].y = circuit.rectangles[13].p1.y;
		pointsElectrode[47].x = circuit.rectangles[13].p1.x;
		pointsElectrode[47].y = circuit.rectangles[13].p1.y;


		// calculate the perimeter
		double perimeter;
		double distanceP;
		Point previousP, nextC, beforeC;	// previousP = either the last integration point saved or the last corner turned 
											// nextC = the nextC corner
											// beforeC = hte corner before 
		int pe=0;
		int t;
		
		// for all electrodes, calculates their perimeter. divides by the number of points we need to get distanceP. then calculates the points
		for(int e=0; e<nElectrodes; e++){

			perimeter=0;
			
			for(int p =0; p<nPointsElectrode[e]-1; p++){

				if(pointsElectrode[ pe + p+1].x == pointsElectrode[ pe + p].x){
					perimeter += abs(pointsElectrode[ pe + p+1].y - pointsElectrode[ pe + p].y);
				} else {
					perimeter += abs(pointsElectrode[ pe + p+1].x - pointsElectrode[ pe + p].x);
				}
			}
			if(pointsElectrode[ pe + 0].x == pointsElectrode[ pe + nPointsElectrode[e]-1].x){
				perimeter += abs(pointsElectrode[ pe + 0].y - pointsElectrode[ pe + nPointsElectrode[e]-1].y);
			} else {
				perimeter += abs(pointsElectrode[ pe + 0].x - pointsElectrode[ pe + nPointsElectrode[e]-1].x);
			}
			//cout << "perimeter " << perimeter << "  interval " <<  perimeter /  POINTSTAKEN << endl;
			
			t=0;
			integrationPoints[e*POINTSTAKEN+0].x = pointsElectrode[pe].x - sqrt(2)*guassDelta;
			integrationPoints[e*POINTSTAKEN+0].y = pointsElectrode[pe].y + sqrt(2)*guassDelta;
			previousP.x = pointsElectrode[pe].x; 
			previousP.y = pointsElectrode[pe].y; 

			//cout << "Point 0" << " " << integrationPoints[e*POINTSTAKEN+0].x << " " << integrationPoints[e*POINTSTAKEN+0].y << endl;						

			//  follows the border clockwise, keeping track if it's going down/up/left/right, and the distance to the last integration point saved
			//  distanceP and previous can either refer to the last integration point or the last corner turned, which ever is closer
			//  this is because the saved point is the the last point + distance
			for(int pt =1; pt< POINTSTAKEN; pt++){		// pt - point taken 
				
				distanceP = perimeter /  POINTSTAKEN;	// distance between each point is 
														// currect point of the electrode
				while(distanceP!=0){
	 
					nextC.x = pointsElectrode[ pe + (t+1)%nPointsElectrode[e]].x;
					nextC.y = pointsElectrode[ pe + (t+1)%nPointsElectrode[e]].y;
					beforeC.x = pointsElectrode[ pe + t%nPointsElectrode[e]].x;
					beforeC.y = pointsElectrode[ pe + t%nPointsElectrode[e]].y;				
					
					if(nextC.x == beforeC.x){					
						if(nextC.y > beforeC.y){
	 																// going up - increasing y. integration surface is to the left
							if( previousP.y + distanceP > nextC.y){
								distanceP = distanceP - ( nextC.y - previousP.y);
								previousP.x = nextC.x;
								previousP.y = nextC.y;
								t++;					
							} else if( previousP.y + distanceP == nextC.y){
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x - sqrt(2)*guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y + distanceP + sqrt(2)*guassDelta;
								//previousP.x same
								previousP.y += distanceP;
								//cout <<  "A2 Point "<< pt << " " << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							} else {
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x  - guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y + distanceP;
								//previousP.x same
								previousP.y += distanceP;
								//cout <<  "A Point "<< pt << " " << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							}
									
						} else { 									// going down - decreasing y.  integration surface is to the right
							if( previousP.y - distanceP < nextC.y){
								distanceP = distanceP - (previousP.y - nextC.y); 
								previousP.x = nextC.x;
								previousP.y = nextC.y;					
								t++;
							} else if(previousP.y - distanceP == nextC.y){
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x + sqrt(2)*guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y - distanceP - sqrt(2)*guassDelta;
								//previousP.x same
								previousP.y -= distanceP;
								//cout << "B2 Point "<< pt << " "  << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							} else {
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x + guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y - distanceP;
								//previousP.x same
								previousP.y -= distanceP;
								//cout << "B Point "<< pt << " "  << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							}
						}	
					} else {
						if(nextC.x > beforeC.x){
																		// going right - increasing x.  integration surface is above
							if( previousP.x + distanceP > nextC.x){
								distanceP = distanceP - (nextC.x - previousP.x);
								previousP.x = nextC.x;
								previousP.y = nextC.y;					
								t++;
							} else if(previousP.x + distanceP == nextC.x){
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x + distanceP + sqrt(2)*guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y + sqrt(2)*guassDelta;
								previousP.x +=  distanceP;
								//previousP.y same
								//cout << "C2 Point "<< pt << " "  << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							} else {
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x + distanceP;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y  + guassDelta;
								previousP.x +=  distanceP;
								//previousP.y same
								//cout << "C Point "<< pt << " "  << integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							}		
					
						} else { 
																		// going left - decreasing x.  integration surface is below 	
							if( previousP.x - distanceP < nextC.x){
								distanceP = distanceP - (previousP.x - nextC.x);
								previousP.x = nextC.x;
								previousP.y = nextC.y;	 				
								t++;
							} else if(previousP.x - distanceP == nextC.x){
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x - distanceP - sqrt(2)*guassDelta;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y - sqrt(2)*guassDelta;
								previousP.x -= distanceP;
								//previousP.y same
								//cout <<  "D2 Point "<< pt << " " <<integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							} else {
								integrationPoints[e*POINTSTAKEN+pt].x = previousP.x - distanceP;
								integrationPoints[e*POINTSTAKEN+pt].y = previousP.y - guassDelta;
								previousP.x -= distanceP;
								//previousP.y same
								//cout <<  "D Point "<< pt << " " <<integrationPoints[e*POINTSTAKEN+pt].x << " " << integrationPoints[e*POINTSTAKEN+pt].y << endl;
								distanceP=0;
							}
						}

					}	
				}
			}

			pe += nPointsElectrode[e];
		}

		delete[] nPointsElectrode;
		delete[] pointsElectrode;

	}