/******
*
*   Tiago Silva
*   Main code of Master Thesis:
*   Parallel Multilevel Monte Carlo Algorithm for 
*   Capacitance Extraction in Non-Manhattan Geometries 
*
*/

#include <math.h>
#include <iostream>
#include <stdio.h>
#include "header.h"
#include <mpi.h>

#define NWARMUP 50000
// max level of descritization, theoretically would never need something higher
#define LM 30
#define ERRORSPLIT 0.75 // division between erros, bias² and variance

int main(int argc, char * argv[]){

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* get myid and # of processors */
  int myid,numproc;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  std::cout.precision(10);

  	// total running time start
	double algStartTime = MPI_Wtime();

  // Initializations
	// particle point
	Point * p = new Point;
	// normal vector of particle point and reflection boundary that surrounds the circuit
	Vector * normal = new Vector;
	normal->x = 0;
	normal->y = 0;
	// target accuracy
	double targetAcc;
	// geometry of the circuit being analysed 
	Circuit circuit;
	// Maximum number of multi levels
	int maxLevel;
	// the distance at which a particle is considered to 'hit' an electrode. the discretization
	double D[LM];    
	// distance between guassian integration surface and eletrodes
	double guassianD;    
	// number of integration points per electrode
	int nintegrationpoints;
	// number of simulations to run per level for estimation of variance and cost 
	int N0=NWARMUP;
	// number of simulations per level
	int NumSimsToDo[LM];
	double NumSimsDone[LM];
	// cost of a simulation on level l    
	double cost[LM];   
	// variance at each level, based on average from all simulations at that level
	double var[LM];
	// voltage at each level, based on average from all simulations at that level
	double volt[LM];
	// constant associated with bias, experimental result 
	double beta=0.776;
	// variable for counting the time passed for cost assessment
	double starttime; 
	// sum of voltages, local= each process. tmp temporary variable for mpi reduce
	double sumVoltage[LM];    
	double localSum[LM];
	double tmpsumVoltage[LM];    
	// sum of the square of voltages, local= each process. tmp temporary variable for mpi reduce
	double sumSqVoltage[LM];    
	double localSqSum[LM];        
	double tmpsumSqVoltage[LM];     
	// cost of each level, based on average from all simulations at that level
	double globalcost[LM];
	// // local= each process. tmp temporary variable for mpi reduce
	double localcost[LM];
	double tmpcost[LM];	
	//  simulations left to run, difference between todo and done
	int dueSimulations;
	// Variance*Cost sum for equation 
	double eqSum;       
	// voltage of the electrode hit with first and second descritization
	double v, v2;
	// distance from particle to closest electrode
	double dist;
	// M defines the decrease in delta at every level, typically 4 is optimal.    
	int M;  
	// loop variable, simulations per process
	int simPerProc;

// setup random number generator to make random angles between 0 and 2*Pi
	  std::random_device rd;
	  std::mt19937 gen(rd());
	  std::uniform_real_distribution<> dis(0, 2.0 * PI);


// input arguments from argv    
	int arg=1;
	// initial point    
	p->x = atof(argv[arg++]);
	p->y = atof(argv[arg++]);
	// degugging option
	bool debug = argv[arg++];

	// target accuracy
	targetAcc = atof(argv[arg++]);

	// maximum level
	maxLevel = atoi(argv[arg++]);

	// M defines the decrease in delta at every level, typically 4 is optimal.
	M = atoi(argv[arg++]);

	guassianD = atof(argv[arg++]);

	nintegrationpoints  = atoi(argv[arg++]);
  
	// geometry of the circuit    
	  
	  // coordinates of the center of the circuit
	  circuit.centre.x = atof(argv[arg++]);
	  circuit.centre.y = atof(argv[arg++]);
	  
	  // radius of the outer boundary we consider to be the end of the circuit
	  circuit.radius = atof(argv[arg++]);

	  // number of rectangles that make up the circuit
	  // one electrode can be made of several rectangles
	  circuit.nrectangles = atoi(argv[arg++]);

	  if(argc != arg + 6* circuit.nrectangles){ // 6 number of arguments per rectangle
		  if(myid==0)
			cout << "Incorrect number of arguments expected " << arg << " + " << 6* circuit.nrectangles << " received " << argc << endl;          
		  MPI_Finalize();
		  exit(1);
		}
		  
	  // circuit input. see circuit.h structures
	  Rectangle * r;
	  for(long i=0; i<circuit.nrectangles; i++){
		r = new Rectangle();

		r->p1.x = atof(argv[arg++]);
		r->p1.y = atof(argv[arg++]);
		r->p2.x = atof(argv[arg++]);
		r->p2.y = atof(argv[arg++]);
		r->corners.push_back( atof(argv[arg]));
		r->corners.push_back( atof(argv[arg]));
		r->corners.push_back( atof(argv[arg]));
		r->corners.push_back( atof(argv[arg++]));

		r->voltage = atof(argv[arg++]);
		r->electrode = i+1;       
		circuit.rectangles.push_back(*r);
	  }
  
  // initial calculations
	  
  // there are two sources of error: BIAS and VARIANCE, like so:
  // mean error <= beta*delta + 2*sqrt( alpha / NSimulations)
  
  // the target accuracy is the maximum acceptable amount of error
  // we let each source of error be half of the maximum error,   
  // the bias error relates directly to the discretization.   
  // half error for the variance relates to the number of simulations ran   
  // so half error for the bias:
	
	/* beta*delta = targetAcc / 2 <=> */
	  D[maxLevel+1] = targetAcc / ( 2*beta );

  for (int i = maxLevel; i >=0; i--){
	D[i] = D[i+1]*M;      // each level has M times smaller discretization
	/*if(myid==0 && debug)
	  cout <<" d[" << i << "]= " << D[i];*/
  }

  /*if(myid==0 && debug)
	cout << endl;
  */
  if(D[0] > guassianD){  // the discretization chosen is so big, the starting position of the particle is already hitting an electrode
	if(myid==0)
	  cout << "Number of levels too high for the accuracy requested" << endl;
	  MPI_Finalize(); 
	  exit(1);
  }
  
  int totalPoints = 9*nintegrationpoints; // TODO:  input electrodes
  Point * integrationPoints = new Point[totalPoints];

  if(nintegrationpoints==0){
	integrationPoints = new Point;
	totalPoints = 1;
	integrationPoints[0].x = p->x;
	integrationPoints[0].y = p->y;
  } else 
	points(circuit, integrationPoints, guassianD);

  for(int l=0; l<=maxLevel; l++) {
	NumSimsToDo[l]   = N0;			// initial number of runs
	NumSimsDone[l] 	= 0;
	globalcost[l]	= 0.0;	
	sumVoltage[l]   = 0.0;
	sumSqVoltage[l] = 0.0;
  }

  // begin algorithm
  do {

	for (int l=0; l<=maxLevel; l++) {
				
		localSum[l] 	= 0.0;
		localSqSum[l] 	= 0.0;
		tmpcost[l]		= 0.0;
		tmpsumVoltage[l] = 0.0;
		tmpsumSqVoltage[l] = 0.0;
	
	  if (NumSimsToDo[l]>0) {        

		simPerProc = NumSimsToDo[l]/numproc;    // divide work equally between processes
		starttime = MPI_Wtime();

		for(long i=0; i<simPerProc; i++){     // run NumSim[l] simulations, nwalks per process
		  p->x = integrationPoints[0].x;   // (x,y) where walk begins
		  p->y = integrationPoints[0].y; 
		  p->firstL = NOCOLLISION;      // first level of discretization collision
		  p->secondL = NOCOLLISION;     // second level 
		  p->electrode = NOELECTRODE;

		  while(p->secondL != COLLIDED){              // while particle isn't within delta distance of an electrode, 
			dist = distance(p, D[l], D[l+1], circuit, normal);    // calculate distance to nearest electrode
									// detects if particle hits electrode at each discretization                                  
			integrator(p, dist, normal->x, normal->y, // particle Walk On Spheres step of size dist
					 dis(gen), D[l]);
		  if(p->firstL == JUSTCOLLIDED){; // if particle just hit electrode in the last step
			  if(l==0){                   // if standard monte carlo
				p->secondL = COLLIDED;  // exit loop, record only 1 voltage of 1 electrode hit
				v2=0;
			  } else {                 // if multilevel monte carlo 
				v2 = voltage(p, circuit);       // save voltage of the electrode hit
				p->firstL = COLLIDED; // continue until particle collides with smaller discretization
			  } 
			}
		  } 
			
		  v = voltage(p, circuit);    // voltage of electrode touched
		  v = v - v2;           // difference between rough and fine simulation, notice v2=0 for standard MC
		
		  localSum[l]  += v;           // sum of voltages
		  localSqSum[l]+= pow(v,2);    // sum of the square of voltages, for calculating the variance          
		}

		localcost[l] = MPI_Wtime() - starttime;
		NumSimsDone[l] += simPerProc*numproc;
	  }
	}

	MPI_Allreduce(&localcost,&tmpcost,LM,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&localSum,&tmpsumVoltage,LM,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&localSqSum,&tmpsumSqVoltage,LM,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	// sum cost, voltage and sqVoltage per level, of this run of L levels, which the previous runs
	
	for (int l=0; l<=maxLevel; l++) {
		globalcost[l] += tmpcost[l];
		sumVoltage[l] += abs(tmpsumVoltage[l]);
		sumSqVoltage[l] += tmpsumSqVoltage[l];
	}

	eqSum = 0;
	// calculate at each level
	for (int l=0; l<=maxLevel; l++) {

	  // average voltage, variance, cost
	  volt[l] = abs(sumVoltage[l]/NumSimsDone[l]);
	  var[l] = fmax(sumSqVoltage[l]/NumSimsDone[l] - pow(volt[l],2), 0);
	  cost[l] = globalcost[l] / NumSimsDone[l];
	  // sum needed for number of simulations equation
	  eqSum += sqrt(var[l]*cost[l]);
	}
	
	// estimate number of simulations per level that should be run to reach the target accuracy
	for (int l=0; l<=maxLevel; l++){
	  NumSimsToDo[l] = ceil(eqSum*sqrt(var[l]/cost[l])/(ERRORSPLIT*pow(targetAcc,2)));
	  NumSimsToDo[l] -= NumSimsDone[l];
	  if(NumSimsToDo[l] < 0)
		NumSimsToDo[l] = 0;
	}
  
  // check if all needed simulations have been run. if less than 1% more are needed, ignore
	dueSimulations = 0;
	for (int l=0; l<=maxLevel; l++)
	  dueSimulations += fmax(0.0, (double)NumSimsToDo[l]-0.01*NumSimsDone[l]);

  } while(dueSimulations!=0);

  double result = 0;

  for (int l=0; l<=maxLevel; l++){
	result += volt[l];
	//if (debug && myid==0) 
	//  cout << " N" << l << ": " << NumSimsDone[l] << " C"<<l<<"l: " << cost[l]<<endl;   
  }
  if (myid==0) 
	cout << "result: " << result << " time: " << MPI_Wtime() - algStartTime << endl;
  
  MPI_Finalize();
}