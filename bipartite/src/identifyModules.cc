// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// identifyModules - detects modules in the graph (based on fitHRG of Aaron Clauset)
// Copyright (C) 2005-2009 Aaron Clauset
// Copyright (C) 2010-2011 Rouven Strauss
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 26 October 2005
// Modified     : many, many times
//		: 27 December 2007 (cleaned up for public consumption)
//
// Modified by	: Rouven Strauss
// Collaborators: Carsten F. Dormann
// Project      : Module Identification in Bipartite Networks
// Location     : University of Munich, Dept. of Computer Science AND Helmholtz Centre for Environmental Research - UFZ, Dept. of Computational Landscape Ecology
// Modified     : 
//		March  4, 2010:		modified methods for usage with bipartite networks
//		March 18, 2010:		functionality for exiting program after maxconverge time steps during which no better dendrogram was found
//					added input parameter handling for quantitative networks
//					modified functions for usage with quantitative networks
//		April 12, 2010:		enhanced program by holding a copy of current best dendrogram instead of writing to file immediately
//		April 13, 2010:		added further logging information
//		June 3-4, 2010:		cleaning up
//		February 13, 2011:	modification of output
//
// ****************************************************************************************************
// *** PROGRAM USAGE NOTES ****************************************************************************
// 
// The input to the algorithm must be a text file whose lines are formatted in the following way:
// a0	b0	w0
// a1	b1	w1
// a2	b2	w2
// and so on, where 1. ax, bx (0 <= x < |lines in file|) are non-negative integers with ax < by for all x,y
//		    2. wx is a positive float or double value representing the weight of the edge between ax and bx
// Each line is terminated by a carriage return and ax,bx and wx are separated by tabs.
// Multi-edges may appear, but will be stripped out automatically.
//
// For instance, here is a bipartite network with a0 connected to b0,b1,b2,b3 and a1 connected to b3,b4,b5, where all edge weights are 1.0 (this corresponds to a binary network):
//
// 0	2	1.0
// 0	3	1.0
// 0	4	1.0
// 0	5	1.0
// 1	5	1.0
// 1	6	1.0
// 1	7	1.0
//
// PLEASE NOTE that if the input .pairs file is formatted incorrectly, the program will crash.
//
// ****************************************************************************************************

#include <R.h>
//#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdio.h>
#include <fstream>
#include <string>
#include "stdlib.h"
#include "time.h"
#include "string.h"

#include "dendro.h"
#include "graph.h"
#include "rbtree.h"

using namespace std;

extern "C" { // wrapper for R

// ******** Function Prototypes ***************************************************************************

bool		markovChainMonteCarlo();
const char*	num2str(const unsigned int);
bool		parseCommandLine(int argc, char * argv[]);
bool		readInputFile();
void		recordModules();
void		recordNamesLUT();

// ******** Structures and Constants **********************************************************************

struct ioparameters {
	int		n_a;			// number of A vertices in input graph
	int		n_b;			// number of B vertices in input graph
	int		n;			// total number of vertices in input graph (i.e. n_a + n_b)
	int		m;			// number of edges in input graph
	unsigned int	maxconverge;		// maximum number of steps during which no better dendrogram was found
	double		temperature;		// SA temperature for computing modules
	double		tolerance;		// threshold for changes in likelihood
						// 	(new likelihood is considered "better" only if it is
						// 		greater than the old likelihood plus the threshold)
	string		d_dir;			// working directory
	string		f_in;			// name of input file (*.pairs)
	string		f_dg;			// name of output file
	string		f_dg_info;		// name of output information file
	string		f_ordA;			// name of output order file for A vertices
	string		f_ordB;			// name of output order file for B vertices
	string		f_modules;		// name of output file for modules
	string		f_pairs;		// name of output random graph file
	string		f_namesLUT;		// name of output names LUT file
	string		s_scratch;		// filename sans extension
	string		s_tag;			// user defined filename tag
	string		start_time;		// time simulation was started
	int		timer;			// timer for reading input
	bool		flag_filename;		// flag indicating whether -filename invoked
	bool		flag_steps;		// flag indicating whether -steps invoked
	bool		flag_onlyEdgeWeights;	// flag indicating whether -onlyEdgeWeights invoked
	bool		flag_tolerance;		// flag indicating whether -tolerance invoked
};

// ******** Global Variables ******************************************************************************

ioparameters	ioparm;				// program parameters
rbtree*		namesLUT;			// look-up table; translates input file vertex names to graph indices
rbtree*		reverseNamesLUT;		// reverse look-up table; translates graph indices to input file vertex names
dendro*		d;				// hrg data structure
dendro*		bestDendro;			// dendrogram with best modularity found so far
unsigned int	t;				// number of time steps
double		temperature;			// current SA temperature (default: 1e-1)
double		dTemperature;			// SA temperature step
double		averageStartTemperature;	// average temperature at which increase of
int		averageDenominator;		// auxiliary variable
unsigned int	converge;			// current number of steps without increase of modularity
double		minTemperature;			// minimum SA temperature (default: 1/100 of start temperature)
short int	billionCount;			// counts number of billion steps
double		bestM;				// best modularity found so far
int		nrOfRecordBreakings;		// counts number of improvements of modularity (after having reached a legal module dendrogram)
unsigned int	period;				// number of MCMC moves to do before writing stuff out (default: 10000)
double		currentM;			// holds current modularity
bool		firstLegalDendrogram;		// indicates whether the better dendrogram found is the first legal one
char*		method = new char[20];		// method to use (Newman's modularity or Strauss' hierarchical modularity function)
MTRand		mtr;				// Mersenne Twister random number generator instance
bool		graphIsConnected;		// indicates whether the graph is connected

// ******** Main Loop *************************************************************************************

void identifyModules(int* r_argc, char* argv[]) {
	
	namesLUT = new rbtree();
	reverseNamesLUT = new rbtree();
	ioparm.n_a			= 0;
	ioparm.n_b			= 0;
	ioparm.temperature		= 1e-5;
	ioparm.tolerance		= 1e-10;
	ioparm.timer   			= 1;
	ioparm.flag_filename		= false;
	ioparm.flag_steps		= false;
	ioparm.flag_onlyEdgeWeights	= false;
	ioparm.flag_tolerance		= false;
	ioparm.s_tag   			= "";
	ioparm.maxconverge		= 0;
	minTemperature			= 0;
	string input   			= "";
	t				= 1;
	billionCount			= 0;
	nrOfRecordBreakings		= 0;
	period				= 10000;
	firstLegalDendrogram		= true;
	
	strcpy(method, "Newman");
	
	time_t t1			= time(&t1);
	int argc			= *r_argc;

	if (parseCommandLine(argc, argv)) {
		d = new dendro(method);								// make the dendro-graph structure for computing
		ioparm.start_time = asctime(localtime(&t1));

		if (!readInputFile()) {								// read input file
			 Rprintf("!! ERROR: Malformed input file.\n");
			 //return 0;
		}

 		bestDendro		= d->deepCopy();					// make the dendro-graph structure holding a copy of the current dendrogram
 		bestM			= d->getModularity();					// store current modularity as best modularity
 
		temperature		= ioparm.temperature;						    // initialize temperature and 
		dTemperature	= (temperature - minTemperature) / (double)(ioparm.maxconverge);    // size of temperature steps

		Rprintf("identifyModules: start building legal dendrogram\n");

		if(!(d->g->nrOfComponents == min(ioparm.n_a, ioparm.n_b))) {
		    Rprintf("\n#steps\tM\tbest M\ttemperature\n");

		    while(converge < ioparm.maxconverge || bestM < 0) {				// while modularity did increase during the last maxconverge steps
                if (!(markovChainMonteCarlo())) {}// return 0; }			// make MCMC move
		    }
		}

		Rprintf("\nidentifyModules: finding best dendrogram complete\n\n");

		bestDendro->refreshModularity();

		if(!strcmp(method, "Newman")) {
			Rprintf("identifyModules: modularity = %g\n\n", bestDendro->getModularity() / 2);
		}
		else {
			Rprintf("identifyModules: modularity = %g\n\n", bestDendro->getModularity());
		}
		recordModules();								// write dendrogram to file

		delete d->g;									// delete graph (same instance for both d and bestDendro)
		delete d;
		//delete bestDendro->g;       // CFD trial
        delete bestDendro;
        namesLUT = 0;
        delete namesLUT;
        reverseNamesLUT = 0;
        delete reverseNamesLUT;
		
		//return 1;
	}

	 //return 0;
}

// ******** Function Definitions **************************************************************************

bool markovChainMonteCarlo() {

	double  dM;
	bool    flag_taken;
	time_t t1 = time(&t1);
	time_t t2 = time(&t2);

	// Because moves in the dendrogram space are chosen (Monte Carlo) so that we sample dendrograms 
	// with probability proportional to their modularity, a modularity-proportional sampling of 
	// the dendrogram models would be equivalent to a uniform sampling of the walk itself. We would
	// still have to decide how often to sample the walk (at most once every n steps is recommended)
	// but for simplicity, the code here simply runs the MCMC itself. To actually compute something
	// over the set of sampled dendrogram models (in a Bayesian model averaging sense), you'll need
	// to code that yourself.
	
	// do 'period' MCMC moves before doing anything else
	for (unsigned int i = 0; i <= period - 1; i++) {
		
		if (!(d->monteCarloMove(dM, flag_taken, temperature, bestM))) {				// make a MCMC move
			Rprintf("!! ERROR: failed to make monte carlo move");
			return false;
		}

		currentM = d->getModularity();								// get modularity of the current dendrogram

		if ((!ioparm.flag_tolerance && currentM > bestM) || (currentM > (bestM + ioparm.tolerance))) {	// in case that the modularity of the current dendrogram is higher than the best one found so far...
			d->refreshModularity();									// correct floating-point errors O(n)
			currentM = d->getModularity();								// get modularity of the current dendrogram
		}

		if ((!ioparm.flag_tolerance && currentM > bestM) || (currentM > (bestM + ioparm.tolerance))) {	// in case that the modularity of the current dendrogram is higher than the best one found so far...

			if(currentM >= 0 && bestM < 0) {
				if(!strcmp(method, "Newman")) {
					if(billionCount > 0) {
						Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM / 2, bestM / 2, temperature);
					}
					else {
						Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM / 2, bestM / 2, temperature);
					}
				}
				else {
					if(billionCount > 0) {
						Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM, bestM, temperature);
					}
					else {
						Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM, bestM, temperature);
					}
				}
			}

			bestM = currentM;								// ... store the current best modularity

			if(bestM >= 0) {
			    if(averageDenominator == 0) {							// ... update the SA temperature
				    averageStartTemperature = temperature;
				    averageDenominator	    = 1;
			    }
			    else {
				    averageStartTemperature *= averageDenominator;
				    averageStartTemperature += temperature;
				    averageDenominator++;
				    averageStartTemperature /= (double)(averageDenominator);
			    }

			    temperature	    = averageStartTemperature + mtr.randExc()*(ioparm.temperature - averageStartTemperature);
			    dTemperature    = (temperature - minTemperature) / (double)(ioparm.maxconverge);

    			    if(firstLegalDendrogram) {
					Rprintf("\nidentifyModules: building of legal dendrogram finished\n");
					Rprintf("identifyModules: start finding best dendrogram\n\n");
					if(!strcmp(method, "Newman")) {
						if(billionCount > 0) {
							Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM / 2, bestM / 2, temperature);
						}
						else {
							Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM / 2, bestM / 2, temperature);
						}
					}
					else {
						if(billionCount > 0) {
							Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM, bestM, temperature);
						}
						else {
							Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM, bestM, temperature);
						}
					}
					firstLegalDendrogram = false;
			    }

			    delete bestDendro;
			    bestDendro = d->deepCopy();							// ... copy d to bestDendro
			}

			nrOfRecordBreakings++;								// ... increment count of record-breakings

			converge = 0;									// ... and reset convergence counter
		}
		else {											// else...
			if(bestM >= 0 && temperature - dTemperature >= minTemperature) temperature -= dTemperature;	// ... decrease SA temperature

			converge++;									// ... and increase convergence counter
		}
	
		// check timer and write some stuff to standard-out to describe the current state of things if necessary
		t2 = time(&t2);
		if (t2 - t1 >= ioparm.timer || i == period - 1) {
			if(!strcmp(method, "Newman")) {
				if(billionCount > 0) {
					Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM / 2, bestM / 2, temperature);
				}
				else {
					Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM / 2, bestM / 2, temperature);
				}
			}
			else {
				if(billionCount > 0) {
					Rprintf("[%d%ld]\t%g\t\t(%g)\t\t%g\n", billionCount, t, currentM, bestM, temperature);
				}
				else {
					Rprintf("[%ld]\t%g\t\t(%g)\t\t%g\n", t, currentM, bestM, temperature);
				}
			}
			t1 = t2;
		}

		t++;
		if (t >= 1000000000) { billionCount++; t = 0; }						// rollover step count
	}

	d->refreshModularity();										// corrects floating-point errors O(n)

	return true;
}

// ********************************************************************************************************

const char* num2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	const char* str  = "";
	if (input == 0) {
		str = "0";
	}
	else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}

// ********************************************************************************************************

bool parseCommandLine(int argc, char * argv[]) {

	int argct = 1;
	string temp, ext;
	string::size_type pos;

	if (argc==1) {
		Rprintf("\n  -- Hierarchical Module Identification --\n");
		Rprintf("  by Rouven Strauss (copyright 2010-2011)\n\n");
		Rprintf("  based on the algorithm \n");
		Rprintf("\n  -- Hierarchical Random Graphs --\n");
		Rprintf("  by Aaron Clauset (copyright 2005-2009)\n\n");
		Rprintf("  Flags:\n");
		Rprintf("  -filename <file>               (required) input .pairs graph file\n");
		Rprintf("  -steps <integer>               (required) maximum number of steps without\n");
		Rprintf("                                            increase of best modularity\n");
		Rprintf("  -method <string>               (optional) method to use (Strauss or Newman)");
		Rprintf("  -label <string>                (optional) label of output of this run\n");
		Rprintf("  -temperature <float>	     (optional) SA start temperature [default: %g]\n", ioparm.temperature);
		Rprintf("  -tolerance <double>            (optional) tolerance for changes in likelihood\n");
		Rprintf("  -onlyEdgeWeights               (optional) use only edge weights\n\n");
		Rprintf("  examples:\n");
		Rprintf("  ./identifyModules -filename graph.pairs -steps 1000000\n");
		Rprintf("  ./identifyModules -filename graph.pairs -steps 1000000 -label test\n");
		Rprintf("  ./identifyModules -filename graph.pairs -steps 1000000 -tolerance 1e-14\n\n");
		return false;
		
	}
	else {
		// cout << "\n" << endl;
		while (argct < argc) {
			temp = argv[argct];
			
			if (temp == "-label") {
				argct++;
				ioparm.s_tag = argv[argct];
			}
			else if (temp == "-filename") {
				ioparm.flag_filename = true;
				argct++;
				temp = argv[argct];
				ext = ".pairs";
				pos = temp.find(ext,0);
				if (pos == string::npos) {
					// cout << "!! ERROR: Input file must claim to be .pairs format.\n";
					return false;
				}
				ioparm.f_in = temp;
				ext = "/";
				pos = string::npos;
				for (unsigned int i=0; i < temp.size(); i++) {
					if (temp[i] == '/') {
						pos = i;
					}
				}
				if (pos != string::npos) {
					ioparm.d_dir = temp.substr(0, pos+1);
					temp = temp.substr(pos+1,temp.size()-pos-1);
				}
				// now grab the filename sans extension for building outputs files
				for (unsigned int i=0; i < temp.size(); i++) {
					if (temp[i] == '.') {
						pos = i;
					}
				}
				ioparm.s_scratch = temp.substr(0,pos);
				
			}
			else if (temp == "-steps") {
				ioparm.flag_steps = true;
				argct++;
				if(atof(argv[argct]) < 0) {          // use atof here?? atoi converts to integers, atof to double
					Rprintf("!! ERROR: -steps argument has to be >= 0!\n");
					return false;
				}
				else {									// hand over an integer representing the maximum number
					ioparm.maxconverge = atof(argv[argct]);				// of steps without increase of best modularity before exiting
                    // use atof here?? atoi converts to integers, atof to double
				}
			}
			else if (temp == "-tolerance") {
				ioparm.flag_tolerance = true;
				argct++;
				if(atof(argv[argct]) < 0) {
					Rprintf("!! Error: -tolerance argument has to be >= 0!\n");
					return false;
				}
				else {
					ioparm.tolerance = atof(argv[argct]);
				}
			}
			else if (temp == "-method") {
				argct++;
				temp = argv[argct];
				if(!strcmp(temp.c_str(), "Newman")) {
					strcpy(method, "Newman");
				}
				else if(!strcmp(temp.c_str(), "Strauss")) {
					strcpy(method, "Strauss");
				}
				else if(strcmp(temp.c_str(), "Strauss")) {
				    Rprintf("!! ERROR: -method argument has to be 'Strauss' or 'Newman'\n");
				    return false;
				}
			}
			else if (temp == "-temperature") {
				argct++;
				if(atof(argv[argct]) <= 0) {
				    Rprintf("!! ERROR: -temperature argument has to be > 0!\n");
				    return false;
				}
				else {									// hand over a double representing the simulated annealing temperature
					ioparm.temperature = atof(argv[argct]);
				}
			}
			else if (temp == "-onlyEdgeWeights") {
				ioparm.flag_onlyEdgeWeights = true;
			}
			else if (temp == "-period") {
				argct++;
				if(atoi(argv[argct]) > 0) { period = atoi(argv[argct]); }
			}
			else {
				Rprintf("!! Warning: ignored argument nr. %d\n", argct);
			}
			argct++;
		}
	}

	if(!ioparm.flag_filename) {
		Rprintf("!! ERROR: flag -filename required!\n");
		return false;
	}

	if(!ioparm.flag_steps) {
		Rprintf("!! ERROR: -steps has to be invoked with appropriate parameters!\n");
		return false;
	}

	ioparm.f_namesLUT = ioparm.d_dir + ioparm.s_scratch + ".lut";
	
	if (ioparm.s_tag != "")    {
		ioparm.s_scratch = ioparm.s_tag;
	}

	if(ioparm.flag_onlyEdgeWeights) {
		Rprintf("identifyModules: only edge weights are being used (no expected edge weights).\n");
	}

	return true;
}

// ********************************************************************************************************

bool readInputFile() {

	FILE *fp;

	int n_a, n_b, n, m, vertex_i, vertex_j, virtualVertex_i, virtualVertex_j, countBVertices;
    double edgeWeight; //, sumEdgeWeight;
	n_a = n_b = n = m = countBVertices = 0;
	// sumEdgeWeight = 0;
	elementrb *item;
	time_t t1 = time(&t1);
	time_t t2 = time(&t2);

	// First, we scan through the input file to create a list of unique vertex names
	// (which we store in the namesLUT), and a count of the number of edges.
	// cout << "\n>> input file scan ( " << ioparm.f_in << " )" << endl;

	// cout << ">> starting to count edges" << endl;

	fp = fopen(ioparm.f_in.c_str(), "r");       				// check whether file exists at all
	if (fp != NULL) fclose(fp);									// exists: fine, now close it again
	else return false;										    // does not exist: return a "false"; does it need to be closed as well?

	ifstream fscan0(ioparm.f_in.c_str(), ios::in);							// read input

	if(fscan0.good()) {

		while (fscan0 >> vertex_i >> vertex_j >> edgeWeight) {					// read friendship pair (vertex_i,vertex_j)
			if (vertex_i == vertex_j) {
				// cout << "!! ERROR: corrupt input file! (same vertex number in one line)" << endl;
			}
			else {
				m++;									// count number of edges
				// sumEdgeWeight += edgeWeight;						// compute total sum of edge weights
				if (namesLUT->findItem(vertex_i) == NULL) {
					namesLUT->insertItem(vertex_i, n_a);
					n_a++;								// increment number of A vertices
				}
			}

    			t2=time(&t2);
			if (t2 - t1 >= ioparm.timer) {							// check timer and display if necessary
				// cout << ">> edges: ["<<m<<"]"<<endl;
				t1 = t2;
			}
		}
        }

	fscan0.close();

	countBVertices = n_a;
	
	ifstream fscan1(ioparm.f_in.c_str(), ios::in);

	if(fscan1.good()) {

		while (fscan1 >> vertex_i >> vertex_j >> edgeWeight) {
			if (vertex_i != vertex_j) {
				if (namesLUT->findItem(vertex_j) == NULL) {
					    namesLUT->insertItem(vertex_j, countBVertices);
					    countBVertices++;
					    n_b++;
				    }
			}
		}

	}

	fscan1.close();
	
	// cout << ">> total amount of edges: ["<<m<<"]"<<endl;
	// cout << ">> sum of edge weights: ["<<sumEdgeWeight<<"]"<<endl;

	d->g = new graph(n_a, n_b, method, ioparm.flag_onlyEdgeWeights);				// create new graph with (n_a + n_b) vertices

	// Finally, we reparse the file and add edges to the graph
	m			= 0;
	bool aToB		= true;

	ifstream fin(ioparm.f_in.c_str(), ios::in);

	if(fin.good()) {

	    while (fin >> vertex_i >> vertex_j >> edgeWeight) {
		    m++;
		    item = namesLUT->findItem(vertex_i); virtualVertex_i = item->value;
		    item = namesLUT->findItem(vertex_j); virtualVertex_j = item->value;
		    if (!(d->g->doesLinkExist(virtualVertex_i, virtualVertex_j))) {
			    if (!(d->g->addLink(virtualVertex_i, virtualVertex_j, edgeWeight, aToB))) {
				    // cout << "!! ERROR: couldn't insert edge (" << vertex_i << " " << vertex_j << " " << edgeWeight <<")" << endl;
				    return false;
			    }
		    }
		    if (!(d->g->doesLinkExist(virtualVertex_j, virtualVertex_i))) {
			    if (!(d->g->addLink(virtualVertex_j, virtualVertex_i, edgeWeight, !aToB))) {
				    // cout << "!! ERROR: couldn't insert edge (" << vertex_i << " " << vertex_j << " " << edgeWeight <<")" << endl;
				    return false;
			    }
		    }
	    }

	}

	fin.close();

	graphIsConnected = d->g->isConnected();
	if(!graphIsConnected) // cout << "\n>> graph is not connected\n" << endl;

	ioparm.m	= d->g->getNumLinks();								// store number of edges created
	ioparm.n_a	= d->g->getNumAVertices();							// store number of A vertices used
	ioparm.n_b	= d->g->getNumBVertices();							// store number of B vertices used
	ioparm.n	= d->g->getNumVertices();							// store total number of vertices used
	// cout << ">> number of A vertices: ["<<ioparm.n_a<<"]"<<endl;
	// cout << ">> number of B vertices: ["<<ioparm.n_b<<"]"<<endl;
	// cout << ">> total number of vertices: ["<<ioparm.n<<"]" << endl;
	
	recordNamesLUT();										// record names LUT to file for future reference

	if(!d->buildDendrogram()) return false;								// create dendrogram

	return true;
}

// ********************************************************************************************************

void recordModules() {

	FILE* infoFile;

	time_t t1;

	// write to files
	ioparm.f_dg = ioparm.d_dir + ioparm.s_scratch + ".den";
	ioparm.f_ordA = ioparm.d_dir + ioparm.s_scratch + ".ordA";
	ioparm.f_ordB = ioparm.d_dir + ioparm.s_scratch + ".ordB";
	ioparm.f_modules = ioparm.d_dir + ioparm.s_scratch + ".mod";

	if(!bestDendro->recordOrderAndModules(*reverseNamesLUT, ioparm.f_ordA, ioparm.f_ordB, ioparm.f_modules)) {
		// cout << "!! ERROR: failed to record order and/or module files" << endl;
        Rprintf("!! ERROR: failed to record order and/or module files");
		return;
	}

	bestDendro->recordDendrogramStructure(ioparm.f_dg);
	
	// write statistics about hrg to file
	ioparm.f_dg_info = ioparm.d_dir + ioparm.s_scratch + ".info";

	t1 = time(&t1); 
	infoFile = fopen(ioparm.f_dg_info.c_str(), "w"); 				//fopen_s removed: Windows thing (largely)
	fprintf(infoFile, 	"--- HIERARCHICAL MODULE IDENTIFICATION ---\n\n" \
				"StartTime                      : %s", ioparm.start_time.c_str());
	fprintf(infoFile, 	"EndTime                        : %s", asctime(localtime(&t1)));
	fprintf(infoFile, 	"InputFile                      : %s\n", ioparm.f_in.c_str());
	fprintf(infoFile, 	"Directory                      : %s\n", (ioparm.d_dir == "" ? "same as executable" : ioparm.d_dir.c_str()));
	fprintf(infoFile, 	"\n--- Information about input parameters ---\n\n" \
				"Method                         : %s\n", method);
	fprintf(infoFile, 	"Label                          : %s\n", (ioparm.s_tag != "" ? ioparm.s_tag.c_str() : "-"));
	fprintf(infoFile, 	"Maximal number of steps\n" \
				"without increase of modularity : %d\n", ioparm.maxconverge);
	fprintf(infoFile, 	"SA temperature                 : %f\n", ioparm.temperature);
	fprintf(infoFile, 	"Only edge weights              : %s\n", (ioparm.flag_onlyEdgeWeights ? "yes" : "no"));
	fprintf(infoFile, 	"\n--- Information about input graph ---\n\n" \
				"A vertices                     : %d\n", ioparm.n_a);
	fprintf(infoFile, 	"B vertices                     : %d\n", ioparm.n_b);
	fprintf(infoFile, 	"Total number of vertices       : %d\n", ioparm.n);
	fprintf(infoFile, 	"Number of edges                : %d\n", ioparm.m/2);
	fprintf(infoFile, 	"\n--- Information about hierarchical module dendrogram ---\n\n");
	if(!strcmp(method, "Newman")) {
		fprintf(infoFile, 	"Modularity                     : %f\n", bestM / 2);
	}
	else {
		fprintf(infoFile, 	"Modularity                     : %f\n", bestM);
	}
	if(billionCount > 0) fprintf(infoFile, 	"Number of MCMC steps           : %d %s %d\n", billionCount, (billionCount > 1 ? "billions" : "billion"), t);
	else fprintf(infoFile, 	"Number of MCMC steps           : %d\n", t);
	fprintf(infoFile, 	"Number of improvements         : %d\n", nrOfRecordBreakings);
	fprintf(infoFile, 	"\n--- Information about created files ---\n\n" \
				"Created files                  : %s.mod\n", ioparm.s_scratch.c_str());
	fprintf(infoFile, 	"                                 %s.ordA\n", ioparm.s_scratch.c_str());
	fprintf(infoFile, 	"                                 %s.ordB\n", ioparm.s_scratch.c_str());
	fprintf(infoFile, 	"                                 %s.den\n", ioparm.s_scratch.c_str());
	fprintf(infoFile, 	"                                 %s-names.lut\n", ioparm.s_scratch.c_str());
	fprintf(infoFile, 	"                                 %s.info\n", ioparm.s_scratch.c_str());
	
	fclose(infoFile);
	
	// cout << ">> recorded general information into file: " << ioparm.f_dg_info << endl;
	
	return;
}

// ********************************************************************************************************

void recordNamesLUT() {

	FILE* lutFile;
	
	keyValuePair *head, *prev;
	
	head = namesLUT->returnTreeAsList();
	while (head != NULL) {
		reverseNamesLUT->insertItem(head->y, head->x);
		prev = head;
		head = head->next;
		delete prev;
	}
	head = NULL; prev = NULL;
	
	elementrb *item;
	
	lutFile = fopen(ioparm.f_namesLUT.c_str(), "w"); //fopen_s
	
	fprintf(lutFile, "virtual\treal\n");
	for (int i=0; i<ioparm.n; i++) {
		item = reverseNamesLUT->findItem(i);
		fprintf(lutFile, "%d\t%d\n", i, item->value);
	}
	
	fclose(lutFile);
	
	// cout << ">> recorded names look-up table to file: " << ioparm.f_namesLUT << endl;
	
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

} // end wrapper for R
