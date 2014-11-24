// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph.h - graph data structure for hierarchical random graphs
// Copyright (C) 2005-2008 Aaron Clauset
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
// Created      : fitHRG:	8 November 2005
//		: consensusHRG:	8 November 2005
//		: predictHRG:	21 June 2006
//
// Modified     : fitHRG:	23 December 2007 (cleaned up for public consumption)
//		: consensusHRG:	23 December 2007 (cleaned up for public consumption)
//		: predictHRG:	23 December 2007 (cleaned up for public consumption)
//
// Modified by Rouven Strauss:
//		March 13, 2010:		merged graph.h files of fitHRG, consensusHRG and predictHRG into one
//					merged graph_simp.h into this one
//		March 18, 2010:		modified functions for usage with quantitative networks
//		March 21, 2010:		added datastructure for computation of expected value
//
// ****************************************************************************************************
// 
// Graph data structure for maximum likelihood hrgs. The basic structure is an adjacency list of
// edges; however, many additional pieces of metadata are stored as well. Each vertex stores its
// external name and its degree. Each edge stores a histogram of the probabilities assigned to it 
// by the dendrogram structure. Generally, edges are directional, an adjacency (i,j) should also be 
// stored as (j,i) in order to maintain the undirectedness of the graph.
// 
// ****************************************************************************************************

#if !defined(graph_INCLUDED)
#define graph_INCLUDED

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include "math.h"

#include "rbtree.h"

using namespace std;

// ******** Basic Structures ******************************************************************************

#if !defined(edge_INCLUDED)
#define edge_INCLUDED
class edge {
public:

	edge();
	~edge();
	
	int	x;					// index of edge terminator
	double	weight;					// weight of edge
	double	originalWeight;				//
	edge*	next;					// pointer to next elementd
};

edge::edge()  {
	x		= -1;
	weight		= 1;
	originalWeight	= 1;
	next		= NULL;
}

edge::~edge() {
    next = NULL;
}

#endif

#if !defined(vert_INCLUDED)
#define vert_INCLUDED
class vert {
public:
	vert();
	~vert();

	string	name;					// (external) name of vertex
	int	degree;					// degree of this vertex

};

vert::vert()  {
	name = "";
	degree = 0;

}

vert::~vert() {}
#endif

#if !defined(block_INCLUDED)
#define block_INCLUDED
struct block { double x; int y; };
#endif

// ******** Graph Class with Edge Statistics *************************************************************

class graph {
public:
	int	nrOfComponents;								// number of disjoint subgraphs
	block*	componentNr;								// array indicating for each vertex to which disjoint subgraph it belongs

	graph(const int, const int, const char*, const bool);
	~graph();

	bool	addLink(const int, const int, const double, const bool);		// add edge (i,j) with weight to graph
	bool	doesLinkExist(const int, const int);					// true if edge (i,j) is already in graph
	bool	isConnected();								// checks whether graph is connected
	void	visit(int, int);							// breadth-first visit of vertices
	double	getOriginalEdgeWeight(const int, const int);				// returns weight of edge (i,j)
	double	getExpectedEdgeWeight(const int, const int);				// returns expected value of edge weight between i and j;
	void	updateEdgeWeights();							// sets edge weights according to method "Strauss"
	edge*	getNeighborList(const int);						// returns edge list of vertex i
	int	getNumLinks();								// returns m
	int	getNumAVertices();							// returns n_a
	int	getNumBVertices();							// returns n_b
	int	getNumVertices();							// returns n
	double	getSumEdgeWeight();							// returns sumEdgeWeight
	double	getMarginTotal(const int);						// returns marginTotal[i]
	void	printPairs();								// prints all edges in graph
	
private:
	edge**		vertexLink;		// linked list of neighbors to vertex
	edge**		vertexLinkTail;		// pointers to tail of neighbor list
	int		n_a;			// number of A vertices
	int		n_b;			// number of B vertices
	int		n;			// number of vertices
	const char*	method;			// method ("Strauss" or "Newman")
	double		sumEdgeWeight;		// total sum of edge weights
	int		m;			// number of directed edges
	double*		marginTotal;		// contains the margin totals needed for computing expectedEdgeWeight
	bool		onlyEdgeWeights;	// flag indicating whether only edge weights should be used (no expected edge weights)
};

// ******** Constructor / Destructor **********************************************************************

graph::graph(const int sizeOfA, const int sizeOfB, const char* usedMethod, const bool flag_onlyEdgeWeights)  {
	nrOfComponents		= 1;
	n_a			= sizeOfA;
	n_b			= sizeOfB;
	n			= sizeOfA + sizeOfB;
	method			= usedMethod;
	sumEdgeWeight		= 0;
	m			= 0;
	onlyEdgeWeights		= flag_onlyEdgeWeights;
	vertexLink		= new edge*[n];
	vertexLinkTail		= new edge*[n];
	marginTotal		= new double[n];
	componentNr		= new block[n];

	for(int i = 0; i < n; i++) {
		vertexLink[i]	    = NULL;
		vertexLinkTail[i]   = NULL;
		marginTotal[i]	    = 0;
		componentNr[i].x    = -1;
		componentNr[i].y    = i;
	}
}

graph::~graph() {

	edge* currentEdge;
	edge* toDelete;

	for(int i = 0; i < n; i++) {
	    currentEdge = vertexLink[i];

	    while(currentEdge != NULL) {
		toDelete	= currentEdge;
		currentEdge	= currentEdge->next;
		delete toDelete;
	    }
	}

	currentEdge = NULL;
	toDelete    = NULL;

	delete [] vertexLink;	    vertexLink	    = NULL;
	delete [] vertexLinkTail;   vertexLinkTail  = NULL;
	delete [] marginTotal;	    marginTotal	    = NULL;
	delete [] componentNr;	    componentNr	    = NULL;
}

// ********************************************************************************************************

bool graph::addLink(const int i, const int j, const double weight, const bool aToB) {	// adds the directed edge (i,j,weight) to the adjacency list for v_i
	if (i >= 0 && i < n && j >= 0 && j < n && ((i < n_a && j >= n_a) || (j < n_a && i >= n_a))) {

		edge* newedge	        = new edge;
		newedge->x		= j;
		newedge->weight		= weight;
		newedge->originalWeight	= weight;

		if(aToB) {
			if (i < n_a && j >= n_a) {
				sumEdgeWeight	+= weight;
				marginTotal[i]	+= weight;

				if(i != j) {
					marginTotal[j]	+= weight;
				}
			}
			else return false;
		}

		if (vertexLink[i] == NULL) {				// first neighbor
			vertexLink[i]	    = newedge;
			vertexLinkTail[i]   = newedge;
		}
		else {							// subsequent neighbor
			vertexLinkTail[i]->next	= newedge;
			vertexLinkTail[i]	= newedge;
		}

		m++;							// increment edge count
		newedge = NULL;

		return true;
	}

	return false;
}

// ********************************************************************************************************

bool graph::doesLinkExist(const int i, const int j) {	// determines if the edge (i,j) already exists in the adjacency list of v_i
	if (i >= 0 && i < n && j >= 0 && j < n && ((i < n_a && j >= n_a) || (j < n_a && i >= n_a))) {

		edge* curr = vertexLink[i];

		while (curr != NULL) {
			if (curr->x == j) {
				curr = NULL;
				return true;
			}
			curr = curr->next;
		}
		curr = NULL;
	}
	return false;
}

// ********************************************************************************************************

bool graph::isConnected() {

    visit(0, nrOfComponents);

    for(int i = 0; i < n; i++) {
	if(componentNr[i].x == -1) {
	    nrOfComponents++;
	    visit(i, nrOfComponents);
	}
    }

    if(nrOfComponents == 1) return true;
    else return false;
}

// ********************************************************************************************************

void graph::visit(int v, int component) {

    componentNr[v].x = component-1;

    edge* curr = vertexLink[v];

    while(curr != NULL) {
	if(componentNr[curr->x].x == -1) {
	    visit(curr->x, component);
	}
        curr = curr->next;
    }

    curr = NULL;
}

// ********************************************************************************************************

double graph::getOriginalEdgeWeight(const int i, const int j) {	// returns the weight of edge (i,j)
	if (i >= 0 && i < n && j >= 0 && j < n && ((i < n_a && j >= n_a) || (j < n_a && i >= n_a))) {

		double weight;
		edge* curr = vertexLink[i];

		while (curr != NULL) {
			if (curr->x == j) {
				weight = curr->originalWeight;
				curr = NULL;
				return weight;
			}
			curr = curr->next;
		}
		curr = NULL;
	}
	return 0;
}

// ********************************************************************************************************

double	graph::getExpectedEdgeWeight(const int i, const int j) {

	if (0 <= i && i < n && 0 <= j && j < n) {

		if(onlyEdgeWeights) return 0;

		if ((i < n_a && j >= n_a) || (j < n_a && i >= n_a)) {
			if(!strcmp(method, "Strauss")) return marginTotal[i] * marginTotal[j] / (double)(sumEdgeWeight) / (getMarginTotal(i) + getMarginTotal(j) - getOriginalEdgeWeight(i,j));
			else return marginTotal[i] * marginTotal[j] / (double)(sumEdgeWeight);
		}
		else {
			return 0.0;
		}
	}
	else {
		//cout << "!! ERROR: trying to get expected edge weight between vertices of which at least one does not exist" << endl;
		return 0.0;
	}
}

// ********************************************************************************************************

void graph::updateEdgeWeights() {
    
    for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++) {

	    edge* curr = vertexLink[i];

	    while (curr != NULL) {
		if (curr->x == j) {
		    if(i != j) curr->weight = curr->weight / (getMarginTotal(i) + getMarginTotal(curr->x) - curr->weight);
		    else curr->weight = curr->weight / getMarginTotal(i);
		    break;
		}
		curr = curr->next;
	    }
	    curr = NULL;
	}
    }

    return;
}

// ********************************************************************************************************

// NOTE: The following method returns addresses; deallocation of returned object is dangerous
edge* graph::getNeighborList(const int i) {
	if (i >= 0 && i < n) {
		return vertexLink[i];
	}
	else {
		return NULL;
	}
}

// ********************************************************************************************************

int	graph::getNumLinks()		{ return m; }
int	graph::getNumAVertices()	{ return n_a; }
int	graph::getNumBVertices()	{ return n_b; }
int	graph::getNumVertices()		{ return n; }
double	graph::getSumEdgeWeight()	{ return sumEdgeWeight; }

// ********************************************************************************************************

double graph::getMarginTotal(int i) {
    if(0 <= i && i < n) return marginTotal[i];
    else return -1.0;
}

// ********************************************************************************************************

void graph::printPairs() {
	edge* curr;
	int edgeCount = 0;
	for (int i=0; i<n; i++) {
		//cout << "[" << i << "]\t";
		curr = vertexLink[i];
		while (curr != NULL) {
			//cout << curr->x << "\t";
			edgeCount++;
			curr = curr->next;
		}
		//cout << "\n";
	}
	//cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

#endif
