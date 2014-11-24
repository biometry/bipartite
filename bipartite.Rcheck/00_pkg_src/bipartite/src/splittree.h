// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// consensus_splittree.h - a binary search tree data structure for storing dendrogram split frequencies
// Copyright (C) 2006-2008 Aaron Clauset
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
// Author	: Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project	: Hierarchical Random Graphs
// Location	: University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created	: 19 April 2006
// Modified	: 19 May 2007
//		: 20 May 2008 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// Data structure for storing the split frequences in the sampled dendrograms. Data is stored 
// efficiently as a red-black binary search tree (this is a modified version of the rbtree.h file).
//
// ****************************************************************************************************

#if !defined(splittree_INCLUDED)
#define splittree_INCLUDED

#include <fstream>
#include <iostream>

using namespace std;

// ******** Basic Structures ******************************************************************************

#if !defined(slist_INCLUDED)
#define slist_INCLUDED
class slist {
public:
	string	x;				// stored elementd in linked-list
	slist*	next;				// pointer to next elementd
	slist(); ~slist();
};
slist::slist()  { x = ""; next = NULL; }
slist::~slist() {}
#endif

class keyValuePairSplit {
public:
	string	x;				// elementsp split (string)
	double	y;				// stored weight   (double)
	int	c;				// stored count    (int)
	keyValuePairSplit* next;		// linked-list pointer
	keyValuePairSplit(); ~keyValuePairSplit();
};
keyValuePairSplit::keyValuePairSplit()  { x = ""; y = 0.0; c = 0; next = NULL; }
keyValuePairSplit::~keyValuePairSplit() {}

// ******** Tree elementsp Class ****************************************************************************

class elementsp {
public:
	string	split;				// split represented as a string
	double	weight;				// total weight of this split
	int	count;				// number of observations of this split
	
	bool	color;				// F: BLACK
						// T: RED
	short int mark;				// marker

	elementsp   *parent;			// pointer to parent node
	elementsp   *left;			// pointer for left subtree
	elementsp   *right;			// pointer for right subtree
	
	elementsp(); ~elementsp();
};
elementsp::elementsp()  { split = ""; weight = 0.0; count = 0; color = false; mark = 0; parent  = NULL; left  = NULL; right  = NULL; }
elementsp::~elementsp() {}

// ******** Red-Black Tree Class **************************************************************************
/*   This vector implementation is a red-black balanced binary tree data structure.
 *   It provides find a stored elementsp in time O(log n), find the maximum elementsp in time O(1),
 *	delete an elementsp in time O(log n), and insert an elementsp in time O(log n).
 *
 *	Note that the split="" is assumed to be a special value, and thus you cannot insert such an item. 
 *	Beware of this limitation.
 */

class splittree {
private:
	elementsp*		root;					// binary tree root
	elementsp*		leaf;					// all leaf nodes
	int			support;				// number of nodes in the tree
	double			total_weight;				// total weight stored
	int			total_count;				// total number of observations stored

	void			rotateLeft(elementsp*);			// left-rotation operator
	void			rotateRight(elementsp*);		// right-rotation operator
	void			insertCleanup(elementsp*);		// house-keeping after insertion
	void			deleteCleanup(elementsp*);				// house-keeping after deletion
	keyValuePairSplit*	returnSubtreeAsList(elementsp*, keyValuePairSplit*);
	void			printSubTree(elementsp*);				// display the subtree rooted at z
	void			deleteSubTree(elementsp*);		// delete subtree rooted at z
	elementsp*		returnMinKey(elementsp*);		// returns minimum of subtree rooted at z
	elementsp*		returnSuccessor(elementsp*);		// returns successor of z's key
	
public:
	splittree(); ~splittree();					// default constructor/destructor

	double			returnValue(const string);		// returns value associated with searchKey
	elementsp*		findItem(const string);			// returns T if searchKey found, and
									// points foundNode at the corresponding node
	void			finishedThisRound();			// update total_count and total_weight
	bool			insertItem(string, double);		// insert a new key with stored value
	void			clearTree();
	void			deleteItem(string);			// delete a node with given key
	void			deleteTree();				// delete the entire tree
	string*			returnArrayOfKeys();			// return array of keys in tree
	slist*			returnListOfKeys();			// return list of keys in tree
	keyValuePairSplit*	returnTreeAsList();			// return the tree as a list of keyValuePairSplits
	keyValuePairSplit	returnMaxKey();				// returns the maximum key in the tree
	keyValuePairSplit	returnMinKey();				// returns the minimum key in the tree
	int			returnNodecount();			// returns number of items in tree
	keyValuePairSplit*	returnTheseSplits(const int);			// returns list of splits with given number of Ms
	double			returnTotal();					// returns sum of stored values

	void			printTree();					// displays tree (in-order traversal)
	void			printTreeAsList();				// list keys (in-order) with values
	void			printTreeAsShortList();				// short list keys (in-order) with values,
	void			recordTreeAsList(const string, const double);	// write list of keys and values to file
};

// ******** Red-Black Tree Methods ************************************************************************

splittree::splittree() {
	root = new elementsp;
	leaf = new elementsp;

	leaf->parent   = root;

	root->left	= leaf;
	root->right    = leaf;
	support		= 0;
	total_weight	= 0.0;
	total_count	= 0;
}

splittree::~splittree() {
	if (root != NULL && (root->left != leaf || root->right != leaf)) { deleteSubTree(root); }
	support      = 0;
	total_weight = 0.0;
	total_count  = 0;
	delete leaf;
	root		   = NULL;
	leaf		   = NULL;
}

void splittree::deleteTree() { if (root != NULL) { deleteSubTree(root); } return; }

void splittree::deleteSubTree(elementsp *z) {

	if (z->left  != leaf) { deleteSubTree(z->left);  }
	if (z->right != leaf) { deleteSubTree(z->right); }
	delete z;
	z = NULL;
	return;
}

// ******** Reset Functions *******************************************************************************

void splittree::clearTree() { // O(n lg n)
	string *array = returnArrayOfKeys();
	for (int i=0; i<support; i++) { deleteItem(array[i]); }
	delete [] array;
	return;
}

// ******** Search Functions ******************************************************************************

// public search function - if there exists a elementsp in the tree with key=searchKey,
// it returns TRUE and foundNode is set to point to the found node; otherwise, it sets
// foundNode=NULL and returns FALSE
elementsp* splittree::findItem(const string searchKey) {

	elementsp *current;    current = root;
	if (current->split=="") { return NULL; }						// empty tree; bail out
	while (current != leaf) {
		if (searchKey < current->split) {						// left-or-right?
			if (current->left  != leaf) { current = current->left;  }		// try moving down-left
			else { return NULL; }							//   failure; bail out
		} else {									// 
			if (searchKey > current->split) {					// left-or-right?
				if (current->right  != leaf) { current = current->right;  }	// try moving down-left
				else { return NULL; }						//   failure; bail out
			} else { return current; }						// found (searchKey==current->split)
		}
	}
	return NULL;
}

double splittree::returnValue(const string searchKey) {
	elementsp* test = findItem(searchKey);
	if (test == NULL) { return 0.0; } else { return test->weight; }
}


// ******** Return Item Functions *************************************************************************
// public function which returns the tree, via pre-order traversal, as a linked list

string* splittree::returnArrayOfKeys() {
	string* array;
	array = new string [support];
	bool flag_go = true;
	int index = 0;
	elementsp *curr;
	
	if (support == 1) { array[0] = root->split; }
	else if (support == 2) {
		array[0] = root->split;
		if (root->left == leaf) { array[1] = root->right->split; } 
		else { array[1] = root->left->split; }
	} else {
		for (int i=0; i<support; i++) { array[i] = -1; }
		// non-recursive traversal of tree structure
		curr		 = root;
		curr->mark = 1;
		while (flag_go) {
			
			if (curr->mark == 1 and curr->left == leaf) {		// - is it time, and is left child the leaf node?
				curr->mark = 2;					// 
			}
			if (curr->mark == 2 and curr->right == leaf) {		// - is it time, and is right child the leaf node?
				curr->mark = 3;					// 
			}
			if (curr->mark == 1) {					// - go left
				curr->mark = 2;					// 
				curr       = curr->left;			// 
				curr->mark = 1;					// 
			} else if (curr->mark == 2) {				// - else go right
				curr->mark = 3;					// 
				curr       = curr->right;			// 
				curr->mark = 1;					// 
			} else {						// - else go up a level
				curr->mark = 0;					// 
				array[index++] = curr->split;			// 
				curr = curr->parent;				// 
				if (curr == NULL) { flag_go = false; }		// 
			}
		}
	}
	
	return array;
} // This does not leak memory (unlike returnListOfKeys)

slist* splittree::returnListOfKeys() {
	keyValuePairSplit	*curr, *prev;
	slist			*head, *tail, *newlist;

	curr = returnTreeAsList();
	while (curr != NULL) {
		newlist    = new slist;
		newlist->x = curr->x;
		if (head == NULL) { head       = newlist; tail = head;    }
		else              { tail->next = newlist; tail = newlist; }
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	return head;
}

keyValuePairSplit* splittree::returnTreeAsList() { // pre-order traversal
	keyValuePairSplit  *head, *tail;

	head    = new keyValuePairSplit;
	head->x = root->split;
	head->y = root->weight;
	head->c = root->count;
	tail    = head;

	if (root->left  != leaf) { tail = returnSubtreeAsList(root->left,  tail); }
	if (root->right != leaf) { tail = returnSubtreeAsList(root->right, tail); }
	
	if (head->x == "") { return NULL; /* empty tree */ } else { return head; }
}

keyValuePairSplit* splittree::returnSubtreeAsList(elementsp *z, keyValuePairSplit *head) {
	keyValuePairSplit *newnode, *tail;
	
	newnode    = new keyValuePairSplit;
	newnode->x = z->split;
	newnode->y = z->weight;
	newnode->c = z->count;
	head->next = newnode;
	tail       = newnode;
	
	if (z->left  != leaf) { tail = returnSubtreeAsList(z->left,  tail); }
	if (z->right != leaf) { tail = returnSubtreeAsList(z->right, tail); }
	
	return tail;
}

keyValuePairSplit splittree::returnMaxKey() {
	keyValuePairSplit themax;
	elementsp *current;
	current = root;
	while (current->right != leaf) {		// search to bottom-right corner of tree
		current = current->right; }		// 
	themax.x = current->split;			// store the data found
	themax.y = current->weight;			// 
	
	return themax;						// return that data
}

keyValuePairSplit splittree::returnMinKey() {
	keyValuePairSplit themin;
	elementsp *current;
	current = root;
	while (current->left != leaf) {		// search to bottom-left corner of tree
		current = current->left; }		// 
	themin.x = current->split;			// store the data found
	themin.y = current->weight;			// 
	
	return themin;						// return that data
}

// private functions for deleteItem() (although these could easily be made public, I suppose)
elementsp* splittree::returnMinKey(elementsp *z) {
	elementsp *current;

	current = z;
	while (current->left != leaf) {		// search to bottom-right corner of tree
		current = current->left; }		// 
	return current;					// return pointer to the minimum
}

elementsp* splittree::returnSuccessor(elementsp *z) {
	elementsp *current, *w;
	
	w = z;
	if (w->right != leaf) {				// if right-subtree exists, return min of it
		return returnMinKey(w->right); }
	current = w->parent;				// else search up in tree
	while ((current!=NULL) && (w==current->right)) {
		w       = current;
		current = current->parent;		// move up in tree until find a non-right-child
	}
	return current;
}

int splittree::returnNodecount() { return support; }

keyValuePairSplit* splittree::returnTheseSplits(const int target) {
	keyValuePairSplit *head, *curr, *prev, *newhead, *newtail, *newpair;
	int count, len;
	
	head = returnTreeAsList();
	prev = newhead = newtail = newpair = NULL;
	curr = head;
	
	while (curr != NULL) {
		count = 0;
		len   = curr->x.size();
		for (int i=0; i<len; i++) { if (curr->x[i] == 'M') { count++; } }
		if (count == target and curr->x[1] != '*') {
			newpair       = new keyValuePairSplit;
			newpair->x    = curr->x;
			newpair->y    = curr->y;
			newpair->next = NULL;
			if (newhead == NULL) { newhead = newpair; newtail = newpair; }
			else { newtail->next = newpair; newtail = newpair; }
		}
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	
	return newhead;
}

double splittree::returnTotal() { return total_weight; }

// ******** Insert Functions ******************************************************************************

void splittree::finishedThisRound() {
	// We need to also keep a running total of how much weight has been added to the histogram.
	if (total_count == 0) { total_weight  = 1.0; total_count = 1; }
	else				  { total_weight += 1.0; total_count++; }
	return;
}	

// public insert function
bool splittree::insertItem(string newKey, double newValue) {
	
	// first we check to see if newKey is already present in the tree; if so, we do nothing;
	// if not, we must find where to insert the key
	elementsp *newNode, *current;
	
	current = findItem(newKey);						// find newKey in tree; return pointer to it O(log k)
	if (current != NULL) {
		// Add weight to the existing item's weight
		if (current->weight > total_weight) {
//			cout << "ERROR: " << current->weight << " > " << total_weight << endl;
//          cout << "for split " << current->split << endl;
			return false;
		}
		current->weight += 1.0;
		// And finally, we keep track of how many observations went into the histogram
		current->count++;
		return true;
	} else {
		newNode			= new elementsp;			// elementsp for the splittree
		newNode->split		= newKey;				//  store newKey
		newNode->weight		= newValue;  				//  store newValue
		newNode->color		= true;					//  new nodes are always RED
		newNode->parent		= NULL;					//  new node initially has no parent
		newNode->left		= leaf;					//  left leaf
		newNode->right		= leaf;					//  right leaf
		newNode->count		= 1;					// 
		support++;							// increment node count in splittree
		
		// must now search for where to insert newNode, i.e., find the correct parent and
		// set the parent and child to point to each other properly
		current = root;
		if (current->split=="") {									// insert as root
			delete root;										//   delete old root
			root		= newNode;								//   set root to newNode
			leaf->parent	= newNode;								//   set leaf's parent
			current		= leaf;									//   skip next loop
		}
		
		while (current != leaf) {									// search for insertion point
			if (newKey < current->split) {								// left-or-right?
				if (current->left != leaf) { current = current->left; }				// try moving down-left
				else {										// else found new parent
					newNode->parent		= current;					// set parent
					current->left		= newNode;					// set child
					current			= leaf;						// exit search
				}
			} else {										// 
				if (current->right != leaf) { current = current->right; }   			// try moving down-right
				else {										// else found new parent
					newNode->parent		= current;					//    set parent
					current->right		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			}
		}

		// now do the house-keeping necessary to preserve the red-black properties
		insertCleanup(newNode);										// do house-keeping to maintain balance
		
	}
	return true;
}

// private house-keeping function for insertion
void splittree::insertCleanup(elementsp *z) {
	
	if (z->parent==NULL) {										// fix now if z is root
		z->color = false; return; }
	elementsp *temp;
	while (z->parent!=NULL && z->parent->color) {							// while z is not root and z's parent is RED
		if (z->parent == z->parent->parent->left) {  						// z's parent is LEFT-CHILD
			temp = z->parent->parent->right;						// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  				// color z's parent BLACK	(Case 1)
				temp->color			= false;  				// color z's uncle BLACK	(Case 1)
				z->parent->parent->color = true;   					// color z's grandparent RED	(Case 1)
				z = z->parent->parent;							// set z = z's grandparent	(Case 1)
			} else {
				if (z == z->parent->right) {						// z is RIGHT-CHILD
					z = z->parent;							// set z = z's parent		(Case 2)
					rotateLeft(z);							// perform left-rotation	(Case 2)
				}
				z->parent->color = false;  						// color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   					// color z's grandparent RED	(Case 3)
				rotateRight(z->parent->parent);    					// perform right-rotation	(Case 3)
			}
		} else {										// z's parent is RIGHT-CHILD
			temp = z->parent->parent->left;							// grab z's uncle
			if (temp->color) {
				z->parent->color = false;  						// color z's parent BLACK	(Case 1)
				temp->color = false;							// color z's uncle BLACK	(Case 1)
				z->parent->parent->color = true;   					// color z's grandparent RED	(Case 1)
				z = z->parent->parent;							// set z = z's grandparent	(Case 1)
			} else {
				if (z == z->parent->left) {						// z is LEFT-CHILD
					z = z->parent;								// set z = z's parent	(Case 2)
					rotateRight(z);							// perform right-rotation	(Case 2)
				}
				z->parent->color = false;  						// color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   					// color z's grandparent RED	(Case 3)
				rotateLeft(z->parent->parent);						// perform left-rotation	(Case 3)
			}
		}
	}

	root->color = false;										// color the root BLACK
	return;
}

// ******** Delete Functions ******************************************************************************
// public delete function
void splittree::deleteItem(string killKey) {
	elementsp *x, *y, *z;
	
	z = findItem(killKey);
	if (z == NULL) { return; }						// item not present; bail out

	if (support==1) {							// -- attempt to delete the root
		root->split	= "";						// restore root node to default state
		root->weight    = 0.0;						// 
		root->color    = false;						// 
		root->parent   = NULL;						// 
		root->left	= leaf;						// 
		root->right    = leaf;						// 
		support--;								// set support to zero
		total_weight   = 0.0;						// set total weight to zero
		total_count--;								// 
		return;									// exit - no more work to do
	}
	
	if (z != NULL) {
		support--;								// decrement node count
		if ((z->left == leaf) || (z->right==leaf)) {		// case of less than two children
			  y = z; }							//    set y to be z
		else { y = returnSuccessor(z); }				//    set y to be z's key-successor
		
		if (y->left!=leaf) { x = y->left; }			// pick y's one child (left-child)
		else			    { x = y->right; }			//				  (right-child)
		x->parent = y->parent;						// make y's child's parent be y's parent

		if (y->parent==NULL) { root = x; }				// if y is the root, x is now root
		else {									// 
			if (y == y->parent->left) {				// decide y's relationship with y's parent
				y->parent->left  = x;				//   replace x as y's parent's left child
			} else {								// 
				y->parent->right = x; }				//   replace x as y's parent's left child
		}										// 

		if (y!=z) {								// insert y into z's spot
			z->split		= y->split;				// copy y data into z
			z->weight		= y->weight;				// 
			z->count		= y->count;				// 
		}										// 

		if (y->color==false) { deleteCleanup(x); }		// do house-keeping to maintain balance
		delete y;									// deallocate y
		y = NULL;									// point y to NULL for safety
	}											// 
		
	return;
}

void splittree::deleteCleanup(elementsp *x) {
	elementsp *w, *t;
	while ((x != root) && (x->color==false)) {			// until x is the root, or x is RED
		if (x==x->parent->left) {					// branch on x being a LEFT-CHILD
			w = x->parent->right;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color = false;					// color w BLACK				(case 1)
				x->parent->color = true;				// color x's parent RED			(case 1)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 1)
				w = x->parent->right;				// make w be x's right sibling	(case 1)
			}
			if ((w->left->color==false) && (w->right->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x = x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->right->color==false) {			// 
					w->left->color = false;			// color w's left child BLACK		(case 3)
					w->color = true;				// color w RED					(case 3)
					t = x->parent;					// store x's parent
					rotateRight(w);				// right rotation on w			(case 3)
					x->parent = t;					// restore x's parent
					w = x->parent->right;			// make w be x's right sibling	(case 3)
				}								// 
				w->color			= x->parent->color; // make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->right->color	= false;			// color w's right child BLACK	(case 4)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 4)
				x = root;							// finished work. bail out		(case 4)
			}									// 
		} else {									// x is RIGHT-CHILD
			w = x->parent->left;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color			= false;			// color w BLACK				(case 1)
				x->parent->color    = true;			// color x's parent RED			(case 1)
				rotateRight(x->parent);				// right rotation on x's parent	(case 1)
				w				= x->parent->left;  // make w be x's left sibling		(case 1)
			}
			if ((w->right->color==false) && (w->left->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x= x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->left->color==false) {			// 
					w->right->color	= false;		// color w's right child BLACK	(case 3)
					w->color			= true;		// color w RED					(case 3)
					t				= x->parent;   // store x's parent
					rotateLeft(w);					// left rotation on w			(case 3)
					x->parent			= t;			// restore x's parent
					w = x->parent->left;			// make w be x's left sibling		(case 3)
				}								// 
				w->color = x->parent->color;			// make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->left->color		= false;			// color w's left child BLACK		(case 4)
				rotateRight(x->parent);				// right rotation on x's parent    (case 4)
				x				= root;			// x is now the root			(case 4)
			}
		}
	}
	x->color = false;								// color x (the root) BLACK		(exit)

	return;
}

// ******** Rotation Functions ****************************************************************************

void splittree::rotateLeft(elementsp *x) {
	elementsp *y;
	// do pointer-swapping operations for left-rotation
	y               = x->right;					// grab right child
	x->right        = y->left;					// make x's RIGHT-CHILD be y's LEFT-CHILD
	y->left->parent = x;						// make x be y's LEFT-CHILD's parent
	y->parent       = x->parent;					// make y's new parent be x's old parent

	if (x->parent==NULL) { root = y; }				// if x was root, make y root
	else {									// 
		if (x == x->parent->left)				// if x is LEFT-CHILD, make y be x's parent's
			{ x->parent->left  = y; }			//    left-child
		else { x->parent->right = y; }				//    right-child
	}										// 
	y->left   = x;							// make x be y's LEFT-CHILD
	x->parent = y;							// make y be x's parent
	
	return;
}

void splittree::rotateRight(elementsp *y) {
	elementsp *x;
	// do pointer-swapping operations for right-rotation
	x                = y->left;					// grab left child
	y->left          = x->right;					// replace left child yith x's right subtree
	x->right->parent = y;						// replace y as x's right subtree's parent
	
	x->parent        = y->parent;					// make x's new parent be y's old parent
	if (y->parent==NULL) { root = x; }				// if y was root, make x root
	else {
		if (y == y->parent->right)				// if y is RIGHT-CHILD, make x be y's parent's
			{ y->parent->right  = x; }			//    right-child
		else { y->parent->left   = x; }				//    left-child
	}
	x->right  = y;							// make y be x's RIGHT-CHILD
	y->parent = x;							// make x be y's parent
	
	return;
}

// ******** Display Functions *****************************************************************************
// public
void splittree::printTree() {
//	cout << "\nTREE SIZE = " << support << "\t" << total_weight << "\t" << total_count << endl;
//	cout << "# "; printSubTree(root);
	return;
}

// private
void splittree::printSubTree(elementsp *z) {
	if (z==leaf) { return; }
	else {
//		cout << "(" << z->split << " " << z->weight << " " << z->count << " " << z->color << ")"<<endl;
//		cout << "L "; printSubTree(z->left); cout << endl;
//		cout << "R "; printSubTree(z->right); cout << endl;
	}
	return;
}

// public
void	splittree::printTreeAsList() {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
	int temp = 0;
	while (curr != NULL) {
		if (curr->y > 5) {
//			cout << curr->x;
//			if (curr->y >= 0.5*total_weight) { cout << "\t* "; temp++; } else { cout << "\t  "; }
//			cout << curr->y << "\t" << curr->c / total_weight << "\n";
		}
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
//	cout << "total_count   = " << total_count  << endl;
//	cout << "total_weight  = " << total_weight << endl;
//	cout << "total_size    = " << support << endl;
//	cout << "weight >= 1/2 = " << temp << endl;

	return;
}

void	splittree::printTreeAsShortList() {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
//	cout << "numSplits = " << support << endl;
	int temp = 0;
	while (curr != NULL) {
		if (curr->y / total_weight >= 0.5) {
//			cout << curr->x << "\t* ";
			temp++;
//			cout << curr->y << "\t" << curr->y / total_weight << "\n";
		}
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
//	cout << "maximum = " << total_weight << endl;
//	cout << "weight >= 1/2 = " << temp << endl;
	
	return;
}

void splittree::recordTreeAsList(const string file_out, const double threshold) {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
	
	ofstream fout(file_out.c_str(), ios::trunc);
	fout.precision(8);
//	for (int i=0; i<len; i++) { fout << "*"; } fout << "\t" << total_weight << "\n";
	while (curr != NULL) {
		if (curr->y / total_weight >= threshold) { fout << curr->x << "\t" << curr->y / total_weight << "\n"; }
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
	fout.close();
	
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

#endif
