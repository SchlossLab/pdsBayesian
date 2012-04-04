#ifndef TAXONOMYNODE
#define TAXONOMYNODE

/*
 *  taxonomynode.h
 *  
 *
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "bayesian.h"

/**************************************************************************************************/

class TaxonomyNode {
	
public:
	TaxonomyNode();
	TaxonomyNode(string, int);
	void setName(string);
	string getName();


	void setParent(int);
	int getParent();
	
	void makeChild(string, int);
	map<string, int> getChildren();
	int getChildIndex(string);
	int	getNumKids();
	int getNumSeqs();
	int getLevel();
	
	virtual void loadSequence(string) = 0;
    virtual void checkTheta() = 0;
	virtual void printTheta() = 0;
	virtual double getPxGivenkj_D_j(string query) = 0;	//P(x | k_j, D, j)
	virtual double getSimToConsensus(string query) = 0;
	virtual void setNumUniqueKmers(int)	{};
	virtual int getNumUniqueKmers()		{	return 0;	};

private:
//	string name;
	int parent;
	map<string, int> children;
	int numChildren;
	int level;
	
protected:
	int numSeqs;
	string name;
	
};

/**************************************************************************************************/

#endif
