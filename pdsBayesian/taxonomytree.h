#ifndef TAXONOMYTREE
#define TAXONOMYTREE

/*
 *  taxonomytree.h
 *  
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */


#include "taxonomynode.h"

class TaxonomyTree {
	
public:
	TaxonomyTree(){};
	virtual void classifyQuery(string, string, string&) = 0;//, string&) = 0;

protected:
	double getLogExpSum(vector<double>, int&);
	void classifyGeneric(string, string, float, string&, string&, string);
				   
	int numTaxa;
	int numLevels;
};

#endif
