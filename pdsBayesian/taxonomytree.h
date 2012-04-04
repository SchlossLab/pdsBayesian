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
	~TaxonomyTree();
	virtual void addTaxonomyToTree(string, string) = 0;
	void classifyQuery(string, string, string&, string&);

protected:
	double getOutlierLogProbability(string);
	double getLogExpSum(vector<double>, int&);
	int getMinRiskIndex(string, vector<int>, vector<double>);
	void sanityCheck(vector<vector<int> >, vector<int>, int&);
				   
	vector<TaxonomyNode*> tree;
	int numTaxa;
	int numLevels;
};

#endif
