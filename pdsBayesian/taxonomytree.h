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
	TaxonomyTree(string, string, string);
	~TaxonomyTree();
	void addTaxonomyToTree(string, string);
	void classifyQuery(string, string, string&, string&);

private:
	double getOutlierLogProbability(string);
	double getLogExpSum(vector<double>, int&);
	int getMinRiskIndex(string, vector<int>, vector<double>);
	void sanityCheck(vector<vector<int> >, vector<int>, int&);
				   
	vector<TaxonomyNode*> tree;
	int numTaxa;
	int numLevels;
	string method;
};

#endif
