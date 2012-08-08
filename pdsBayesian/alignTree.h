//
//  alignTree.h
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pdsBayesian_alignTree_h
#define pdsBayesian_alignTree_h

class AlignNode;

class AlignTree : public TaxonomyTree {

public:
	AlignTree(string, string);
	~AlignTree();
	void addTaxonomyToTree(string&, string&);
	void classifyQuery(string, string, string&);//, string&);
	
private:
	double getOutlierLogProbability(string&);
	int getMinRiskIndexAlign(string&, vector<int>&, vector<double>&);
	void aggregateThetas();
	int sanityCheck(vector<vector<int> >&, vector<int>&);

	int numSeqs;
	vector<AlignNode*> tree;
};

#endif
