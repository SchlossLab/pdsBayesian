//
//  kmerTree.h
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pdsBayesian_kmerTree_h
#define pdsBayesian_kmerTree_h

class KmerNode;

class KmerTree : public TaxonomyTree {
	
public:
	KmerTree(string, string, int);
	~KmerTree();
	void addTaxonomyToTree(string, vector<int>&);
	void classifyQuery(string, string, string&);//, string&);

private:
	string deGap(string);
	vector<int> ripKmerProfile(string);
	int getMinRiskIndexKmer(vector<int>&, vector<int>&, vector<double>&);
	void aggregateThetas();
	int sanityCheck(vector<vector<int> >&, vector<int>&);

	int kmerSize;
	int numPossibleKmers;
	vector<KmerNode*> tree;

};

#endif
