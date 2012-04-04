#ifndef KMERNODE
#define KMERNODE

/*
 *  kmerNode.h
 *  bayesian
 *
 *  Created by Pat Schloss on 10/11/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "taxonomynode.h"

/**********************************************************************************************************************/

class KmerNode : public TaxonomyNode {
	
public:
	KmerNode(string, int, int);
	void loadSequence(string);
	void printTheta();
	double getPxGivenkj_D_j(string query);
	double getSimToConsensus(string query);;
	void checkTheta(){};
	void setNumUniqueKmers(int num)	{	numUniqueKmers = num;	}
	int getNumUniqueKmers()			{	return numUniqueKmers;	}


private:
	vector<int> ripKmerProfile(string);
	string getKmerBases(int);
	
	int kmerSize;								//	value of k
	int numPossibleKmers;						//	4^kmerSize
	int numUniqueKmers;							//	number of unique kmers seen in a group ~ O_kj
	int numKmers;								//	number of kmers in a sequence
	vector<int> kmerVector;						//	counts of kmers across all sequences in a node
};

/**********************************************************************************************************************/

#endif

