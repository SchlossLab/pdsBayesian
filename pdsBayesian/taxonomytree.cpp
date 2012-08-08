/*
 *  taxonomytree.cpp
 *  
 *
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "taxonomytree.h"
#include "alignNode.h"
#include "kmerNode.h"

/**************************************************************************************************/

double TaxonomyTree::getLogExpSum(vector<double> probabilities, int& maxIndex){
	
//	http://jblevins.org/notes/log-sum-exp
	
	double maxProb = probabilities[0];
	maxIndex = 0;

	int numProbs = (int)probabilities.size();
	
	for(int i=1;i<numProbs;i++){
		if(probabilities[i] >= maxProb){
			maxProb = probabilities[i];
			maxIndex = i;
		}
	}
	
	double probSum = 0.0000;
	
	for(int i=0;i<numProbs;i++){
		probSum += exp(probabilities[i] - maxProb);		
	}

	probSum = log(probSum) + maxProb;
	
	return probSum;
	
}

/**************************************************************************************************/
