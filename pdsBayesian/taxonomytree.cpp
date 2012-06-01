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

TaxonomyTree::~TaxonomyTree(){
	
	for(int i=0;i<tree.size();i++){
		delete tree[i];
	}
	
}	

/**************************************************************************************************/

double TaxonomyTree::getLogExpSum(vector<double> probabilities, int& maxIndex){
	
//	http://jblevins.org/notes/log-sum-exp
	
	double maxProb = probabilities[0];
	maxIndex = 0;

	int numProbs = (int)probabilities.size();
	
	for(int i=1;i<numProbs;i++){
		if(probabilities[i] > maxProb){
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

int TaxonomyTree::getMinRiskIndex(string sequence, vector<int> taxaIndices, vector<double> probabilities){
	
	int numProbs = (int)probabilities.size();
	
	vector<double> G(numProbs, 0.2);	//a random sequence will, on average, be 20% similar to any other sequence
	vector<double> risk(numProbs, 0);

	for(int i=1;i<numProbs;i++){ //use if you want the outlier group
		G[i] = tree[taxaIndices[i]]->getSimToConsensus(sequence);
	}
	
	double minRisk = 1e6;
	int minRiskIndex = 0;
	
	for(int i=0;i<numProbs;i++){
		
		for(int j=0;j<numProbs;j++){
			if(i != j){
				risk[i] += probabilities[j] * G[j];
			}			
		}
		
		if(risk[i] < minRisk){
			minRisk = risk[i];
			minRiskIndex = i;
		}
	}
	
	return minRiskIndex;
}

/**************************************************************************************************/

void TaxonomyTree::sanityCheck(vector<vector<int> > indices, vector<int> maxIndices, int& maxLevel){

	for(int position=1;position<=maxLevel;position++){
		int predictedParent = tree[indices[position][maxIndices[position]]]->getParent();
		int actualParent = indices[position-1][maxIndices[position-1]];
		
		if(predictedParent != actualParent){
			maxLevel = position - 1;
			return;
		}
	}
}

/**************************************************************************************************/

void TaxonomyTree::classifyGeneric(string seqName, string querySequence, float logPOutlier, string& taxonProbabilityString, string& levelProbabilityString){
	

	vector<vector<double> > pXgivenKj_D_j(numLevels);
	vector<vector<int> > indices(numLevels);
	for(int i=0;i<numLevels;i++){
		pXgivenKj_D_j[i].push_back(logPOutlier);
		indices[i].push_back(-1);
	}
	
	
	for(int i=0;i<numTaxa;i++){
		pXgivenKj_D_j[tree[i]->getLevel()].push_back(tree[i]->getPxGivenkj_D_j(querySequence));
		indices[tree[i]->getLevel()].push_back(i);
	}
	
	vector<double> sumLikelihood(numLevels, 0);
	vector<double> bestPosterior(numLevels, 0);
	vector<int> maxIndex(numLevels, 0);
	int maxPosteriorIndex;
	
	
	//let's find the best level and taxa within that level
	for(int i=0;i<numLevels;i++){ //go across all j's - from the root to genus

		int numTaxaInLevel = (int)indices[i].size();
		
		vector<double> posteriors(numTaxaInLevel, 0);		
		sumLikelihood[i] = getLogExpSum(pXgivenKj_D_j[i], maxPosteriorIndex);
		
		maxPosteriorIndex = 0;
		for(int j=0;j<numTaxaInLevel;j++){
			posteriors[j] = exp(pXgivenKj_D_j[i][j] - sumLikelihood[i]);
			if(posteriors[j] > posteriors[maxPosteriorIndex]){	
				maxPosteriorIndex = j;
			}
				
		}
		
		
		//maxIndex[i] = getMinRiskIndex(querySequence, indices[i], posteriors);	//should perhaps fix this in the future
		maxIndex[i] = maxPosteriorIndex;
		bestPosterior[i] = posteriors[maxIndex[i]];	
	}

	vector<double> pX_level(numLevels, 0);

	for(int i=0;i<numLevels;i++){
		pX_level[i] = pXgivenKj_D_j[i][maxIndex[i]] - tree[indices[i][maxIndex[i]]]->getNumSeqs();
	}

	int max_pLevel_X_index = -1;
	double pX_level_sum = getLogExpSum(pX_level, max_pLevel_X_index);
	double max_pLevel_X = exp(pX_level[max_pLevel_X_index] - pX_level_sum);
	
	vector<double> pLevel_X(numLevels, 0);
	for(int i=0;i<numLevels;i++){
		pLevel_X[i] = exp(pX_level[i] - pX_level_sum);
	}
	
	
	
	
	sanityCheck(indices, maxIndex, max_pLevel_X_index);

	
	stringstream taxonProbabilityOutput;
	taxonProbabilityOutput.setf(ios::fixed, ios::floatfield);
	taxonProbabilityOutput.setf(ios::showpoint);
    
    stringstream levelProbabilityOutput;
	levelProbabilityOutput.setf(ios::fixed, ios::floatfield);
	levelProbabilityOutput.setf(ios::showpoint);

    
	taxonProbabilityOutput << seqName << '(' << max_pLevel_X_index << ';' << max_pLevel_X << ')' << '\t';
	levelProbabilityOutput << seqName << '(' << max_pLevel_X_index << ';' << max_pLevel_X << ')' << '\t';

	for(int i=1;i<numLevels;i++){
		if(indices[i][maxIndex[i]] != -1){
			taxonProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(6) << bestPosterior[i] << ");";
			levelProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(6) << pLevel_X[i] << ");";
		}
		else{
			taxonProbabilityOutput << "incertae_sedis" << '(' << setprecision(6) << bestPosterior[i] << ");";
			levelProbabilityOutput << "incertae_sedis" << '(' << setprecision(6) << pLevel_X[i] << ");";
		}

	}
	
	taxonProbabilityString = taxonProbabilityOutput.str();
	levelProbabilityString = levelProbabilityOutput.str();

}

/**************************************************************************************************/
