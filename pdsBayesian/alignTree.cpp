//
//  alignTree.cpp
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include <iostream>
#include "alignNode.h"
#include "taxonomytree.h"
#include "alignTree.h"

/**************************************************************************************************/

AlignTree::AlignTree(string referenceFileName, string taxonomyFileName){
	
	AlignNode* newNode = new AlignNode("Root", 0);
	tree.push_back(newNode);			//	the tree is stored as a vector of elements of type TaxonomyNode

	string refName, refTaxonomy, refSequence;
	
	map<string, string> taxonomyData;
	ifstream taxonomyFile(taxonomyFileName.c_str());
	while(taxonomyFile){
		taxonomyFile >> refName >> refTaxonomy;
		taxonomyData[refName] = refTaxonomy;		//	create a map that links the reference sequence names to the
		gobble(taxonomyFile);						//	 taxonomy string for the sequence
	}
	taxonomyFile.close();
	
	ifstream referenceFile(referenceFileName.c_str());
	while(referenceFile){
		referenceFile >> refName >> refSequence;	//	read in fasta-formatted data - should replace with mothur command
		refName = refName.substr(1);				//	remove the ">" character in front of the sequence name
		refTaxonomy = taxonomyData[refName];		//	lookup the taxonomy string for the current reference sequence
		addTaxonomyToTree(refTaxonomy, refSequence);
		gobble(referenceFile);						//	removes extra whitespace
	}
	
	numTaxa = (int)tree.size();

	
	numLevels = 0;
	for(int i=0;i<numTaxa;i++){
		int level = tree[i]->getLevel();
		if(level > numLevels){	numLevels = level;	}
	}
	numLevels++;

	
	aggregateThetas();
	
	
//	tree[0]->printTheta();
//	cout << endl;
//	tree[1]->printTheta();
//	cout << endl;
	
	
	referenceFile.close();
	
	int dbSize = tree[0]->getNumSeqs();
	
	for(int i=0;i<numTaxa;i++){
        tree[i]->checkTheta();
		tree[i]->setTotalSeqs(dbSize);
	}
	
	cout << "done building" << endl;
}

/**************************************************************************************************/

AlignTree::~AlignTree(){
	
	for(int i=0;i<tree.size();i++){
		delete tree[i];
	}
	
}	

/**************************************************************************************************/

void AlignTree::addTaxonomyToTree(string& taxonomy, string& sequence){
	
	AlignNode* newNode;
	string taxonName = "";
	int treePosition = 0;							//	the root is element 0
	
	int level = 1;
	
	for(int i=0;i<taxonomy.length();i++){			//	step through taxonomy string...
		
		if(taxonomy[i] == ';'){						//	looking for semicolons...
			
			int newIndex = tree[treePosition]->getChildIndex(taxonName);	//	look to see if your current node already
			//	has a child with the new taxonName
			if(newIndex != -1)	{	treePosition = newIndex;	}		//	if you've seen it before, jump to that
			else {														//	 position in the tree
				int newChildIndex = (int)tree.size();						//	otherwise, we'll have to create one...
				tree[treePosition]->makeChild(taxonName, newChildIndex);
				
				newNode = new AlignNode(taxonName, level);
				
				newNode->setParent(treePosition);
				
				tree.push_back(newNode);
				treePosition = newChildIndex;
			}
			
			//	sequence data to that node to update that node's theta - seems slow...				
			taxonName = "";								//	clear out the taxon name that we will build as we look 
			level++;
		}												//	for a semicolon
		else{
			taxonName += taxonomy[i];					//	keep adding letters until we reach a semicolon
		}
	}
	tree[treePosition]->loadSequence(sequence);	//	now that we've gotten to the correct node, add the
	
}

/**************************************************************************************************/

void AlignTree::aggregateThetas(){
	
	vector<vector<int> > levelMatrix(numLevels+1);

	for(int i=0;i<tree.size();i++){
		levelMatrix[tree[i]->getLevel()].push_back(i);
	}
		
	for(int i=numLevels-1;i>0;i--){
		for(int j=0;j<levelMatrix[i].size();j++){
			
			AlignNode* holder = tree[levelMatrix[i][j]];
			
			tree[holder->getParent()]->addThetas(holder->getTheta(), holder->getNumSeqs());				
		}
		
	}
	
}

/**************************************************************************************************/

double AlignTree::getOutlierLogProbability(string& sequence){
	
	double count = 0;
	
	for(int i=0;i<sequence.length();i++){
		
		if(sequence[i] != '.'){	count++;	}
		
	}
	
	return count * log(0.2);
}

/**************************************************************************************************/

int AlignTree::getMinRiskIndexAlign(string& sequence, vector<int>& taxaIndices, vector<double>& probabilities){
	
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

int AlignTree::sanityCheck(vector<vector<int> >& indices, vector<int>& maxIndices){
	
	int finalLevel = (int)indices.size()-1;
	
	for(int position=1;position<indices.size();position++){
		int predictedParent = tree[indices[position][maxIndices[position]]]->getParent();
		int actualParent = indices[position-1][maxIndices[position-1]];
		
		if(predictedParent != actualParent){
			finalLevel = position - 1;
			return finalLevel;
		}
	}
	return finalLevel;
	
}

/**************************************************************************************************/

void AlignTree::classifyQuery(string seqName, string querySequence, string& taxonProbabilityString){//, string& levelProbabilityString){
	
	double logPOutlier = getOutlierLogProbability(querySequence);
	
	vector<vector<double> > pXgivenKj_D_j(numLevels);
	vector<vector<int> > indices(numLevels);
	for(int i=0;i<numLevels;i++){
		pXgivenKj_D_j[i].push_back(logPOutlier);
		indices[i].push_back(-1);
	}
	
	
	for(int i=0;i<numTaxa;i++){
//		cout << i << '\t' << tree[i]->getName() << '\t' << tree[i]->getLevel() << '\t' << tree[i]->getPxGivenkj_D_j(querySequence) << endl;
		
		pXgivenKj_D_j[tree[i]->getLevel()].push_back(tree[i]->getPxGivenkj_D_j(querySequence));
		indices[tree[i]->getLevel()].push_back(i);
	}
	
	vector<double> sumLikelihood(numLevels, 0);
	vector<double> bestPosterior(numLevels, 0);
	vector<int> maxIndex(numLevels, 0);
	int maxPosteriorIndex;
	
	
//	cout << "before best level" << endl;
	
	//let's find the best level and taxa within that level
	for(int i=0;i<numLevels;i++){ //go across all j's - from the root to genus
		
		int numTaxaInLevel = (int)indices[i].size();
		
//		cout << "numTaxaInLevel:\t" << numTaxaInLevel << endl;
		
		vector<double> posteriors(numTaxaInLevel, 0);		
		sumLikelihood[i] = getLogExpSum(pXgivenKj_D_j[i], maxPosteriorIndex);
		
		maxPosteriorIndex = 0;
		for(int j=0;j<numTaxaInLevel;j++){
			posteriors[j] = exp(pXgivenKj_D_j[i][j] - sumLikelihood[i]);
			
			if(posteriors[j] > posteriors[maxPosteriorIndex]){	
				maxPosteriorIndex = j;
			}
			
		}

		maxIndex[i] = getMinRiskIndexAlign(querySequence, indices[i], posteriors);
		
		maxIndex[i] = maxPosteriorIndex;
		bestPosterior[i] = posteriors[maxIndex[i]];	
	}
	
//	vector<double> pX_level(numLevels, 0);
//	
//	for(int i=0;i<numLevels;i++){
//		pX_level[i] = pXgivenKj_D_j[i][maxIndex[i]] - tree[indices[i][maxIndex[i]]]->getNumSeqs();
//	}
//	
//	int max_pLevel_X_index = -1;
//	double pX_level_sum = getLogExpSum(pX_level, max_pLevel_X_index);
//	double max_pLevel_X = exp(pX_level[max_pLevel_X_index] - pX_level_sum);
//	
//	vector<double> pLevel_X(numLevels, 0);
//	for(int i=0;i<numLevels;i++){
//		pLevel_X[i] = exp(pX_level[i] - pX_level_sum);
//	}
	
	
	
	
	int saneDepth = sanityCheck(indices, maxIndex);
	
	stringstream taxonProbabilityOutput;
	taxonProbabilityOutput.setf(ios::fixed, ios::floatfield);
	taxonProbabilityOutput.setf(ios::showpoint);
    
	//	stringstream levelProbabilityOutput;
	//	levelProbabilityOutput.setf(ios::fixed, ios::floatfield);
	//	levelProbabilityOutput.setf(ios::showpoint);
	
    
	taxonProbabilityOutput << seqName << '\t';
	//	taxonProbabilityOutput << seqName << '(' << max_pLevel_X_index << ';' << max_pLevel_X << ')' << '\t';
	//	levelProbabilityOutput << seqName << '(' << max_pLevel_X_index << ';' << max_pLevel_X << ')' << '\t';
	
	for(int i=1;i<=saneDepth;i++){
		if(indices[i][maxIndex[i]] != -1){
			taxonProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(6) << bestPosterior[i] << ");";
			//			levelProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(6) << pLevel_X[i] << ");";
		}
		else{
			taxonProbabilityOutput << "unclassified" << '(' << setprecision(6) << bestPosterior[i] << ");";
			//			levelProbabilityOutput << "unclassified" << '(' << setprecision(6) << pLevel_X[i] << ");";
		}
	}
	
	for(int i=saneDepth+1;i<numLevels;i++){
		taxonProbabilityOutput << "unclassified" << '(' << setprecision(6) << "0.00" << ");";
	}
	
	taxonProbabilityString = taxonProbabilityOutput.str();
//	levelProbabilityString = levelProbabilityOutput.str();
	
}


/**************************************************************************************************/
