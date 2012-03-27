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

TaxonomyTree::TaxonomyTree(string referenceFileName, string taxonomyFileName, string type) : method(type) {

	TaxonomyNode* newNode;
	if(method == "align"){
		newNode = new AlignNode("Root", 0);
	}
	else if(method == "kmer"){
		newNode = new KmerNode("Root", 0, 8);
	}
	
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
	referenceFile.close();
	
	numTaxa = tree.size();
	
	numLevels = 0;
	for(int i=0;i<numTaxa;i++){
		int level = tree[i]->getLevel();
		if(level > numLevels){	numLevels = level;	}
        tree[i]->checkTheta();
	}
	numLevels++;
	
    
}

/**************************************************************************************************/

TaxonomyTree::~TaxonomyTree(){
	
	for(int i=0;i<tree.size();i++){
		delete tree[i];
	}
	
}	

/**************************************************************************************************/

void TaxonomyTree::addTaxonomyToTree(string taxonomy, string sequence){
	
	TaxonomyNode* newNode;
	string taxonName = "";
	int treePosition = 0;							//	the root is element 0
	
	tree[treePosition]->loadSequence(sequence);		//	add sequence to the node to build thetas

	int level = 1;
	
	for(int i=0;i<taxonomy.length();i++){			//	step through taxonomy string...
		
		if(taxonomy[i] == ';'){						//	looking for semicolons...
			
			int newIndex = tree[treePosition]->getChildIndex(taxonName);	//	look to see if your current node already
			//	has a child with the new taxonName
			if(newIndex != -1)	{	treePosition = newIndex;	}		//	if you've seen it before, jump to that
			else {														//	 position in the tree
				int newChildIndex = tree.size();						//	otherwise, we'll have to create one...
				tree[treePosition]->makeChild(taxonName, newChildIndex);

				if(method == "align"){
					newNode = new AlignNode(taxonName, level);
				}
				else if(method == "kmer"){
					newNode = new KmerNode(taxonName, level, 8);
				}
				
				newNode->setParent(treePosition);
				
				tree.push_back(newNode);
				treePosition = newChildIndex;
			}
			
			tree[treePosition]->loadSequence(sequence);	//	now that we've gotten to the correct node, add the
			//	sequence data to that node to update that node's theta - seems slow...				
			
			taxonName = "";								//	clear out the taxon name that we will build as we look 
			level++;
		}												//	for a semicolon
		else{
			taxonName += taxonomy[i];					//	keep adding letters until we reach a semicolon
		}
	}
}

/**************************************************************************************************/

double TaxonomyTree::getOutlierLogProbability(string sequence){
	
	double count = 0;
	
	for(int i=0;i<sequence.length();i++){
		
		if(sequence[i] != '.'){	count++;	}
		
	}
	
	return count * log(0.2);
}

/**************************************************************************************************/

double TaxonomyTree::getLogExpSum(vector<double> probabilities, int& maxIndex){
	
//	http://jblevins.org/notes/log-sum-exp
	
	double maxProb = probabilities[0];
	maxIndex = 0;

	int numProbs = probabilities.size();
	
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
	
	int numProbs = probabilities.size();
	
	vector<double> G(numProbs, 0.2);	//a random sequence will, on average, be 20% similar to any other sequence
	vector<double> risk(numProbs, 0);

	for(int i=1;i<numProbs;i++){ //use if you want the outlier group
//	for(int i=0;i<numProbs;i++){
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

void TaxonomyTree::classifyQuery(string seqName, string querySequence, string& taxonProbabilityString, string& levelProbabilityString){
	
	double logPOutlier = getOutlierLogProbability(querySequence);

	vector<vector<double> > pXgivenKj_D_j(numLevels);
	vector<vector<int> > indices(numLevels);
	for(int i=0;i<numLevels;i++){               //comment out to remove outlier group
		pXgivenKj_D_j[i].push_back(logPOutlier);
		indices[i].push_back(-1);
	}
	
	
	for(int i=0;i<numTaxa;i++){
		pXgivenKj_D_j[tree[i]->getLevel()].push_back(tree[i]->getPxGivenkj_D_j(querySequence));
		indices[tree[i]->getLevel()].push_back(i);
	}
	
	vector<double> levelProbability(numLevels, 0);
	vector<double> bestPosterior(numLevels, 0);
	vector<int> maxIndex(numLevels, 0);
	

	for(int i=0;i<numLevels;i++){
		
		levelProbability[i] = getLogExpSum(pXgivenKj_D_j[i], maxIndex[i]);
		
		int numTaxaInLevel = indices[i].size();
		
		vector<double> posteriors(numTaxaInLevel, 0);
		for(int j=0;j<numTaxaInLevel;j++){
			posteriors[j] = exp(pXgivenKj_D_j[i][j] - levelProbability[i]);
		}
		
		maxIndex[i] = getMinRiskIndex(querySequence, indices[i], posteriors);

		bestPosterior[i] = posteriors[maxIndex[i]];
	}
	int maxLevel = 0;

	double allLevelSum = getLogExpSum(levelProbability, maxLevel);

	for(int i=0;i<numLevels;i++){	levelProbability[i] -= log(indices[i].size());	}
	
	allLevelSum = getLogExpSum(levelProbability, maxLevel);

	double levelPosterior = exp(levelProbability[maxLevel] - allLevelSum);
	

	sanityCheck(indices, maxIndex, maxLevel);
	
	stringstream taxonProbabilityOutput;
	taxonProbabilityOutput.setf(ios::fixed, ios::floatfield);
	taxonProbabilityOutput.setf(ios::showpoint);
    
    stringstream levelProbabilityOutput;
	levelProbabilityOutput.setf(ios::fixed, ios::floatfield);
	levelProbabilityOutput.setf(ios::showpoint);

    
	taxonProbabilityOutput << seqName << '(' << maxLevel << ';' << levelPosterior << ')' << '\t';
	levelProbabilityOutput << seqName << '(' << maxLevel << ';' << levelPosterior << ')' << '\t';

	for(int i=1;i<numLevels;i++){
		if(indices[i][maxIndex[i]] != -1){
			taxonProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(1) << bestPosterior[i] * 100 << ");";
			levelProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(1) << exp(levelProbability[i] - allLevelSum) << ");";
		}
		else{
			taxonProbabilityOutput << "incertae_sedis" << '(' << setprecision(1) << bestPosterior[i] * 100 << ");";
			levelProbabilityOutput << "incertae_sedis" << '(' << setprecision(1) << exp(levelProbability[i] - allLevelSum) << ");";
		}

	}
	
	taxonProbabilityString = taxonProbabilityOutput.str();
	levelProbabilityString = levelProbabilityOutput.str();

}

/**************************************************************************************************/
