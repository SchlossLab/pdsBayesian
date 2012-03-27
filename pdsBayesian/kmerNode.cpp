/*
 *  kmerNode.cpp
 *  bayesian
 *
 *  Created by Pat Schloss on 10/11/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "kmerNode.h"


/**********************************************************************************************************************/

KmerNode::KmerNode(string s, int l, int n) : TaxonomyNode(s, l), kmerSize(n) {
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };

	numPossibleKmers = power4s[kmerSize];
	numUniqueKmers = 0;
	
	kmerVector.assign(numPossibleKmers, 0);

}

/**********************************************************************************************************************/

void KmerNode::loadSequence(string sequence){
	
	vector<int> kmerProfile = ripKmerProfile(sequence);
	
	for(int i=0;i<numPossibleKmers;i++){
		if(kmerVector[i] == 0 && kmerProfile[i] != 0)	{	numUniqueKmers++;	}
		
		kmerVector[i] += kmerProfile[i];
	}

	numSeqs++;
}	

/**********************************************************************************************************************/

void KmerNode::printTheta(){
	
	for(int i=0;i<numPossibleKmers;i++){
		cout << kmerVector[i] << '\t';
	}
	cout << endl;	
	
}

/**********************************************************************************************************************/

vector<int> KmerNode::ripKmerProfile(string alignSequence){

	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };

	string unalignSequence;
	
	int alignLength = alignSequence.length();	//	this function runs through the alignment and increments the frequency
												//	of each base for a particular taxon.  we are building the thetas
	
	for(int i=0;i<alignLength;i++){
		if(alignSequence[i] != '.' && alignSequence[i] != '-'){
			unalignSequence += alignSequence[i];			
		}
	}

	
	int nKmers = unalignSequence.length() - kmerSize + 1;
	
	vector<int> kmerProfile(numPossibleKmers + 1, 0);
	
	for(int i=0;i<nKmers;i++){
		int kmer = 0;
		for(int j=0;j<kmerSize;j++){
			if(toupper(unalignSequence[j+i]) == 'A')		{	kmer += (0 * power4s[kmerSize-j-1]);	}
			else if(toupper(unalignSequence[j+i]) == 'C')	{	kmer += (1 * power4s[kmerSize-j-1]);	}
			else if(toupper(unalignSequence[j+i]) == 'G')	{	kmer += (2 * power4s[kmerSize-j-1]);	}
			else if(toupper(unalignSequence[j+i]) == 'U')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
			else if(toupper(unalignSequence[j+i]) == 'T')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
			else											{	kmer = power4s[kmerSize]; j = kmerSize;	}
		}
		kmerProfile[kmer] = 1;
	}
		
	return kmerProfile;	
}

/**************************************************************************************************/

double KmerNode::getSimToConsensus(string query){
	
	vector<int> queryKmerProfile = ripKmerProfile(query);

	double sum = 0;
	double length = 0;
	
	for(int i=0;i<numPossibleKmers;i++){
		if(queryKmerProfile[i] != 0){
		
			sum += (double)kmerVector[i] / (double)numSeqs;
			length++;
			
		}
	}	
	
	return sum / length;
}

/**********************************************************************************************************************/

double KmerNode::getPxGivenkj_D_j(string query)	{	

	vector<int> queryKmerProfile = ripKmerProfile(query);
	
	double sumLogProb = 0.0000;
	
	double correction = pow((1.0/(double)numUniqueKmers), numSeqs) + 0.0001;
	
	for(int i=0;i<numPossibleKmers;i++){
		
		if(queryKmerProfile[i] != 0){			
			sumLogProb += log(kmerVector[i] / numSeqs + correction);
		}
		
	}
	
	sumLogProb /= (numSeqs + (double)numUniqueKmers * correction);
	
	return sumLogProb;
}

/**********************************************************************************************************************/
