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
	
	cout << name << endl;
	for(int i=0;i<numPossibleKmers;i++){
		if(kmerVector[i] != 0){
			cout << getKmerBases(i) << '\t' << kmerVector[i] << endl;
		}
	}
	cout << endl;	
	
}

/**********************************************************************************************************************/

string KmerNode::getKmerBases(int kmerNumber){
	
	//	Here we convert the kmer number into the kmer in terms of bases.
	//
	//	Example:	Score = 915 (for a 6-mer)
	//				Base6 = (915 / 4^0) % 4 = 915 % 4 = 3 => T	[T]
	//				Base5 = (915 / 4^1) % 4 = 228 % 4 = 0 => A	[AT]
	//				Base4 = (915 / 4^2) % 4 = 57 % 4 = 1 => C	[CAT]
	//				Base3 = (915 / 4^3) % 4 = 14 % 4 = 2 => G	[GCAT]
	//				Base2 = (915 / 4^4) % 4 = 3 % 4 = 3 => T	[TGCAT]
	//				Base1 = (915 / 4^5) % 4 = 0 % 4 = 0 => A	[ATGCAT] -> this checks out with the previous method
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
	
	string kmer = "";
	
	if(kmerNumber == power4s[kmerSize]){//pow(4.,7)){	//	if the kmer number is the same as the maxKmer then it must
		for(int i=0;i<kmerSize;i++){					//	have had an N in it and so we'll just call it N x kmerSize
			kmer += 'N';
		}
	}
	else{
		for(int i=0;i<kmerSize;i++){
			int nt = (int)(kmerNumber / (float)power4s[i]) % 4;		//	the '%' operator returns the remainder 
			if(nt == 0)		{	kmer = 'A' + kmer;	}				//	from int-based division ]
			else if(nt == 1){	kmer = 'C' + kmer;	}
			else if(nt == 2){	kmer = 'G' + kmer;	}
			else if(nt == 3){	kmer = 'T' + kmer;	}
		}
	}
	return kmer;
}

/**********************************************************************************************************************/

vector<int> KmerNode::ripKmerProfile(string sequence){

	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };

//	assume all input sequences are unaligned
//	string unalignSequence;
//	
//	int alignLength = (int)alignSequence.length();	//	this function runs through the alignment and increments the frequency
//												//	of each base for a particular taxon.  we are building the thetas
//	
//	for(int i=0;i<alignLength;i++){
//		if(alignSequence[i] != '.' && alignSequence[i] != '-'){
//			unalignSequence += alignSequence[i];			
//		}
//	}

	
	int nKmers = (int)sequence.length() - kmerSize + 1;
	
	vector<int> kmerProfile(numPossibleKmers + 1, 0);
	
	for(int i=0;i<nKmers;i++){
		int kmer = 0;
		for(int j=0;j<kmerSize;j++){
			if(toupper(sequence[j+i]) == 'A')		{	kmer += (0 * power4s[kmerSize-j-1]);	}
			else if(toupper(sequence[j+i]) == 'C')	{	kmer += (1 * power4s[kmerSize-j-1]);	}
			else if(toupper(sequence[j+i]) == 'G')	{	kmer += (2 * power4s[kmerSize-j-1]);	}
			else if(toupper(sequence[j+i]) == 'U')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
			else if(toupper(sequence[j+i]) == 'T')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
			else									{	kmer = power4s[kmerSize]; j = kmerSize;	}
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
	double alpha = 1.0 / (double)totalSeqs;

	for(int i=0;i<numPossibleKmers;i++){
		
		if(queryKmerProfile[i] != 0){									//numUniqueKmers needs to be the value from Root;
			sumLogProb += log((kmerVector[i] + alpha) / (numSeqs + numUniqueKmers * alpha));
		}
		
	}
	return sumLogProb;

}

/**********************************************************************************************************************/
