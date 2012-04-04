/**************************************************************************************************/

/*
 *  alignBayesian.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/3/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "bayesian.h"
#include "taxonomynode.h"
#include "taxonomytree.h"
#include "aligntree.h"
#include "kmertree.h"


/**************************************************************************************************/

int main(int argc, char *argv[]){
	
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);

	string taxonomyFileName;	//	the name of the file that contains the taxonomy data
	string referenceFileName;	//	the name of the file that contains the reference sequences
	string queryFileName;		//	the name of the file that contains the query sequences, aligned on the
								//	  same basis as the reference sequences
	int kmerSize = 8;
	string method = "align";
	
	if(argc > 1) {
		for(char **p=argv+1;p<argv+argc;p++) {
			if(strcmp(*p,"-tax")==0) {
				if(++p>=argv+argc){}
				istringstream f(*p);
				if(!(f >> taxonomyFileName)){}
				if(taxonomyFileName=="") {
					cerr << "Error: must provide taxonomy file." << endl;
				}
			}
			else if(strcmp(*p,"-ref")==0) {
				if(++p>=argv+argc){}
				istringstream f(*p);
				if(!(f >> referenceFileName)){}
				if(referenceFileName=="") {
					cerr << "Error: must provide reference file." << endl;
				}
			}
			else if(strcmp(*p,"-query")==0) {
				if(++p>=argv+argc){}
					istringstream f(*p);
				if(!(f >> queryFileName)){}
				if(queryFileName=="") {
					cerr << "Error: must provide query file." << endl;
				}
			}
			else if(strcmp(*p,"-method")==0) {
				if(++p>=argv+argc){}
				istringstream f(*p);
				if(!(f >> method)){}
				if(queryFileName=="") {
					cerr << "Error: must provide query file." << endl;
				}
			}
			else if(strcmp(*p,"-ksize")==0) {
				if(++p>=argv+argc){}
				istringstream f(*p);
				if(!(f >> kmerSize)){}
				if(referenceFileName=="") {
					cerr << "Error: must provide reference file." << endl;
				}
			}
		}
	}
		
	cout << "Building database..." << endl;
	TaxonomyTree* database;
	
	if(method == "align"){
		database = new AlignTree(referenceFileName, taxonomyFileName);	//build the tree structure that holds the reference
	}
	else if(method == "kmer"){
		database = new KmerTree(referenceFileName, taxonomyFileName, kmerSize);	//build the tree structure that holds the reference
	}
	
	string taxProbFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".tprob.taxonomy";
	ofstream taxProbFile(taxProbFileName.c_str());

    string levelProbFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".lprob.taxonomy";
	ofstream levelProbFile(levelProbFileName.c_str());

	int index = 1;
	
	cout << "Classifying sequences..." << endl;
	ifstream queryFile(queryFileName.c_str());
	while(queryFile){
			
		string seqName;
		string sequence;
		string taxonProbOutput, levelProbOutput;
		
		queryFile >> seqName >> sequence;

		seqName = seqName.substr(1);
		
		gobble(queryFile);

		database->classifyQuery(seqName, sequence, taxonProbOutput, levelProbOutput);
		
		taxProbFile << taxonProbOutput << endl;
        levelProbFile << levelProbOutput << endl;

		if(index % 100 == 0){	cout << index << endl;	}
		index++;
		
	}
	queryFile.close();
	taxProbFile.close();
	levelProbFile.close();
}

/**************************************************************************************************/
