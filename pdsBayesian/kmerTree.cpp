//
//  kmerTree.cpp
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include <iostream>
#include "taxonomytree.h"
#include "kmerTree.h"
#include "kmerNode.h"

/**************************************************************************************************/

KmerTree::KmerTree(string referenceFileName, string taxonomyFileName, int k) : kmerSize(k){
	
	TaxonomyNode* newNode = new KmerNode("Root", 0, kmerSize);
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
	
	numTaxa = (int)tree.size();
	
	numLevels = 0;
	for(int i=0;i<numTaxa;i++){
		int level = tree[i]->getLevel();
		if(level > numLevels){	numLevels = level;	}
        tree[i]->checkTheta();
		tree[i]->setNumUniqueKmers(tree[0]->getNumUniqueKmers());
	}
	numLevels++;

}

/**************************************************************************************************/

void KmerTree::addTaxonomyToTree(string taxonomy, string sequence){
	
	TaxonomyNode* newNode;
	string taxonName = "";
	int treePosition = 0;							//	the root is element 0
	
	tree[treePosition]->loadSequence(sequence);		//	add sequence to the node to build thetas
	
	int level = 1;
	
	for(int i=0;i<taxonomy.length();i++){			//	step through taxonomy string...
		
		if(taxonomy[i] == ';'){						//	looking for semicolons...
			
			int newIndex = tree[treePosition]->getChildIndex(taxonName);//	look to see if your current node already
																		//	   has a child with the new taxonName
			if(newIndex != -1)	{	treePosition = newIndex;	}		//	if you've seen it before, jump to that
			else {														//	   position in the tree
				int newChildIndex = (int)tree.size();					//	otherwise, we'll have to create one...
				tree[treePosition]->makeChild(taxonName, newChildIndex);
				
				newNode = new KmerNode(taxonName, level, kmerSize);
				
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

