//
//  kmerTree.h
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pdsBayesian_kmerTree_h
#define pdsBayesian_kmerTree_h


class KmerTree : public TaxonomyTree {
	
public:
	KmerTree(string, string, int);
	void addTaxonomyToTree(string, string);
	
private:
	int kmerSize;
	
};

#endif
