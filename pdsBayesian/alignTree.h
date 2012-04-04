//
//  alignTree.h
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pdsBayesian_alignTree_h
#define pdsBayesian_alignTree_h


class AlignTree : public TaxonomyTree {

public:
	AlignTree(string, string);
	void addTaxonomyToTree(string, string);
	
private:
	
};

#endif
