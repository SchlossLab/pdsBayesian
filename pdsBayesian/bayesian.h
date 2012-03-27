#ifndef BAYESIAN
#define BAYESIAN

/*
 *  bayesian.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/3/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

/**************************************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string> 
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>

/**************************************************************************************************/

using namespace std;

/**************************************************************************************************/

inline void gobble(istream& f){
	
	char d;
    while(isspace(d=f.get()))		{;}
	f.putback(d);
	//
}

/**************************************************************************************************/

#endif
