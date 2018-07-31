/*
 *  definitions.h
 *  twpost
 *
 *  Created by Daniel Gordon on 10/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <valarray>
#include <string>

#define THIS_MACHINE_LITTLE_ENDIAN
//#define THIS_MACHINE_BIG_ENDIAN

using namespace std;

typedef complex<float> complx;
typedef int64_t int64;
const float pi = 3.1415926;
const complx ii = complx(0,1);
