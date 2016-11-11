/*
 * main.cpp
 *
 *  Created on: 11.08.2016
 *      Author: klyshko
 */

#include <cmath>
#include <stdio.h>
#include <vector>

#define MT_fragment_len 5.0

struct coord{
	
}


class MT{
	public:
		float x, y, z;
		int phase;

		
		MT(float x, float y, float z) : x(x), y(y), z(z){
			phase = 1.0;
		}


};

class Pole{
	public:
		float x, y, z, r;
		vector <MT> tubes;
		int num_tubes; 
		
		Pole(){
			x = 0.0;
			y = 0.0;
			z = 0.0;
			r = 1.0;
			num_tubes = 0;
		}
};






int main(int argc, char *argv[]){

    return 0;
}










