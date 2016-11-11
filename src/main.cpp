/*
 * main.cpp
 *
 *  Created on: 11.08.2016
 *      Author: klyshko
 */

#include <cmath>
#include <stdio.h>
#include <vector>

#define MT_FRAGMENT_LEN 5.0

struct {
	float x;
	float y;
	float z;
} Coord;

Coord normalize(Coord r){
	float len = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);	
	r.x /= len;
	r.y /= len;
	r.z /= len;  
	return r;
}

class MT{
	public:
		Coord dir;
		int phase;
		float length;

		
		MT(float x, float y, float z) : dir.x(x), dir.y(y), dir.z(z){
			phase = 1.0;
			length = 0.0;
		}
		
		void changeLength(){
			if (phase = 1.0){
				length += MT_FRAGMENT_LEN;
`			} else {
				length -= MT_FRAGMENT_LEN;
			}	
		}

		Coord getTipCoord(){
			Coord r;
			r.x = self.dir.x * length;
			r.y = self.dir.y * length;
			r.z = self.dir.z * length;
			return r;
		}

};

class Pole{
	public:
		float radius;
		Coord coord;
		vector <MT*> tubes;
		int num_tubes; 
		
		Pole(){
			coord.x = 0.0;
			coord.y = 0.0;
			coord.z = 0.0;
			radius = 1.0;
			num_tubes = 1;
		}
};


int main(int argc, char *argv[]){

    return 0;
}










