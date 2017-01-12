#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pdbio.h"
#include "xyzio.h"
#include "dcdio.h"
#include "md.h"
#include "ran2.h"

int rseed = -1223894;

#define T 300 //K
#define Kb 8.314462e-3 //In kJ/K*mol
#define tau 0.001
#define gamma 50.0

int main(int argc, char* argv[]){
	XYZ xyz;
	readXYZ("Atomscoor.xyz", &xyz);
	
	MDSystem mds;
	
	mds.N = xyz.atomCount;
	mds.r = (float3*)calloc(mds.N, sizeof(float3));
	mds.v = (float3*)calloc(mds.N, sizeof(float3));
	mds.f = (float3*)calloc(mds.N, sizeof(float3));
	mds.m = (float*)calloc(mds.N, sizeof(float));

	int i;
	for(i = 0; i < mds.N; i++){
		mds.r[i].x = xyz.atoms[i].x;
		mds.r[i].y = xyz.atoms[i].y;
		mds.r[i].z = xyz.atoms[i].z;
	}

	for(i = 0; i < mds.N; i++){
		mds.m[i] = 1.0;
	}
	
	for(i = 0; i < mds.N; i++){
		float mult = sqrtf(Kb*T/mds.m[i]);
		mds.v[i].x = mult*gasdev(&rseed);
		mds.v[i].y = mult*gasdev(&rseed);
		mds.v[i].z = mult*gasdev(&rseed);
	}

	DCD dcd;	
	createDCD(&dcd, mds.N, 11, 0, 1.0, 1, 0, 0, 0, 0);
	dcdOpenWrite(&dcd, "dynam.dcd");
	dcdWriteHeader(dcd);	

	float var = sqrtf(Kb*T*2.0*gamma/tau);
	
	long int step;
	for(step = 0; step < 1000000; step++){

		for(i = 0; i < mds.N; i++){

			mds.f[i].x = 0.0;
			mds.f[i].y = 0.0;
			mds.f[i].z = 0.0;

			mds.f[i].x += -gamma*mds.v[i].x + var*gasdev(&rseed);
			mds.f[i].y += -gamma*mds.v[i].y + var*gasdev(&rseed);
			mds.f[i].z += -gamma*mds.v[i].z + var*gasdev(&rseed);		

			mds.v[i].x += tau*mds.f[i].x/mds.m[i];
			mds.v[i].y += tau*mds.f[i].y/mds.m[i];
			mds.v[i].z += tau*mds.f[i].z/mds.m[i];

			mds.r[i].x += mds.v[i].x*tau;
			mds.r[i].y += mds.v[i].y*tau;
			mds.r[i].z += mds.v[i].z*tau;
		}
		
		if(step % 10000 == 0){
			for(i = 0; i < mds.N; i++){
				dcd.frame.X[i] = 10.0*mds.r[i].x;
				dcd.frame.Y[i] = 10.0*mds.r[i].y;
				dcd.frame.Z[i] = 10.0*mds.r[i].z;
			}
			dcdWriteFrame(dcd);
		}
	}


}
