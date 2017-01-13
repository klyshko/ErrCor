#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "pdbio.h"
#include "psfio.h"
#include "dcdio.h"
#include "md.h"
#include "ran2.h"

using namespace std;


#define TOTAL_STEPS 100000
#define DCD_STRIDE	100 

#define MAX_MT_LENGTH	512	
#define MT_STATE_POL	1
#define MT_STATE_DEPOL	2
#define MT_NUM			750 
#define MT_RADIUS 		12.0	//nm
#define POLE_RADIUS		400.0	//nm
#define POLE_COUNT		1

#define MT_POL_RATE				0.007//7.e-8 //0.007	for 10 miliseconds	
#define MT_DEPOL_RATE			0.//0.17//1.7e-6//0.17		
#define MT_CATASTROPHE_FREQ 	45e-5//4.5e-10//45e-5	
#define MT_RESCUE_FREQ			0.0015//1.5e-8//0.0015

#define T				300.0	//K
#define KB 				8.314462e-3	   //kJ/K/mol
#define K_SPRING		10.0
#define K_THETA			155812.5
#define MEMBRANE_STIFF	1000.0

#define VISCOSITY		600 //xz which units and why?

#define ELLIPSE_A		8000
#define ELLIPSE_B		5000
#define ELLIPSE_C		5000

#define KIN_RADIUS			300
#define KIN_HEIGHT			500
#define KIN_DIST_0			800

//////////////////////////////////////////////////////////////////////////////////////

float dt = 100;	//ns

//////////////////////////////////////////////////////////////////////////////////////

int rseed = -123213;

//////////////////////////////////////////////////////////////////////////////////////


typedef struct{
	float3* r;
	int n;
	int state;
} MT;

typedef struct{
	float3* r;
	int n;
	int* state;
	float3 center;
} KT;

//////////////////////////////////////////////////////////////////////////////////////

MT* pol[2];
MDSystem mds;

//////////////////////////////////////////////////////////////////////////////////////

int writeXYZ(const char* filename, float3* data, int N, const char* modifier);
int initPSF();

//////////////////////////////////////////////////////////////////////////////////////

int computeGamma(float r);

//////////////////////////////////////////////////////////////////////////////////////

void kineticStep();
void membraneKinetics();
void membraneMechanics(float3* r, float3* f);
void computeHarmonic(float3* r, float3* f);
void computeAngles(float3* r, float3* f);
void integrateCPU(float3* r, float3* f, float dt, int N);
vector<float3> generateKinetochore();

//////////////////////////////////////////////////////////////////////////////////////

int main(){

	vector <float3> kinetochore = generateKinetochore();
	
	mds.N = 2 * MT_NUM*MAX_MT_LENGTH + kinetochore.size() - 1;
	mds.r = (float3*)calloc(mds.N, sizeof(float3));
	mds.f = (float3*)calloc(mds.N, sizeof(float3));
	pol[0] = (MT*)calloc(MT_NUM, sizeof(MT));
	pol[1] = (MT*)calloc(MT_NUM, sizeof(MT));
	int m, i, j;
	
	for (j = 0;  j < 2; j++){
		for(m = 0; m < MT_NUM; m++){
			pol[j][m].r = &mds.r[j * MT_NUM*MAX_MT_LENGTH + MAX_MT_LENGTH*m];
			pol[j][m].n = 2;
			pol[j][m].state = MT_STATE_POL;
		}	
	} 

	KT kt;
	kt.n = kinetochore.size() - 1;
	kt.r = (float3*)malloc(kt.n * sizeof(float3));
	kt.state = (int*)malloc(kt.n * sizeof(int));
	kt.center.x = kinetochore[kinetochore.size()-1].x;
	kt.center.y = kinetochore[kinetochore.size()-1].y;
	kt.center.z = kinetochore[kinetochore.size()-1].z;

	for (i = 0; i < kinetochore.size() - 1; i++){
		mds.r[2 * MT_NUM*MAX_MT_LENGTH + i].x = kinetochore[i].x;
		mds.r[2 * MT_NUM*MAX_MT_LENGTH + i].y = kinetochore[i].y;
		mds.r[2 * MT_NUM*MAX_MT_LENGTH + i].z = kinetochore[i].z;

		if (kt.r[i].x < kt.center.x){
			kt.state[i] = 0; //left part
		} else {
			kt.state[i] = 1; //right part
		}
	}
	kt.r = &mds.r[2 * MT_NUM*MAX_MT_LENGTH];

	
	float3 r[2];
	r[0].x = -3000; 
	r[0].y = 0; 
	r[0].z = 0;

	r[1].x = 3000;
	r[1].y = 0;
	r[1].z = 0;


	for (j = 0; j < 2; j++){
		for(m = 0; m < MT_NUM; m++){
		
			float theta = M_PI*ran2(&rseed);
			float phi = 2.0*M_PI*ran2(&rseed);
			pol[j][m].r[0].x = r[j].x + POLE_RADIUS*sin(theta)*cos(phi);
			pol[j][m].r[0].y = r[j].y + POLE_RADIUS*sin(theta)*sin(phi);
			pol[j][m].r[0].z = r[j].z + POLE_RADIUS*cos(theta);
			pol[j][m].r[1].x = r[j].x + (POLE_RADIUS + 2.0*MT_RADIUS)*sin(theta)*cos(phi);
			pol[j][m].r[1].y = r[j].y + (POLE_RADIUS + 2.0*MT_RADIUS)*sin(theta)*sin(phi);
			pol[j][m].r[1].z = r[j].z + (POLE_RADIUS + 2.0*MT_RADIUS)*cos(theta);
			float mindist = 100.0;
			int m1;
			for(m1 = 0; m1 < m; m1++){
				float dx = pol[j][m1].r[0].x - pol[j][m].r[0].x;
				float dy = pol[j][m1].r[0].y - pol[j][m].r[0].y;
				float dz = pol[j][m1].r[0].z - pol[j][m].r[0].z;
				float dr = sqrt(dx*dx + dy*dy + dz*dz);
				if(dr < mindist){
					mindist = dr;
				}
			}
			if(mindist < MT_RADIUS*2.0){
				m--;
			}
		}
		
	}

	initPSF();

	//writeXYZ("MT.xyz", mds.r, mds.N, "w");

	DCD dcd;	
	createDCD(&dcd, mds.N, 11, 0, 1.0, 1, 0, 0, 0, 0);
	dcdOpenWrite(&dcd, "output/dynam.dcd");
	dcdWriteHeader(dcd);	
	
	dcd.frame.X = (float*)calloc(mds.N, sizeof(float));
	dcd.frame.Y = (float*)calloc(mds.N, sizeof(float));
	dcd.frame.Z = (float*)calloc(mds.N, sizeof(float));
	
	long int step;

	for(i = 0; i < mds.N; i++){		
		dcd.frame.X[i] = mds.r[i].x;
		dcd.frame.Y[i] = mds.r[i].y;
		dcd.frame.Z[i] = mds.r[i].z;
	}
	dcdWriteFrame(dcd);

	for(step = 0; step < TOTAL_STEPS; step++){


		kineticStep();
		//membraneMechanics(mds.r, mds.f);
		membraneKinetics();
		computeHarmonic(mds.r, mds.f);
		computeAngles(mds.r, mds.f);
		integrateCPU(mds.r, mds.f, dt, mds.N);


		if(step % DCD_STRIDE == 0){

			for(i = 0; i < mds.N; i++){		
				dcd.frame.X[i] = mds.r[i].x;
				dcd.frame.Y[i] = mds.r[i].y;
				dcd.frame.Z[i] = mds.r[i].z;
			}
			dcdWriteFrame(dcd);
		}
		printf("%ld\n", step);
	}
	

	return 0;

}

vector<float3> generateKinetochore(){
	/*
	float xc = 2 * ran2(&rseed) * ELLIPSE_A - ELLIPSE_A;
	float yc = 2 * ran2(&rseed) * ELLIPSE_B - ELLIPSE_B;
	float zc = 2 * ran2(&rseed) * ELLIPSE_C - ELLIPSE_C;

	float3 dir;
	dir.x = 2 * ran2(&rseed) - 1.0;
	dir.y = 2 * ran2(&rseed) - 1.0;
	dir.z = 2 * ran2(&rseed) - 1.0;

	int nh = KIN_HEIGHT / (2 * MT_RADIUS);
	int nr = KIN_RADIUS / (2 * MT_RADIUS);
	int n = 0;

*/
	float xc = 0.0;
	float yc = 0.0;
	float zc = 0.0;

	float3 center;
	center.x = xc;
	center.y = yc;
	center.z = zc;

	float xinit_neg = xc - KIN_DIST_0 / 2;
	float xinit_pos = xc + KIN_DIST_0 / 2;
	float yinit = yc;// + KIN_DIST_0 / 2;
	float zinit = zc;// - KIN_HEIGHT / 2;

	printf("%f %f %f\n", xc, yc, zc);

	
	vector <float3> r;
	for (float z = zinit - KIN_HEIGHT / 2; z < zinit + KIN_HEIGHT / 2; z += 2 * MT_RADIUS ){

		for (float x = xinit_neg - KIN_RADIUS; x <= xinit_neg; x += 2 * MT_RADIUS){
			for (float y = yc - KIN_RADIUS; y < yc + KIN_RADIUS; y += 2 * MT_RADIUS){
				if ((x - xinit_neg) * (x - xinit_neg) + (y - yc) * (y - yc) < KIN_RADIUS * KIN_RADIUS){
					float3 coord;
					coord.z = z; 
					coord.x = x;
					coord.y = y;
					r.push_back(coord);
					printf("%f %f %f\n", x,y,z);
				}
			}
		}

		for (float x =  xinit_pos; x <= xinit_pos + KIN_RADIUS; x += 2 * MT_RADIUS){
			for (float y = yc - KIN_RADIUS; y < yc + KIN_RADIUS; y += 2 * MT_RADIUS){
				if ((x - xinit_pos) * (x - xinit_pos) + (y - yc) * (y - yc) < KIN_RADIUS * KIN_RADIUS){
					float3 coord;
					coord.z = z; 
					coord.x = x;
					coord.y = y;
					r.push_back(coord);
					printf("%f %f %f\n", x,y,z);
				}
			}
		}
	}
	r.push_back(center);

/*
	float3 *rarray = (float3*)malloc(r.size() * sizeof(float3));
	for (int i = 0; i < r.size(); i++){
		//printf("lol\t%f\t%f\t%f\n", r[i].x, r[i].y, r[i].z);
		rarray[i] = r[i];
	}

	//writeXYZ("cylinder.xyz", rarray, r.size(), "w");
*/
	return r;

}

int computeGamma(float r){
	return 6 * M_PI * VISCOSITY * r;
}

void membraneKinetics(){
	int m,j; 
	float rx, ry, rz, re2;
	for (j = 0; j < 2; j++){
		for(m = 0; m < MT_NUM; m++){

			if (pol[j][m].state != MT_STATE_DEPOL){
				rx = pol[j][m].r[pol[j][m].n - 1].x;
				ry = pol[j][m].r[pol[j][m].n - 1].y;
				rz = pol[j][m].r[pol[j][m].n - 1].z;

				re2 = pow(rx, 2) / pow(ELLIPSE_A, 2) + pow(ry, 2) / pow(ELLIPSE_B, 2) + pow(rz, 2) / pow(ELLIPSE_C, 2); 
				if (re2 > 1.0) {
					pol[j][m].state = MT_STATE_DEPOL;
				}
			}
			
		}
	}
	

}

void kineticStep(){
	int m,j;
	for (j = 0; j < 2; j++){
		float p = ran2(&rseed);
		for(m = 0; m < MT_NUM; m++){
			
			if(pol[j][m].state == MT_STATE_POL){
				if (p <= MT_CATASTROPHE_FREQ){
					//printf("catastrophe!\n");
					pol[j][m].state = MT_STATE_DEPOL;
				}
			} else

			if(pol[j][m].state == MT_STATE_DEPOL){
				if (p <= MT_RESCUE_FREQ){
					pol[j][m].state = MT_STATE_POL;
				}
			}

			if(pol[j][m].state == MT_STATE_POL && pol[j][m].n + 1 < MAX_MT_LENGTH){ // Warning needed
				p = ran2(&rseed);
				if (p <= MT_POL_RATE){
					int n = pol[j][m].n;
					float dx = pol[j][m].r[n-1].x - pol[j][m].r[n-2].x;
					float dy = pol[j][m].r[n-1].y - pol[j][m].r[n-2].y;
					float dz = pol[j][m].r[n-1].z - pol[j][m].r[n-2].z;
					float dr = sqrt(dx*dx + dy*dy + dz*dz);
					dx = 2.0*MT_RADIUS*dx/dr;
					dy = 2.0*MT_RADIUS*dy/dr;
					dz = 2.0*MT_RADIUS*dz/dr;
					pol[j][m].r[n].x = pol[j][m].r[n-1].x + dx;
					pol[j][m].r[n].y = pol[j][m].r[n-1].y + dy;
					pol[j][m].r[n].z = pol[j][m].r[n-1].z + dz;
					pol[j][m].n ++;
				}
			}
			if(pol[j][m].state == MT_STATE_DEPOL && pol[j][m].n > 2){ // Warning needed
				p = ran2(&rseed);
				if (p <= MT_DEPOL_RATE){
					//printf("Depol!\n");
					int n = pol[j][m].n;
					pol[j][m].r[n-1].x = 0.0;
					pol[j][m].r[n-1].y = 0.0;
					pol[j][m].r[n-1].z = 0.0;
					pol[j][m].n --;
				}
			}
		
		}
	}
	
}

void membraneMechanics(float3* r, float3* f){
	int m, i, mi, j;
	float re2, norm;
	for (j = 0; j < 2; j++){
		for(m = 0; m < MT_NUM; m++){
			for(mi = 0; mi < pol[j][m].n-1; mi++){

				i = j * MT_NUM * MAX_MT_LENGTH + m * MAX_MT_LENGTH + mi;
				float3 ri = r[i];
				re2 = pow(ri.x, 2) / pow(ELLIPSE_A, 2) + pow(ri.y, 2) / pow(ELLIPSE_B, 2) + pow(ri.z, 2) / pow(ELLIPSE_C, 2);
				norm = sqrt(pow(ri.x, 2) + pow(ri.y, 2) + pow(ri.z, 2)); 

				if (re2 > 1.0){
					f[i].x -= MEMBRANE_STIFF * (re2 - 1.0) * ri.x / norm;
					f[i].y -= MEMBRANE_STIFF * (re2 - 1.0) * ri.y / norm;		
					f[i].z -= MEMBRANE_STIFF * (re2 - 1.0) * ri.z / norm; 
				}
				
			}
		}
	}
}

void computeHarmonic(float3* r, float3* f){
	int i, j, m, mi, k;

	float r0 = 2.0*MT_RADIUS;
	for (k = 0; k < 2; k++){

		for(m = 0; m < MT_NUM; m++){
			for(mi = 0; mi < pol[k][m].n-1; mi++){
				i = k * MT_NUM * MAX_MT_LENGTH + m * MAX_MT_LENGTH + mi;
				j = i + 1;

				float3 ri = r[i];
				float3 rj = r[j];

				float dx = rj.x - ri.x;
				float dy = rj.y - ri.y;
				float dz = rj.z - ri.z;
		
				float dr = sqrtf(dx*dx + dy*dy + dz*dz);
				float df = K_SPRING*(dr - r0)/dr;
			
				f[i].x += df*dx;
				f[i].y += df*dy;		
				f[i].z += df*dz;
				f[j].x -= df*dx;
				f[j].y -= df*dy;
				f[j].z -= df*dz;
			}
		}
	}
	
}

void computeAngles(float3* r, float3* f){
	int i, j, k, m, mi, l;

	float theta0 = M_PI;
	for (l = 0; l < 2; l++){
		for(m = 0; m < MT_NUM; m++){
			for(mi = 0; mi < pol[l][m].n-2; mi++){
				i = l * MT_NUM * MAX_MT_LENGTH + m * MAX_MT_LENGTH + mi;
				j = i + 1;
				k = i + 2;
			
				float3 r1 = r[i];
				float3 r2 = r[j];
				float3 r3 = r[k];

				float3 dr12, dr32;
				dr12.x = r1.x - r2.x;
				dr12.y = r1.y - r2.y;
				dr12.z = r1.z - r2.z;

				dr32.x = r3.x - r2.x;
				dr32.y = r3.y - r2.y;
				dr32.z = r3.z - r2.z;

				float r12inv = 1.0f/sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);
				float r32inv = 1.0f/sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);
				float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if(costheta > 1.0f){
					costheta = 1.0f;
				} else
				if(costheta < -1.0f){
					costheta = -1.0f;
				}
				float sintheta = sqrtf(1.0f - costheta*costheta);
				float theta = acos(costheta);
				float diff = theta - theta0;
				if(sintheta < 1.e-6){
					if(diff < 0){
						diff *= 2.0f*K_THETA;
					} else {
						diff *= -2.0f*K_THETA;
					}
				} else {
					diff *= (-2.0f*K_THETA) / sintheta;
				}
				float c1 = diff*r12inv;
				float c2 = diff*r32inv;

				float3 f1, f2, f3;
				f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
				f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
				f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
				f2 = f1;
				f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
				f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
				f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
				f2.x += f3.x;
				f2.y += f3.y;
				f2.z += f3.z;

				f2.x = -f2.x;
				f2.y = -f2.y;
				f2.z = -f2.z;

				f[i].x += f1.x;
				f[i].y += f1.y;
				f[i].z += f1.z;

				f[j].x += f2.x;
				f[j].y += f2.y;
				f[j].z += f2.z;

				f[k].x += f3.x;
				f[k].y += f3.y;
				f[k].z += f3.z;
			}
		}
	}
	
}

void integrateCPU(float3* r, float3* f, float dt, int N){

	int i, m, mi, j;

	float gamma = computeGamma(MT_RADIUS);

	double var = sqrtf(KB*T*2.0*dt/gamma);
	for (j = 0; j < 2; j++){
		for(m = 0; m < MT_NUM; m++){
			for(mi = 1; mi < pol[j][m].n; mi++){
				i = j * MT_NUM * MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi;
		
				r[i].x += f[i].x*dt/gamma + var*gasdev(&rseed);
				r[i].y += f[i].y*dt/gamma + var*gasdev(&rseed);
				r[i].z += f[i].z*dt/gamma + var*gasdev(&rseed);	

				f[i].x = 0.0f;
				f[i].y = 0.0f;
				f[i].z = 0.0f;
			}
		}
	}
	

}


int writeXYZ(const char* filename, float3* data, int N, const char* modifier){
	FILE* file = fopen(filename, modifier);
	fprintf(file, "%d\n", N);
	fprintf(file, "POLE_COORDS\n");
	int i;
	for (i = 0; i < N; i++){
		fprintf(file, "MT\t%f\t%f\t%f\n", data[i].x, data[i].y, data[i].z);
	}
	fclose(file);
}


int initPSF(){
	PSF psf;
	psf.natom = mds.N;
	psf.atoms = (PSFAtom*)calloc(psf.natom, sizeof(PSFAtom));
	psf.nbond = MT_NUM*(MAX_MT_LENGTH-1);
	psf.bonds = (PSFBond*)calloc(psf.nbond, sizeof(PSFBond));
	psf.ntheta = 0;
	psf.nphi = 0;
	psf.nimphi = 0;
	psf.ncmap = 0;
	psf.nnb = 0;
	int i, j, m, mi, b, l;
	for(i = 0; i < psf.natom; i++){
		psf.atoms[i].id = i + 1;
		if(i < mds.N){
			sprintf(psf.atoms[i].name, "MT");
			sprintf(psf.atoms[i].type, "MT");
			sprintf(psf.atoms[i].segment, "MT");
		}
		sprintf(psf.atoms[i].resName, "Pole");
		psf.atoms[i].resid = i + 1;
		psf.atoms[i].m = 1.0;
		psf.atoms[i].q = 0.0;
	}

	float r0 = 2.0*MT_RADIUS;

	b = 0;
	/*
	for (l = 0; l < 2 ; l++){
		for(m = 0; m < MT_NUM; m++){
			for(mi = 0; mi < MAX_MT_LENGTH-1; mi++){

				psf.bonds[b].i = l * MT_NUM * MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi + 1;
				psf.bonds[b].j = l * MT_NUM * MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi + 2;
				b++;
			}
		}

	}
	*/
	char filename[1024];

	sprintf(filename, "output/MT.psf");
	writePSF(filename, &psf);
}
