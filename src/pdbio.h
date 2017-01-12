/*
 * pdbio.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#ifndef PDBIO_H_
#define PDBIO_H_
#include <stdio.h>

/*
 * Structures
 */

typedef struct {

  int id;
  char   name[5], chain, resName[4], altLoc;
  int    resid;
  double x, y, z;

  char segment[5];

  double occupancy;
  double beta;

} PDBAtom;

typedef struct {

	int resid1;
	char chain1;

	int resid2;
	char chain2;

} PDBSSBond;

typedef struct {
	int atomCount;
	PDBAtom* atoms;
	int ssCount;
	PDBSSBond* ssbonds;
} PDB;

/*
 * Public methods
 */
void readPDB(const char* filename, PDB* pdbData);
void writePDB(const char* filename, PDB* pdbData);
void printAtom(PDBAtom atomData);
void printAtomToFile(FILE* file, PDBAtom atomData);


#endif /* PDBIO_H_ */
