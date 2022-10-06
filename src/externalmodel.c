/**

@brief
Example external model for the DFTB+ project: www.dftbplus.org
Copyright (C) 2022 B. Hourahine

See the LICENSE file for terms of usage and distribution.

@file
*/

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/externalmodel.h"

// Project version
void dftbp_model_apbi(int* major, int* minor, int* patch)
{
  *major = 0;
  *minor = 1;
  *patch = 0;
}


// Declare capabilities of this model to DFTB+
void dftbp_provided_with(char* modelname,
                         typeof (mycapabilities) *capabilities)
{

  // Name of external model
  sprintf(modelname, "Huckel toy model");

  // flags for capabilities of this specific model
  *capabilities = (mycapabilities) {
    .hamiltonian = true, .overlap = false, .energy = false,
    .derivativeOrder = 0, .selfconsistent = false, .spinchannels = 0
  };

  return;
}


// Initialise this model, reading some settings from DFTB+
int initialise_model_for_dftbp(int* nspecies, char* speciesName[],
                               double* interactionCutoff,
                               double* environmentCutoff,
                               int** nShellsOnSpecies,
                               int** shellLValues, double** shellOccs,
                               intptr_t *state, char* message)
{

  // Allocate structure for internal state, and generate am intptr for
  // return to DFTB+
  struct mystate* internalState = (struct mystate*)
    malloc(sizeof(struct mystate));
  *state = (intptr_t) internalState;

  FILE *input;
  int ii, items;

  // Open input file for some constants for this model, assuming it's
  // in the runtime directory
  input = fopen("input.dat", "r");
  if (!input) {
    sprintf(message, "Library error opening input file.\n");
    return -1;
  }

  // read ancillary input file for model parameters and then store
  // them into the model's internal structure, that in turn will get passed
  // around between calls to the model.
  items = fscanf(input, "%lf %lf %lf", &internalState->onsites[0],
                 &internalState->onsites[1], &internalState->onsites[2]);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file at first line\n");
    return -3;
  }
  if (items != 3) {
    sprintf(message, "Toy library malformed first line of data file: %i\n",
	    items);
    return -3;
  }
  items = fscanf(input, "%lf %lf %lf", &internalState->hopping[0],
		 &internalState->hopping[1], &internalState->hopping[2]);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 2nd line"
	    " (hoping integrals)\n");
    return -3;
  }
  if (items != 3) {
    sprintf(message, "Toy library malformed second line of data file\n");
    return -3;
  }
  items = fscanf(input, "%lf", interactionCutoff);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 3rd line"
	    " (bond cut-off)\n");
    return -3;
  }
  if (items != 1) {
    sprintf(message, "Toy library malformed third line of data file\n");
    return -3;
  }

  for (ii = 0; ii < *nspecies; ii++) {
    // This specific model is only for H and C atoms, so will throw an
    // error otherwise
    if (strcmp(speciesName[ii], "C") != 0 && strcmp(speciesName[ii], "H") != 0)
      {
	sprintf(message,
		"Toy library only knows about C and H atoms, not %s.\n",
		speciesName[ii]);
	return -2;
      }
    if (strcmp(speciesName[ii], "H") == 0) {
      (*internalState).species2params[ii] = 0;
    }
    if (strcmp(speciesName[ii], "C") == 0) {
      (*internalState).species2params[ii] = 1;
    }
  }

  // Surroundings required around atoms and bonds:
  *environmentCutoff = 4.0;

  // Maximum number of shells of orbitals on atoms in this model
  int maxShells = 2;

  /* This particular model is so only a single s shell on H and two s
     shells on C Note count for shells uses Fortran convention
     starting at 1, while L uses physics convention of s=0 */
  *nShellsOnSpecies =  calloc(*nspecies, sizeof(int));
  // Note, this will be a fortran array on the other end, so
  // [nSpecies][maxShells] (hence the 2):
  *shellLValues =  calloc(*nspecies * maxShells, sizeof(int));
  for (ii=0; ii<*nspecies; ii++) {
    if ((*internalState).species2params[ii] == 0)
      {
        // Hydrogen, s
        (*nShellsOnSpecies)[ii] = 1;
        (*shellLValues)[ii] = 0;
      } else
      {
        // Carbon, s and s* shells
        (*nShellsOnSpecies)[ii] = 2;
        (*shellLValues)[ii] = 0;
        (*shellLValues)[ii+1] = 0;
      }
  }

  // each atom is neutral when it's single shell containing only one
  // electron:
  *shellOccs =  calloc(*nspecies * maxShells, sizeof(double));
  for (ii=0; ii<*nspecies; ii++) {
    int jj = ii * maxShells;
    if ((*internalState).species2params[ii] == 0)
      {
        // H atom, the s orbital is occupied
        (*shellOccs)[jj] = 1.0;
      } else {
      // C atom, only the s orbital is occupied
      (*shellOccs)[jj] = 1.0;
    }
  }

  (*internalState).initialised = true;

  printf(" Initial on-site energies : H %f, C %f %f\n",
         (*internalState).onsites[0], (*internalState).onsites[1],
	 (*internalState).onsites[2]);


  (*internalState).nAtomicClusters = 0;
  (*internalState).indexAtomicClusters = NULL;
  (*internalState).atomicClusters = NULL;
  (*internalState).atomicGlobalAtNos = NULL;

  /*
  (*internalState).nBndClusters = 0;
  (*internalState).indexBndClusters = NULL;
  (*internalState).bndClusters = NULL;
  (*internalState).bndGlobalAtNos = NULL;
  */

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


// Update this model, using geometric and other information from DFTB+
int update_model_for_dftbp(intptr_t *state, int* species, int* nAtomicClusters,
                           int* indexAtomicClusters, double* atomicClusters,
                           int* atomicGlobalAtNos, int* nBndClusters,
                           int* indexBndClusters, double* bndClusters,
                           int* bndGlobalAtNos, char* message)
{

  // map pointer back to structure
  struct mystate* internalState = (struct mystate*) *state;

  if (!(*internalState).initialised) {
    sprintf(message, "Model is not properly initialised");
    return -1;
  }


  internalState->globalSpeciesOfAtoms = species;

  internalState->nAtomicClusters = *nAtomicClusters;
  internalState->indexAtomicClusters = indexAtomicClusters;
  internalState->atomicClusters = atomicClusters;
  internalState->atomicGlobalAtNos = atomicGlobalAtNos;

  printf("Number of atomic clusters: %i\n", *nAtomicClusters);

  internalState->nBndClusters = *nBndClusters;
  internalState->indexBndClusters = indexBndClusters;
  internalState->bndClusters = bndClusters;
  internalState->bndGlobalAtNos = bndGlobalAtNos;

  printf("Number of bond clusters: %i\n", *nBndClusters);

  printf(" Initial on-site energies : H %f, C %f %f\n",
         (*internalState).onsites[0], (*internalState).onsites[1],
	 (*internalState).onsites[2]);

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


// Get model predictions back to DFTB+
int predict_model_for_dftbp(intptr_t *state, double *h0, int* h0Index,
                            int* h0IndexStride, int* nElemPerAtom,
                            char* message)
{

  int jj;

  struct mystate* internalState = (struct mystate*) *state;

  // On-site matrix elements for each atom
  for (int ii=0; ii<(*internalState).nAtomicClusters; ii++) {
    jj = ii * *h0IndexStride; // stride on h0Index 2D array
    // Start of on-site block for atom ii in H0 matrix:
    int iStart = h0Index[jj];
    int iAtom = (*internalState).atomicGlobalAtNos[jj];
    int iSpecies = (*internalState).globalSpeciesOfAtoms[iAtom-1];

    // simple use of onsite energies
    if ((*internalState).species2params[iSpecies-1] == 0)
      {
        h0[iStart] = *((*internalState).onsites);
      } else
      {
        // 2x2 symmetric matrix of s and s* onsites, reshaped as a 4
        // element vector:
        h0[iStart] = *((*internalState).onsites+1); // s orbital
        h0[iStart+3] = *((*internalState).onsites+2); // s* orbital
      }
  }

  // Toy crystal field model, demonstrating getting atoms in the halo
  // around each atomic site
  //
  // Model affects s and s* mixing on C atoms, but only from any
  // surrounding H atoms within the environment cutoff distance
  // (ignoring C). H atoms themselves would be unaffected by their
  // surroundings, as there is only a single s orbital on those atoms.
  for (int ii = 0; ii < (*internalState).nAtomicClusters; ii++) {
    jj = ii * *h0IndexStride; // stride on h0Index 2D array
    // Start of on-site block for atom ii in H0 matrix:
    int iAtom = (*internalState).atomicGlobalAtNos[jj];
    int iSpecies = (*internalState).globalSpeciesOfAtoms[iAtom-1];
    int iParam = (*internalState).species2params[iSpecies-1];
    if (iParam == 1) { // C atom at origin of cluster
      int iStart = h0Index[jj];
      // Carbon atom, so check it's neighbours
      int jStart = (*internalState).indexAtomicClusters[ii] - 1;
      int jEnd = (*internalState).indexAtomicClusters[ii+1] - 1;
      for (int iAt = jStart; iAt < jEnd; iAt++) {
        int kk = 3 * iAt;
        int ll = (*internalState).atomicGlobalAtNos[iAt];
        int mm = (*internalState).globalSpeciesOfAtoms[ll-1];
        int jParam = (*internalState).species2params[mm-1];
        if (jParam == 0) { // H atom in suroundings
          // distance squared
          double d2 =
            (*((*internalState).atomicClusters+kk))
	    *(*((*internalState).atomicClusters+kk))
            +(*((*internalState).atomicClusters+kk+1))
	    *(*((*internalState).atomicClusters+kk+1))
            +(*((*internalState).atomicClusters+kk+2))
	    *(*((*internalState).atomicClusters+kk+2));
	  // scale such that final term is ~1E-12 by 4.0 a.u. cutoff
          d2 *= -1.0*(6.0/4.0)*(6.0/4.0);
          h0[iStart+1] += 1000.0*exp(d2);
          h0[iStart+2] += 1000.0*exp(d2);
        }
      }
    }
  }



  // off-site (diatomic) elements
  for (int iClust = 0; iClust < (*internalState).nBndClusters; iClust++) {
    int jAtStart = (*internalState).indexBndClusters[iClust] - 1;
    int jAtEnd = (*internalState).indexBndClusters[iClust+1] - 1;

    //printf("Cluster %i, %i:%i\n", iClust+1, jAtStart, jAtEnd-1);

    int jOrigAt = (*internalState).bndGlobalAtNos[jAtStart] - 1;

    int iHamBnd = jOrigAt * *h0IndexStride + 1; // stride on h0Index
						// 2D array

    printf("iHam neigh of %i: %i, %i\n", jOrigAt, iHamBnd,
	   h0Index[iHamBnd]-1);
    printf("Bonds from atom %i : %i\n", jOrigAt+1, nElemPerAtom[jOrigAt]);

    // coordinates
    for (int iAt = jAtStart; iAt < jAtEnd; iAt++) {
    //  jj = iClust * *h0IndexStride + ;
    //  printf("Atom in bond cluster %i %i %i\n", iClust, iAt, jj);
    //  int iStart = h0Index[jj];
    //  printf("%i\n", iStart);
    //
    //  int ll = h0Index[jj+kk];
    //  h0[ll] += 0.1;
    }
  }

  // blank return message if nothing failing
  sprintf(message, "\n");
  return 0;

}


// Clean up after this model, freeing any memory in the mystate type
void cleanup_model_for_dftbp(intptr_t *state) {

  // DFTB+ only sees integer pointer "state", so need original
  // structure to clean up
  struct mystate* internalState = (struct mystate*) *state;

  printf("Cleaning up\n");

  free(internalState);
  *state = (intptr_t) internalState;

}
