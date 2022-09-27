/**

@brief
Example external model for the DFTB+ project: www.dftbplus.org
Copyright (C) 2022 B. Hourahine

See the LICENSE file for terms of usage and distribution.

@file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "../include/externalmodel.h"


/**
   Combined application programming interface and application binary
   interface semantic version, as supported by this library. Note, the
   format and is externally defined in the DFTB+ project

   @param major Major version, revised on breaking changes
   @param minor Minor version, revised on extensions
   @param patch Patch version, revised on invisible changes

*/
void dftbp_model_apbi(int* major, int* minor, int* patch){
  *major = 0;
  *minor = 1;
  *patch = 0;
}


/**
   Declare capabilities of this model to DFTB+ via the external model
   API.

   @param modelname null terminated string for name of this model
   @param capabilities structure with capabilities of the model

*/
void dftbp_provided_with(char* modelname, typeof (mycapabilities) *capabilities){

  // Name of external model
  sprintf(modelname, "Huckel toy model");

  // flags for capabilities of this specific model
  *capabilities = (mycapabilities) {
    .hamiltonian = true, .overlap = false, .energy = false, .derivativeOrder = 0,
    .selfconsistent = false, .spinchannels = 0
  };

  return;
}

/**
   Set up this model, read some settings from DFTB+ over it's external
   model API and initialise it's data structure for handling via
   DFTB+.

   @param nspecies number of chemical species/types present
   @param species array of null terminated strings labelling chemical species
   @param interactionCutoff array of cutoffs for distance over which
   atoms of each species have interactions (i.e. diatomic matrix
   elements or repulsive contributions)
   @param environmentCutoff Distance over which neighbours influence
   interactions, i,e. 0 if model is environmentally independent, or
   nearest-neighbour/longer if surrounding atoms influence
   interactions. This is used to cut out local clusters of around
   the interacting atom/dimer
   @param nShellsOnSpecies number of shells of atomic orbitals, set to 0 if not
   a hamiltonian model
   @param shellLValues Angular momentum of shells species resolved atomic
   shells, freed on return to DFTB+
   @param shellOccs Reference occupation for neutrality, freed on return to DFTB+
   @param state internal state and data of the model, not checked in
   DFTB+, just passed around
   @param message return message, in event of routine failure (return != 0)

   @return 0 on successful return, non-zero if there is an error
   message to check

*/
int initialise_model_for_dftbp(int* nspecies, char* speciesName[], double* interactionCutoff,
                               double* environmentCutoff, int** nShellsOnSpecies,
                               int** shellLValues, double** shellOccs, intptr_t *state,
                               char* message) {

  // Allocate structure for internal state, and generate intptr for return to DFTB+
  typeof(mystate)* internalState = (typeof(mystate)*) malloc(sizeof(typeof(mystate)));
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
  items = fscanf(input, "%lf %lf", &internalState->onsites[0], &internalState->onsites[1] );
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file at first line\n");
    return -3;
  }
  if (items != 2) {
    sprintf(message, "Toy library malformed first line of data file: %i\n", items);
    return -3;
  }
  items = fscanf(input, "%lf %lf %lf", &internalState->hopping[0], &internalState->hopping[1],
                 &internalState->hopping[2]);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 2nd line (hoping integrals)\n");
    return -3;
  }
  if (items != 3) {
    sprintf(message, "Toy library malformed second line of data file\n");
    return -3;
  }
  items = fscanf(input, "%lf", interactionCutoff);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 3rd line (bond cut-off)\n");
    return -3;
  }
  if (items != 1) {
    sprintf(message, "Toy library malformed third line of data file\n");
    return -3;
  }

  for (ii = 0; ii < *nspecies; ii++) {
    // This specific model is only for H and C atoms, so will throw an error otherwise
    if (strcmp(speciesName[ii], "C") != 0 && strcmp(speciesName[ii], "H") != 0) {
      sprintf(message, "Toy library only knows about C and H atoms, not %s.\n",
              speciesName[ii]);
      return -2;
    }
    (*internalState).species2params[ii] = -1;
    if (strcmp(speciesName[ii], "H") == 0) {
      (*internalState).species2params[ii] = 0;
    }
    if (strcmp(speciesName[ii], "C") == 0) {
      (*internalState).species2params[ii] = 1;
    }
  }

  /* Surroundings required around atoms and bonds: */
  *environmentCutoff = 4.0;

  /* This particular model is Huckel-like, so only a single shells on
     each species containing an s-like orbital, irrespective of
     species. Note count for shells uses Fortran convention starting
     at 1, while L uses physics convention of s=0 */
  *nShellsOnSpecies =  malloc(*nspecies * sizeof(int));
  *shellLValues =  malloc(*nspecies * 1 * sizeof(int));
  for (ii=0; ii<*nspecies; ii++) {
    (*nShellsOnSpecies)[ii] = 1;
    (*shellLValues)[ii] = 0;
  }

  // each atom is neutral when it's single shell containing only one
  // electron:
  *shellOccs =  malloc(*nspecies * 1 * sizeof(double));
  for (ii=0; ii<*nspecies; ii++) {
    if (strcmp(speciesName[ii], "C") == 0) {
      (*shellOccs)[ii] = 1.0;
    }
    if (strcmp(speciesName[ii], "H") == 0) {
      (*shellOccs)[ii] = 1.0;
    }
  }

  (*internalState).initialised = true;

  printf(" Initial on-site energies : H %f, C %f\n", (*internalState).onsites[0],
         (*internalState).onsites[1]);

  (*internalState).nAtomClusters = 0;
  (*internalState).indexAtomicClusters = NULL;
  (*internalState).atomicClusters = NULL;

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


/**
   Update this model, using geometric and other information from DFTB+
   over it's external model API.

   @param state internal state and data of the model, this is not
   checked by DFTB+, just passed around by it

   @param species Species index for atoms in the global structure

   @param nAtomicClusters Number of atom centred clusters

   @param indexAtomicClusters starting index for location of
   coordinates

   @param atomicClusters Geometric clusters centred on atoms for the
   onsite matrix element predictions

   @param clusterGlobalAtNos Numbers of atoms in clusters in the
   global system

   @param message return message, in event of routine failure
   (return != 0)

   @return 0 on successful return, non-zero if there is an error
   message to check

*/
int update_model_for_dftbp(intptr_t *state, int* species, int* nAtomicClusters,
                           int* indexAtomicClusters, double* atomicClusters,
                           int* clusterGlobalAtNos, char* message) {

  // map pointer back to structure
  typeof(mystate)* internalState = (typeof(mystate)*) *state;

  if (!(*internalState).initialised) {
    sprintf(message, "Model is not properly initialised");
    return -1;
  }

  internalState->globalSpeciesOfAtoms = species;

  internalState->nAtomClusters = *nAtomicClusters;
  internalState->indexAtomicClusters = indexAtomicClusters;
  internalState->atomicClusters = atomicClusters;
  internalState->clusterGlobalAtNos = clusterGlobalAtNos;

  printf("Number of atomic clusters: %i\n", *nAtomicClusters);


  /*
  for (int ii=0; ii<*nAtomicClusters; ii++) {
    printf("Atom index %i %i:%i\n", ii+1, indexAtomicClusters[ii],
           indexAtomicClusters[ii+1]-1);
    for (int jj=indexAtomicClusters[ii]-1; jj<indexAtomicClusters[ii+1]-1; jj++) {
      int kk = 3*jj;
      printf("%f %f %f\n", 0.529177249*atomicClusters[kk],
             0.529177249*atomicClusters[kk+1], 0.529177249*atomicClusters[kk+2]);
    }
    }
  */

  printf("On-site energies internally : H %f, C %f\n", (*internalState).onsites[0],
         (*internalState).onsites[1]);

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


/**
    Get model predictions

    @param state internal state and data of the model, this is not
    checke by DFTB+, just passed around by it

    @param h0 hamiltonian

    @param h0Index hamiltonian index for blocks of matrix elements

    @param message return message, in event of routine failure
    (return != 0)

    @return 0 on successful return, non-zero if there is an error
    message to check

*/
int predict_model_for_dftbp(intptr_t *state, double* h0, int* h0Index, int* h0IndexStride,
                              int* nElemPerAtom, char* message) {

  int jj;

  typeof(mystate)* internalState = (typeof(mystate)*) *state;

  //printf("On-site energies : H %f, C %f\n", (*internalState).onsites[0],
  //       (*internalState).onsites[1]);

  /*
  for (int ii=0; ii<(*internalState).nAtomClusters; ii++) {
    // Print out cluster of atoms surrounding each atom
    iStart = *((*internalState).indexAtomicClusters+ii);
    iEnd = *((*internalState).indexAtomicClusters+ii+1)-1;
    printf("Atom index %i %i:%i\n", ii+1, iStart, iEnd-1);
    for (int jj=iStart-1; jj<iEnd; jj++) {
      int kk = 3*jj;
      int ll = *((*internalState).clusterGlobalAtNos+jj);
      int mm = *((*internalState).globalSpeciesOfAtoms+ll-1);
      printf("%i %i %f %f %f\n", ll, mm,
	     0.529177249 * *((*internalState).atomicClusters+kk),
             0.529177249 * *((*internalState).atomicClusters+kk+1),
             0.529177249 * *((*internalState).atomicClusters+kk+2));
    }
  }
  */

  // On-site matrix elements
  for (int ii=0; ii<(*internalState).nAtomClusters; ii++) {
    jj = ii * *h0IndexStride; // stride on h0Index 2D array
    // Start of on-site block for atom ii in H0 matrix:
    int iStart = h0Index[jj];
    int iAtom = *((*internalState).clusterGlobalAtNos+ii);
    int iSpecies = *((*internalState).globalSpeciesOfAtoms+iAtom-1);
    int iParam = (*internalState).species2params[iSpecies-1];
    printf("Index %i %i %i %i %i\n", ii, jj, iStart, iSpecies, iParam);
    h0[iStart] = (*internalState).onsites[iParam];

    /*
    int iAtStart = *((*internalState).indexAtomicClusters+ii);
    int iAtEnd = *((*internalState).indexAtomicClusters+ii+1)-1;
    // Loop over surounding atoms in cluster, and add an (unphysical) shift to onsite energies
    for (int jj=iAtStart; jj<iAtEnd; jj++) {
      int kk = 3*jj;
      int ll = *((*internalState).clusterGlobalAtNos+jj);
      int mm = *((*internalState).globalSpeciesOfAtoms+ll-1);
      printf("%i %i %f %f %f\n", ll, mm,
	     0.529177249 * *((*internalState).atomicClusters+kk),
             0.529177249 * *((*internalState).atomicClusters+kk+1),
             0.529177249 * *((*internalState).atomicClusters+kk+2));
    }
    */
  }


  // blank return message if nothing failing
  sprintf(message, "\n");
  return 0;

}


/**
   Clean up after this model, freeing any memory in the mystate type

   @param state internal state and data of the model. This is not
   checke by DFTB+, just passed around by it, so we need to remove any
   allocated memory here.

   @param message return message, in event of routine failure
   (return != 0)

   @return 0 on successful return, non-zero if there is an error
   message to check

*/
int cleanup_model_for_dftbp(intptr_t *state, char* message) {

  // DFTB+ only sees integer pointer "state", so need original
  // structure to clean up
  typeof(mystate)* internalState = (typeof(mystate)*) *state;

  printf("\nInternal check for cleanup_model_for_dftbp, Model is initialised? ");
  printf((*internalState).initialised ? "true\n" : "false\n");
  printf("Cleaning up\n");

  free(internalState);
  *state = (intptr_t) internalState;

  (*internalState).initialised = false;

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}
