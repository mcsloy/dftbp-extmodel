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
    .hamiltonian = true, .overlap = false, .energy = false, .selfconsistent = false,
    .spinchannels = 0, .atomsubset = false, .mpi = false
  };

  return;
}

/**
   Set up this model, read some settings from DFTB+ over it's external
   model API and initialise it's data structure for handling via
   DFTB+.

   @param nspecies number of chemical species/types present
   @param species array of null terminated strings labelling chemical species
   @param cutoffs array of cutoffs for distance over which atoms of
   each species have hamiltonian interactions
   @param state internal state and data of the model, not checked in
   DFTB+, just passed around
   @param message return message, in event of routine failure (return != 0)

   @return 0 on successful return, non-zero if there is an error
   message to check

*/
int initialise_model_for_dftbp(int* nspecies, char* species[], double* qmCutoff, typeof (mystate) *state,
                               char* message) {

  FILE *input;
  int i, items, natspec;

  /* Open input file for some constants for this model, assuming it's
     in the runtime directory */
  input=fopen("input.dat", "r");
  if (!input) {
    sprintf(message, "Library error opening input file.\n");
    return -1;
  }

  /* read ancillary input file for model parameters and then store
     them into the model's internal structure */
  items = fscanf(input,"%f %f", &state->onsites[0], &state->onsites[1] );
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file at first line\n");
    return -3;
  }
  if (items != 2) {
    sprintf(message, "Toy library malformed first line of data file: %i\n", items);
    return -3;
  }
  items = fscanf(input,"%f %f %f", &state->hopping[0], &state->hopping[1], &state->hopping[2]);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 2nd line\n");
    return -3;
  }
  if (items != 3) {
    sprintf(message, "Toy library malformed second line of data file\n");
    return -3;
  }
  items = fscanf(input,"%f %f %f", &state->cutoffs[0], &state->cutoffs[1], &state->cutoffs[2]);
  if (items == EOF) {
    sprintf(message, "Toy library malformed end of data file before 3rd line\n");
    return -3;
  }
  if (items != 3) {
    sprintf(message, "Toy library malformed third line of data file\n");
    return -3;
  }

  // This specific model is only for H and C atoms, so will throw an error otherwise
  *qmCutoff = 0.0;
  natspec = 0;
  for (i=0;i<*nspecies;i++) {
    if (strcmp(species[i], "C") != 0 && strcmp(species[i], "H") != 0) {
      sprintf(message, "Toy library only knows about C and H atoms, not %s.\n", species[i]);
      return -2;
    }
    if (strcmp(species[i], "H") == 0) {
      if (*qmCutoff < (*state).cutoffs[0]) {
        *qmCutoff = (*state).cutoffs[0];
      }
      natspec++;
    }
    if (strcmp(species[i], "C") == 0) {
      if (*qmCutoff < (*state).cutoffs[1]) {
        *qmCutoff = (*state).cutoffs[1];
      }
      natspec++;
    }
  }
  if (natspec == 2) {
    if (*qmCutoff < (*state).cutoffs[2]) {
      *qmCutoff = (*state).cutoffs[2];
    }
  }

  printf("internal cut %f", *qmCutoff);

  (*state).initialised = true;

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


/**
   Update this model, using geometric and other information from DFTB+
   over it's external model API.

   @param state internal state and data of the model, this is not
   checke by DFTB+, just passed around by it

   @param message return message, in event of routine failure
   (return != 0)

   @return 0 on successful return, non-zero if there is an error
   message to check

*/
int update_model_for_dftbp(typeof (mystate) *state, char* message) {

  printf("\nInternal check, Model is initialised? ");
  printf((*state).initialised ? "true\n" : "false\n");

  printf("On-site energies : H %f, C %f\n", (*state).onsites[0], (*state).onsites[1]);

  // blank return message if nothing happening
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
int cleanup_model_for_dftbp(typeof (mystate) *state, char* message) {

  printf("\nInternal check, Model is initialised? ");
  printf((*state).initialised ? "true\n" : "false\n");

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}
