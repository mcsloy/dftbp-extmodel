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
#include "../api/externalmodel.h"


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
int initialise_model_for_dftbp(int* nspecies, char* species[], double cutoffs[], typeof (mystate) *state,
                               char* message) {

  FILE *input;
  int i;

  /* Open input file for some constants for this model, assuming it's
     in the runtime directory */
  input=fopen("input.dat", "r");
  if (!input) {
    sprintf(message, "Library error opening input file.\n");
    return -1;
  }

  // This specific model is only for H and C atoms, so will throw an error otherwise
  for (i=0;i<*nspecies;i++) {
    if (strcmp(species[i], "C") != 0 && strcmp(species[i], "H") != 0) {
      sprintf(message, "Toy library only knows about C and H atoms, not %s.\n", species[i]);
      return -2;
    }

    if (strcmp(species[i], "C") == 0) {
      cutoffs[i] = 2.35;
    }
    if (strcmp(species[i], "H") == 0) {
      cutoffs[i] = 1.35;
    }

  }

  *state = (mystate) {.initialised = true, .number = 6};

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
  printf("Internal model state cargo %d\n", (*state).number);

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
