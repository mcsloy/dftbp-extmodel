/**

@brief
Example external model for the DFTB+ project: www.dftbplus.org
Copyright (C) 2006 - 2022  B. Hourahine

See the LICENSE file for terms of usage and distribution.

@file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

  /** Declared capabilities of this model for initial setup of the
      DFTB+ code */
  typedef struct {

    /** Can this library deliver a hamiltonian to the calling code? */
    bool hamiltonian;

    /** Can this library deliver an overlap matrix (i.e. is a non-orthogonal model)? */
    bool overlap;

    /** Does this library deliver energy terms (other than the band-structure energy)? */
    bool energy;

    /** Is the model self-consistent */
    bool selfconsistent;

    /** Number of spin channels supported by model (none:0, collinear:1, non-collinear:3) */
    int spinchannels;

    /** Can model results be returned for a subset of atoms? */
    bool atomsubset;

    /** Is this library MPI aware (so requires a communicator)? */
    bool mpi;

  } mycapabilities;


  /** Internal state of this model, including any initialisation
      passed in from DFTB+ */
  typedef struct {

    
  } mystate;


  /**
     Combined application programming interface and application binary
     interface semantic version, as supported by this library. Note,
     the format and is externally defined in the DFTB+ project

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
     Declare capabilities of this model to DFTB+ via the external model API.

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
     Set up this model, read some settings from DFTB+ over it's
     external model API, read a parameter file and initialise it's
     data structure for handling via DFTB+.

     @param nspecies number of chemical species/types present
     @param species array of null terminated strings labelling chemical species
     @param message return message, in event of routine failure (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

   */
  int initialise_model_for_dftbp(int* nspecies, char* species[], char* message) {

    FILE *input;
    int i;

    /* Open input file for some constants for this model, assuming
       it's in the runtime directory */
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
    }

    // blank return message if nothing happening
    sprintf(message, "\n");
    return 0;

  }


#ifdef __cplusplus
}
#endif
