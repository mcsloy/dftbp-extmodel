/*------------------------------------------------------------------------------------------------*/
/*  Example external model for the DFTB+ project: www.dftbplus.org                                */
/*  Copyright (C) 2006 - 2022  B. Hourahine                                                       */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    /* Declared capabilities of this model for initial setup of the
       DFTB+ code */

    // Can this library deliver a hamiltonian?
    bool hamiltonian;

    // Can this library deliver an overlap matrix (i.e. is a non-orthogonal model)?
    bool overlap;

    // Does this library deliver energy terms (other than the band-structure energy)?
    bool energy;

    // Is the model self-consistent
    bool selfconsistent;

    // Can model results be returned for a subset of atoms?
    bool atomsubset;

    // Is this library MPI aware (so requires a communicator)?
    bool mpi;

  } mycapabilities;


  typedef struct {
    /* Internal state of this model, including any initialisation from DFTB+ */

  } mystate;


  void dftbp_provided_with(char* modelname, typeof (mycapabilities) *capabilities){

    /* Declare capabilities of this model to DFTB+ via the API */

    // Name of external model
    sprintf(modelname, "Huckel toy model");

    // flags for capabilities
    *capabilities = (mycapabilities) {
      .hamiltonian = true, .overlap = false, .energy = false, .selfconsistent = false,
      .atomsubset = false, .mpi = false
    };

    return;
  }


  int initialise_model_for_dftbp(int* nspecies, char* species[], char* message) {
    /* Set up model, read some settings from DFTB+, read a parameter file */

    FILE *input;
    int i;

    /* Open input file for some constants for this model, assuming
       it's in the runtime directory */
    input=fopen("input.dat", "r");
    if (!input) {
      sprintf(message, "Library error opening input file.\n");
      return -1;
    }

    for (i=0;i<*nspecies;i++) {
      if (strcmp(species[i], "C") != 0 && strcmp(species[i], "H") != 0) {
	sprintf(message, "Toy library only knows about C and H atoms, not %s.\n", species[i]);
	return -2;
      }
    }

    // blank message if nothing happening
    sprintf(message, "\n");
    return 0;

  }

#ifdef __cplusplus
}
#endif
