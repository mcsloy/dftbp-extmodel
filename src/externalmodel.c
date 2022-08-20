/*------------------------------------------------------------------------------------------------*/
/*  Example external model for the DFTB+ project: www.dftbplus.org                                */
/*  Copyright (C) 2006 - 2022  B. Hourahine                                                       */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    /* Declared capabilities of this model for the DFTB+ code */

    // Can this library deliver a hamiltonian?
    bool hamiltonian;

    // Can this library deliver an overlap matrix (i.e. non-orthogonal)?
    bool overlap;

    // Does this library deliver energy terms other than the band-strucuture energy?
    bool energy;

    // Can model results be returned for a subset of atoms?
    bool atomsubset;

    // Is this library MPI aware (so requires a communicator)?
    bool mpi;

  } mycapabilities;


  void dftbp_provided_with(typeof (mycapabilities) *capabilities){

    /* Declare capabilities of this model to DFTB+ via the API */

    *capabilities = (mycapabilities) {
      .hamiltonian=false, .overlap = true, .energy = false, .atomsubset = false, .mpi = false
    };

  }

#ifdef __cplusplus
}
#endif
