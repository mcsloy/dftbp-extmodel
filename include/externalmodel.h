/**

   @brief
   Headers for example external model for the DFTB+ project:
   www.dftbplus.org Copyright (C) 2022 B. Hourahine

   See the LICENSE file for terms of usage and distribution.

   @file
*/

#ifndef H_EXTERNALMODEL
#define H_EXTERNALMODEL

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

    /** Order of derivatives returned by the model */
    int derivativeOrder;

    /** Is the model self-consistent */
    bool selfconsistent;

    /** Number of spin channels supported by model (none:0, collinear:1, non-collinear:3) */
    int spinchannels;

  } mycapabilities;


  /** Internal state of this model, including any initialisation passed
      in from DFTB+, any data read in, checks for state, etc. */
  typedef struct {

    bool initialised;

    // internal model parameters and state

    float onsites[2]; // H and C alph
    float hopping[3]; // H-H, H-C and C-C beta
    float cutoffs[3]; // H-H, H-C and C-C cutoff distances, in a.u.

  } mystate;


  /**
     Combined application programming interface and application binary
     interface semantic version, as supported by this library. Note,
     the format and is externally defined in the DFTB+ project

     @param major Major version, revised on breaking changes
     @param minor Minor version, revised on extensions
     @param patch Patch version, revised on invisible changes

  */
  void dftbp_model_apbi(int* major, int* minor, int* patch);

  /**
     Declare capabilities of this model to DFTB+ via the external model API.

     @param modelname null terminated string for name of this model
     @param capabilities structure with capabilities of the model

  */
  void dftbp_provided_with(char* modelname, typeof (mycapabilities) *capabilities);

  /**
     Set up this model, read some settings from DFTB+ over it's
     external model API and initialise it's data structure for
     handling via DFTB+.

     @param nspecies number of chemical species/types present
     @param species array of null terminated strings labelling chemical species
     @param interactCutoff Longest cutoff for, i.e. distance over which
     atoms of each species have hamiltonian or repulsive interactions
     @param environmentCutoff Distance over which neighbours influence
     interactions, i,e. 0 if model is environmentally independent, or
     nearest-neighbour/longer if surrounding atoms influence
     interactions. This is used to cut out local clusters of around
     the interacting atom/dimer
     @param nShellsOnSpecies number of shells of atomic orbitals, set to 0 if not
     a hamiltonian model
     @param shells Angular momentum of shells species resolved atomic
     shells, freed on return to DFTB+
     @param shellOccs Reference occupation for neutrality, freed on return to DFTB+
     @param state internal state and data of the model, not checked by
     DFTB+, just passed around
     @param message return message, in event of routine failure (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

   */
  int initialise_model_for_dftbp(int* nspecies, char* species[], double* interactCutoff,
                                 double* environmentCutoff, int* nShellsOnSpecies[], int** shells,
                                 double** shellOccs, typeof(mystate) *state, char* message);


    /**
     Update this model, using geometric and other information from
     DFTB+ over it's external model API.

     @param state internal state and data of the model, this is not
     checke by DFTB+, just passed around by it

     @param message return message, in event of routine failure
     (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

   */
  int update_model_for_dftbp(typeof(mystate) *state, char* message);


  /**
      Get model predictions

      @param state internal state and data of the model, this is not
      checke by DFTB+, just passed around by it

      @param message return message, in event of routine failure
      (return != 0)

      @return 0 on successful return, non-zero if there is an error
      message to check

  */
  int predict_model_for_dftbp(typeof(mystate) *state, char* message);


  /**
     Clean up after this model, freeing any memory in the mystate type

     @param state internal state and data of the model. This is not
     checke by DFTB+, just passed around by it, so we need to remove
     any allocated memory here.

     @param message return message, in event of routine failure
     (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

  */
  int cleanup_model_for_dftbp(typeof (mystate) *state, char* message);
  
#ifdef __cplusplus
}
#endif

#endif
