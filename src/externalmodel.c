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
#include <julia.h>

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
  sprintf(modelname, "ACEhamiltonians");

  // flags for capabilities of this specific model
  *capabilities = (mycapabilities) {
    .hamiltonian = true, .overlap = true, .energy = false,
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
  
  
  /* Todo:
    - There is a non-trivial chance overflow errors arising from type collisions. A
      safeguard should be put inplace to ensure that integers provided by Fortran
      and Julia match up with those use by C; i.e. all 64 or 32 bit with no mixing.
    - Set TLS settings
    - The Julia model should be updated to:
      - identify itself as either a Hamiltonian or overlap model.
      - provide names for the species present.
      - yield neutral occupations.
    - Need to convert units from Bohr to Angstroms.

  */

  // Allocate structure for internal state, and generate am intptr for
  // return to DFTB+
  struct mystate* internalState = (struct mystate*)
    malloc(sizeof(struct mystate));
  *state = (intptr_t) internalState;



  FILE *input;
  int ii, jj, kk;

  // Initialise Julia
  jl_init();
  
  jl_value_t *h_model;
  jl_value_t *s_model;
  
  // Large structures will frequently get deleted from memory due to Julia's garbage
  // collection subroutine. Thus one must "lock" important Julia based objects in
  // memory by creating an artificial reference to prevent this from happening.
  jl_value_t* refs = jl_eval_string("refs = IdDict()");
  jl_function_t* setindex = jl_get_function(jl_base_module, "setindex!");

  jl_eval_string("include(\"/home/ajmhpc/Projects/ACEhamiltonians/Code/Working/DFTB_Interface/ACEhamiltonians.jl-Dev/tools/dftbp_api.jl\")");
  
  jl_function_t *load_model = jl_eval_string("load_model");
  jl_function_t *n_orbs_per_atom = jl_eval_string("n_orbs_per_atom");
  jl_function_t *offers_species = jl_eval_string("offers_species");
  jl_function_t *species_name_to_id = jl_eval_string("species_name_to_id");
  jl_function_t *max_interaction_cutoff = jl_eval_string("max_interaction_cutoff");
  jl_function_t *max_environment_cutoff = jl_eval_string("max_environment_cutoff");
  jl_function_t *shells_on_species = jl_eval_string("shells_on_species!");
  jl_function_t *n_shells_on_species = jl_eval_string("n_shells_on_species");
  jl_function_t *shell_occupancies = jl_eval_string("shell_occupancies!");


  jl_value_t* array_type_i32 = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
  jl_value_t* array_type_f64 = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);

  /* LOADING MODELS 
     Load in the ACEhamiltonians models for the Hamiltonian and overlap matrices
     from the files located at the paths specified by the environmental variables
     "H_MODEL_PATH" and "S_MODEL_PATH" respectively. These are then stored in the
     state entity as `hamiltonian_model` and `overlap_model`.
  */
  if(getenv("H_MODEL_PATH")) {
      h_model = jl_call1(load_model, jl_cstr_to_string(getenv("H_MODEL_PATH")));
      jl_call3(setindex, refs, h_model, h_model);
      internalState->hamiltonian_model = h_model;
    } else {
      // The "printf" command can be removed once DFTB+ has been updated to print out
      // the "message" string feed.
      printf("Environmental variable \"H_MODEL_PATH\" has not been set.\n");
      sprintf(message, "Environmental variable \"H_MODEL_PATH\" has not been set.\n");
      return -1;
  }
  

  if(getenv("S_MODEL_PATH")) {
      s_model = jl_call1(load_model, jl_cstr_to_string(getenv("S_MODEL_PATH")));
      jl_call3(setindex, refs, s_model, s_model);
      internalState->overlap_model = s_model;
    } else {
      // The "printf" command can be removed once DFTB+ has been updated to print out
      // the "message" string feed.
      printf("Environmental variable \"S_MODEL_PATH\" has not been set.\n");
      sprintf(message, "Environmental variable \"S_MODEL_PATH\" has not been set.\n");
      return -1;
  
  }

  *environmentCutoff = jl_unbox_float64(jl_call1(max_environment_cutoff, h_model));
  *interactionCutoff = jl_unbox_float64(jl_call1(max_interaction_cutoff, h_model));
  internalState->n_orbitals = calloc(*nspecies, sizeof(int));
  internalState->species_id = calloc(*nspecies, sizeof(int));
  *nShellsOnSpecies =  calloc(*nspecies, sizeof(int));

  int n_total_shells = 0;
  int max_n_shells = 0;

  bool offered;
  
  jl_value_t *name;
  jl_value_t *species_i_id;

  for (ii = 0; ii < *nspecies; ii++) {
    name = jl_cstr_to_string(speciesName[ii]);
    offered = jl_unbox_bool(jl_call2(offers_species, name, h_model));

    if (offered) {  // If the species is provided by the model
      // Get the id used by the model to refer to the species
      species_i_id = jl_call2(species_name_to_id, name, h_model);

      // Store the id into the species_id reference array
      internalState->species_id[ii] = jl_unbox_int32(species_i_id);

      // Add its orbital count to the n_orbitals array 
      internalState->n_orbitals[ii] = jl_unbox_int32(jl_call2(n_orbs_per_atom, species_i_id, h_model));

      // Work out the number of shells on the species
      (*nShellsOnSpecies)[ii] = jl_unbox_int32(jl_call2(n_shells_on_species, species_i_id, h_model));

      // Keep track of the i) total number of all shells and, ii) the maximum 
      // number of shells on any species for use later on.
      n_total_shells = n_total_shells + (*nShellsOnSpecies)[ii];

      if ((*nShellsOnSpecies)[ii] > max_n_shells) {
        max_n_shells = (*nShellsOnSpecies)[ii];
      }

    } else {  // If it is not offered by the model then raise an error.
      printf("Species %s is not offered by the provided model\n", speciesName[ii]);
      sprintf(message, "Species %s is not offered by the provided model\n", speciesName[ii]);
      return -3;
    };
  }

  // The arrays `shellLValues` and `shellOccs` specify the azimuthal quantum
  // number and occupancy for each shell respectively. While the size of the
  // first array, `shellLValues`, is what one might expect, the second array,
  // `shellOccs`, is typically larger. This is because it is reshaped by DFTB+
  // into n×m Fortran array where `n` is the number of species and `m` is the
  // maximum number of shells present on any atom. 
  *shellLValues = calloc(n_total_shells, sizeof(int));
  *shellOccs =  calloc(*nspecies * max_n_shells, sizeof(double));


  kk = 0;
  jl_array_t *shell_ls, *shell_ocs;
  
  for (ii=0; ii<*nspecies; ii++) {
    
    jj = ii * max_n_shells;

    // Gather the part of the array `shellLValues` that will hold the azimuthal quantum
    // numbers associated with species `ii`. 
    shell_ls = jl_ptr_to_array_1d(array_type_i32, (*shellLValues)+kk, (*nShellsOnSpecies)[ii], 0);
    shell_ocs = jl_ptr_to_array_1d(array_type_f64, (*shellOccs)+jj, (*nShellsOnSpecies)[ii], 0);

    // Fill the azimuthal number array using the Julia function `shells_on_species`.
    jl_call3(
      shells_on_species, (jl_value_t*)shell_ls,
      jl_box_int32((*internalState).species_id[ii]), h_model);

    // Fill the occupation array using the Julia function `shells_on_species`.
    jl_call3(
      shell_occupancies, (jl_value_t*)shell_ocs,
      jl_box_int32((*internalState).species_id[ii]), h_model);

    // Update the counter index which is responsible for tracking the start of the
    // next slice of `shellLValues`.
    kk = kk + (*nShellsOnSpecies)[ii];
  }


  (*internalState).initialised = true;

  (*internalState).nAtomicClusters = 0;
  (*internalState).indexAtomicClusters = NULL;
  (*internalState).atomicClusters = NULL;
  (*internalState).atomicGlobalAtNos = NULL;
  (*internalState).atomic_species_ids = NULL;

  (*internalState).nBndClusters = 0;
  (*internalState).indexBndClusters = NULL;
  (*internalState).bndClusters = NULL;
  (*internalState).bndGlobalAtNos = NULL;
  (*internalState).bondClusterIndex = NULL;
  (*internalState).bond_species_ids = NULL;
  

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;



}


// Update this model, using geometric and other information from DFTB+
int update_model_for_dftbp(intptr_t *state, int* species, int* nAtomicClusters,
                           int* indexAtomicClusters, double* atomicClusters,
                           int* atomicGlobalAtNos, int* nBndClusters,
                           int* indexBndClusters, double* bndClusters,
                           int* bndGlobalAtNos, int* atomClusterIndex,
                           int* bondClusterIndex, char* message)
{
  int ii, count;

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
  internalState->atomClusterIndex = atomClusterIndex;
  
  // Construct an array specifying the species-id of each and every atom present
  // in the atom-clusters. 
  count = indexAtomicClusters[*nAtomicClusters] - 1;
  internalState->atomic_species_ids = calloc(count, sizeof(int));
  for (ii = 0; ii < count; ii++) {
    internalState->atomic_species_ids[ii] = internalState->species_id[species[atomicGlobalAtNos[ii] - 1] - 1];
  }
  
  internalState->nBndClusters = *nBndClusters;
  internalState->indexBndClusters = indexBndClusters;
  internalState->bndClusters = bndClusters;
  internalState->bndGlobalAtNos = bndGlobalAtNos;
  internalState->bondClusterIndex = bondClusterIndex;


  // Construct an array specifying the species-id of each and every atom present
  // in the bond-clusters. 
  count = indexBndClusters[*nBndClusters] - 1;
  internalState->bond_species_ids = calloc(count, sizeof(int));
  for (ii = 0; ii < count; ii++) {
    internalState->bond_species_ids[ii] = internalState->species_id[species[bndGlobalAtNos[ii] - 1] - 1];
  }

  // blank return message if nothing happening
  sprintf(message, "\n");
  return 0;

}


// Make and then get model predictions back to DFTB+
int predict_model_for_dftbp(intptr_t *state, double *h0, double *over,
                            char* message)
{
  struct mystate* internalState = (struct mystate*) *state;

  jl_function_t *build_on_site_atom_block = jl_eval_string("build_on_site_atom_block!");
  jl_function_t *build_off_site_atom_block = jl_eval_string("build_off_site_atom_block!");

  jl_value_t* array_type_i32 = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
  jl_value_t* array_type_f64 = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);

  int i_block, i_start, j_start, i_species, j_species, n_orbs_i, n_orbs_j, n_neighbours;
  int i_atom, j_atom;
  
  jl_value_t *args[4];

  // Loop over all on site atom blocks by index
  for (i_block=0; i_block<(*internalState).nAtomicClusters; i_block++) {
      
      // Index specifying the start of the atom-block `i` within the Hamiltonian
      // and overlap matrix arrays `h0` and `over` respectively.
      i_start = (*internalState).atomClusterIndex[i_block];

      // Index specifying the start of the cluster associated with atom-block `i`
      // in arrays `atomicClusters` & `species_id`. Must less by one to account
      // for indices here using Fortran indexing convention (unlike atomClusterIndex).
      j_start = (*internalState).indexAtomicClusters[i_block] - 1;

      // Index of the atom to which this atom-block pertains
      i_atom = (*internalState).atomicGlobalAtNos[j_start] - 1;

      // Species index; this is just 1 less than the number as used internally by DFTB+
      i_species = (*internalState).globalSpeciesOfAtoms[i_atom] - 1;

      // Number of orbitals on the associated species
      n_orbs_i = internalState->n_orbitals[i_species];

      // Number of atoms present in cluster `i`
      n_neighbours = (*internalState).indexAtomicClusters[i_block+1] - 1 - j_start;

      // Gather the coordinates associate with atom-block `i` 
      args[1] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64,(*internalState).atomicClusters + (3 * j_start),
        3 * n_neighbours, 0);

      // Get the species id's of the atoms present in the atom cluster 
      args[2] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_i32,
        (*internalState).atomic_species_ids + j_start,
        n_neighbours,
        0);

      // Collect the relevant atom-block from the Hamiltonian matrix
      args[0] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64, h0+i_start, n_orbs_i * n_orbs_i, 0);
    
      // Add the Hamiltonian model to the arguments list
      args[3] = (jl_value_t*)(*internalState).hamiltonian_model;

      // Call out to the Julia on-site block constructor method to build the
      // current atom-block of the Hamiltonian matrix      
      jl_call(build_on_site_atom_block, args, 4);

      // Repeat this process for the associated block of the overlap matrix
      args[0] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64, over+i_start, n_orbs_i * n_orbs_i, 0);
      args[3] = (jl_value_t*)(*internalState).overlap_model;
      jl_call(build_on_site_atom_block, args, 4);
    }

  // Loop over all on site atom blocks by index
  for (i_block=0; i_block<(*internalState).nBndClusters; i_block++) {

      // Index specifying the start of the atom-block `i` within the Hamiltonian
      // and overlap matrix arrays `h0` and `over` respectively.
      i_start = (*internalState).bondClusterIndex[i_block];

      // Index specifying the start of the cluster associated with atom-block `i`
      // in arrays `bndClusters` & `species_id`. Must less by one to account
      // for indices here using Fortran indexing convention (unlike bondClusterIndex).
      j_start = (*internalState).indexBndClusters[i_block] - 1;

      // Index of the atoms to which this interaction-block pertains
      i_atom = (*internalState).bndGlobalAtNos[j_start] - 1;
      j_atom = (*internalState).bndGlobalAtNos[j_start + 1] - 1;

      // Species index; this is just one less than the number as used
      // internally by DFTB+
      i_species = (*internalState).globalSpeciesOfAtoms[i_atom] - 1;
      j_species = (*internalState).globalSpeciesOfAtoms[j_atom] - 1;
      
      // Number of orbitals on the associated species
      n_orbs_i = internalState->n_orbitals[i_species];
      n_orbs_j = internalState->n_orbitals[j_species];

      // Number of atoms present in cluster `i`
      n_neighbours = (*internalState).indexBndClusters[i_block+1] - 1 - j_start;

      // Gather the coordinates associate with atom-block `i` 
      args[1] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64,(*internalState).bndClusters + (3 * j_start),
        3 * n_neighbours, 0);

      // Get the species id's of the atoms present in the atom cluster 
      args[2] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_i32,
        (*internalState).bond_species_ids + j_start,
        n_neighbours,
        0);

      // Collect the relevant atom-block from the Hamiltonian matrix
      args[0] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64, h0+i_start, n_orbs_i * n_orbs_j, 0);
    
      // Add the Hamiltonian model to the arguments list
      args[3] = (jl_value_t*)(*internalState).hamiltonian_model;

      // Call out to the Julia on-site block constructor method to build the
      // current atom-block of the Hamiltonian matrix      
      jl_call(build_off_site_atom_block, args, 4);

      // Repeat this process for the associated block of the overlap matrix
      args[0] = (jl_value_t*)jl_ptr_to_array_1d(
        array_type_f64, over+i_start, n_orbs_i * n_orbs_j, 0);
      args[3] = (jl_value_t*)(*internalState).overlap_model;
      jl_call(build_off_site_atom_block, args, 4);
      
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
  jl_atexit_hook(0);

}
