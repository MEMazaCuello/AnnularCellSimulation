/****
  'montecarloRoutines.h':
  --------------------
  Implementation of 'montecarloRoutines'.
  For declaration details, see 'montecarloRoutines.cpp'.
  
  Needs: 'AnnularCell.h'
  
  --------------------
  The 'montecarloRoutines' define how a Monte Carlo step
  is done in the 'AnnularCell' representing the system.
  It also include functions that compute and save the 
  local order parameters -- nematic, tetratic and smectic, 
  of the liquid crystal of 'Rod's inside the 'AnnularCell'. 
        
  --------------------
  Last modified: 2019-06-02 
  By: M. E. Maza Cuello
****/

#ifndef MONTECARLOROUTINES_H_INCLUDED
#define MONTECARLOROUTINES_H_INCLUDED

#include "AnnularCell.h"

/* Methods */

void stepMontecarlo(AnnularCell& cell);
  /**
    'stepMontecarlo': Perform a Monte Carlo step on the 'cell'.
                      I. e. for each 'Rod' inside the 'cell', try
                      to randomly make a (small) change in its 
                      position and orientation.
                      The changes' amplitude is dynamically adjusted
                      so that the acceptation probability is ~ 50 %.
  **/

void stepMontecarloNeighbourhood(AnnularCell& cell);
  /**
    'stepMontecarloNeighbourhood': Perform a Monte Carlo step on the 'cell'.
                                   I. e. for each grid neighborhood inside the 'cell', try
                                   to randomly make a (small) change in its
                                   'Rod's' position and orientation.
                                   The changes' amplitudes are dynamically adjusted
                                   so that the acceptation probability is ~ 50 %.
  **/

void getOrderParameters(AnnularCell& cell, std::string& filename);
  /**
    'getOrderParameters': Compute the local nematic, tetratic and smectic 
                          order parameters of the liquid crystal of 'Rod's. 
                          Save the configuration and the order parameters
                          in a .txt file called 'filename.txt'.
  **/

#endif // MONTECARLOROUTINES_H_INCLUDED
