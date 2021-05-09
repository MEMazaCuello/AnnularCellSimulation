/**
  * "montecarloRoutines.h":
  * --------------------
  * Declaration of Montecarlo methods.
  * For implementation details, see "montecarloRoutines.cpp".
  *
  * Needs: "AnnularCell.h".
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#ifndef MONTECARLOROUTINES_H_INCLUDED
#define MONTECARLOROUTINES_H_INCLUDED

#include "AnnularCell.h"

/**
  * Performs single Montecarlo step (MCS) over the cell.
  * A MCS tries to translate and rotate each single particle by
  * a small amount. If no overlapping occurs, the new position and
  * orientation are saved for this particle.
  * The maximum amounts that a particle can translate or rotate are
  * dynamically controlled so as to obtain (in average) the acceptance
  * rate set in "parameters.cpp".
  * Convergence of acceptance rate usually takes some tens of iterations.
  */

// Rods are tried in random order (random permutation of indexes)
void stepMontecarlo(AnnularCell& cell);

// Rods within a grid box are tried in order (random permutation of boxes)
void stepMontecarloNeighborhood(AnnularCell& cell);

#endif // MONTECARLOROUTINES_H_INCLUDED
