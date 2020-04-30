/**
  * "montecarloRoutines.cpp":
  * --------------------
  * Implementation of Montecarlo methods.
  * For implementation details, see "montecarloRoutines.h".
  *
  * Needs: "montecarloRoutines.h".
  *
  * --------------------
  * Last modified: 2020-04-30
  * By: M. E. Maza-Cuello
  */

#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

#include "montecarloRoutines.h"

// Parameters defined in "parameters.cpp" _____________________________
extern const double LENGTH;
extern const double HALF_PI;
extern const double PI;
extern const double OUTER_RADIUS;
extern const double HALF_DIAGONAL;
extern const double ACCEPTANCE;
extern const int    SEED;

// Methods ____________________________________________________________

// Random number generator
std::mt19937_64 g(SEED);
std::uniform_real_distribution<double> dist(-1.0,1.0);

void stepMontecarlo(AnnularCell& cell)
{
  static const double INV_NUMBER_OF_RODS = 1.0d/NUMBER_OF_RODS;
  static double DELTA_SPACE = 0.1*LENGTH;   // Maximum variation of coordinates
  static double DELTA_ANGLE = 0.1*HALF_PI;  // Maximum variation of orientation
  static std::vector<int> indexes(NUMBER_OF_RODS);

  bool validPosition;
  int  numSuccess = 0;
  std::vector<int> oldCoords;
  std::vector<int> newCoords;

  // Generate random index permutation
  std::iota(indexes.begin(), indexes.end(), 0);
  std::shuffle(indexes.begin(), indexes.end(), g);

  for(int i : indexes)
  {
    cell.m_aux_rod = cell.getRod(i);

    // Old grid coordinates (before moving)
    oldCoords = cell.m_grid.getCoords(cell.m_aux_rod);

    cell.m_aux_rod.translate(dist(g)*DELTA_SPACE,dist(g)*DELTA_SPACE);
    cell.m_aux_rod.rotate(dist(g)*DELTA_ANGLE); // Also ensures that rod's angle is in [-pi/2,pi/2]

    if ( cell.rodIsOutsideWalls(cell.m_aux_rod) )
    {
      continue;
    }

    // New grid coordinates (after moving)
    newCoords = cell.m_grid.getCoords(cell.m_aux_rod);

    cell.m_grid.m_neighbors = cell.m_grid.getNeighbors(newCoords);

    validPosition = true;
    for (int idx : cell.m_grid.m_neighbors)
    {
      if (idx != i)  // Rod is neighbor of itself
      {
        if (cell.m_aux_rod.isTouchingRod(cell.getRod(idx)))
        {
          validPosition = false;
          break;
        }
      }
    }

    if (validPosition)
    {
      // Update grid coordinates of rod
      cell.m_grid.moveIndex(i, oldCoords, newCoords);
      cell.m_bundle[i]= cell.m_aux_rod;
      numSuccess++;
    }
  }

  // Correction of deltas to give (on average) acceptance rate
  DELTA_SPACE *= (1.0d - ACCEPTANCE + double(numSuccess)*INV_NUMBER_OF_RODS);
  DELTA_ANGLE *= (1.0d - ACCEPTANCE + double(numSuccess)*INV_NUMBER_OF_RODS);
}

void stepMontecarloNeighborhood(AnnularCell& cell)
{
  static const double INV_NUMBER_OF_RODS = 1.0d/NUMBER_OF_RODS;
  static double DELTA_SPACE = 0.1*LENGTH;   // Maximum variation of coordinates
  static double DELTA_ANGLE = 0.1*HALF_PI;  // Maximum variation of orientation
  static std::vector<int> gridBoxes(cell.m_grid.m_BOXES*cell.m_grid.m_BOXES);

  std::vector<int> indexes;
  bool validPosition;
  int  numSuccess = 0;
  std::vector<int> oldCoords;
  std::vector<int> newCoords;

  // Random boxes permutation
  std::iota(gridBoxes.begin(), gridBoxes.end(), 0);
  std::shuffle(gridBoxes.begin(), gridBoxes.end(), g);

  int gi, gj;
  for (int gbox : gridBoxes)
  {
    gi = gbox % cell.m_grid.m_BOXES;
    gj = gbox / cell.m_grid.m_BOXES;

    cell.m_grid.m_neighbors = cell.m_grid.getNeighbors(gi, gj);

    indexes = cell.m_grid.getBox(gi, gj);
    for (int i : indexes)
    {
      cell.m_aux_rod = cell.getRod(i);

      // Old grid coordinates (before moving)
      oldCoords = cell.m_grid.getCoords(cell.m_aux_rod);

      cell.m_aux_rod.translate(dist(g)*DELTA_SPACE,dist(g)*DELTA_SPACE);
      cell.m_aux_rod.rotate(dist(g)*DELTA_ANGLE); // Also ensures that rod's angle is in [-pi/2,pi/2]

      if ( cell.rodIsOutsideWalls(cell.m_aux_rod) )
      {
        continue;
      }

      // New grid coordinates (after moving)
      newCoords = cell.m_grid.getCoords(cell.m_aux_rod);

      validPosition = true;
      for (int idx : cell.m_grid.m_neighbors)
      {
        if (idx != i)  // Rod is neighbor of itself
        {
          if (cell.m_aux_rod.isTouchingRod(cell.getRod(idx)))
          {
            validPosition = false;
            break;
          }
        }
      }

      if (validPosition)
      {
        // Update grid coordinates of rod
        cell.m_grid.moveIndex(i, oldCoords, newCoords);
        cell.m_bundle[i]= cell.m_aux_rod;
        numSuccess++;
      }
    }
  }

  // Correction of deltas to give (on average) acceptance rate
  DELTA_SPACE *= (1.0d - ACCEPTANCE + double(numSuccess)*INV_NUMBER_OF_RODS);
  DELTA_ANGLE *= (1.0d - ACCEPTANCE + double(numSuccess)*INV_NUMBER_OF_RODS);
}
