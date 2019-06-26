/****
  'montecarloRoutines.h':
  --------------------
  Declaration of 'montecarloRoutines'.
  For implementation details, see 'montecarloRoutines.h'.
  
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

#include "montecarloRoutines.h"
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

/** PARAMETERS OBTAINED FROM 'parameters.cpp' **/
// Number of Rods
extern const int    NUMBER_OF_RODS;
// Lengths
extern const double LENGTH;
extern const double HALF_DIAGONAL;
extern const double OUTER_RADIUS;
// Angles
extern const double PI;
extern const double HALF_PI;

/* Methods */

// Random number generator
std::mt19937_64 g(98710536);
std::uniform_real_distribution<double> dist(-1.0,1.0);

void stepMontecarlo(AnnularCell& cell)
{
  /**
    'stepMontecarlo': Perform a Monte Carlo step on the 'cell'.
                      I. e. for each 'Rod' inside the 'cell', try
                      to randomly make a (small) change in its 
                      position and orientation.
                      The changes' amplitudes are dynamically adjusted
                      so that the acceptation probability is ~ 50 %.
  **/
  
  // Auxiliary constant variables
  static const double MAX_RADIUS = OUTER_RADIUS-HALF_DIAGONAL;
  static const double INV_NUMBER_OF_RODS = 1.0d/NUMBER_OF_RODS;
  // Auxiliary variables:
  static double DELTA_SPACE = 0.1*LENGTH;   // Maximum variation of coordinates
  static double DELTA_ANGLE = 0.1*HALF_PI;  // Maximum variation of orientation
  // Auxiliary array of indexes of the 'Rod's
  static std::vector<int> indexes(NUMBER_OF_RODS);
  
  // Generate random index permutation
  std::iota(indexes.begin(), indexes.end(), 0);
  std::shuffle(indexes.begin(), indexes.end(), g);
  
  // Boolean used to exit loops and mark new position as valid
  bool validPosition;
  // Number of 'Rod's that have been succesfully moved
  int  numSuccess = 0;
  for(int i : indexes)
  {
    // Get 'Rod' with index 'i'
    cell.m_aux_rod = cell.getRod(i);

    /* Add random displacement */
    // Coordinates (max displacement = DELTA_SPACE)
    cell.m_aux_rod.m_xPos  += dist(g)*DELTA_SPACE;
    cell.m_aux_rod.m_yPos  += dist(g)*DELTA_SPACE;
    // Orientation (max displacement = DELTA_ANGLE)
    cell.m_aux_rod.m_angle += dist(g)*DELTA_ANGLE;
    
    /* Adjust positions to fit inside the 'AnnularCell' */
    // Coordinates in interval [-OUTER_RADIUS, OUTER_RADIUS]
    if(cell.m_aux_rod.m_xPos > OUTER_RADIUS){ cell.m_aux_rod.m_xPos = MAX_RADIUS;}
    else if(cell.m_aux_rod.m_xPos < - OUTER_RADIUS){ cell.m_aux_rod.m_xPos = -MAX_RADIUS;}
    if(cell.m_aux_rod.m_yPos > OUTER_RADIUS){ cell.m_aux_rod.m_yPos = MAX_RADIUS;}
    else if(cell.m_aux_rod.m_yPos < - OUTER_RADIUS){ cell.m_aux_rod.m_yPos = -MAX_RADIUS;}
    // Orientation in [-HALF_PI, HALF_PI]
    while(cell.m_aux_rod.m_angle > HALF_PI)
    {
        cell.m_aux_rod.m_angle -= PI;
    }
    while(cell.m_aux_rod.m_angle < -HALF_PI)
    {
        cell.m_aux_rod.m_angle += PI;
    }

    // Assume the new position is valid
    validPosition = true;
    // If 'm_aux_rod' is touching inner or outer walls of 'cell',
    if(cell.rodIsTouchingInnerWall(cell.m_aux_rod)||cell.rodIsTouchingOuterWall(cell.m_aux_rod))
    {
      // ... new position is not valid
      validPosition = false;
    }else{
      // If not, get neighbours of 'm_aux_rod'      
      cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_aux_rod));
      for(int idx : cell.m_grid.m_neighbours)
      { 
        if(idx != i)
        {
          // ... check if is touching any neighbour...
          if(cell.m_aux_rod.isTouchingRod(cell.getRod(idx)))
          {
            // ... and if it is, new position is not valid 
            validPosition = false;
            break;
          }
        }
      }
    }

    // If position is valid
    if(validPosition)
    {
      // ... move index to new coordinates in the 'm_grid'
      cell.m_grid.moveIndex(i,cell.m_grid.getGridCoords(cell.m_bundle[i]), cell.m_grid.getGridCoords(cell.m_aux_rod));
      // ... save 'm_aux_rod' in the corresponding index
      cell.m_bundle[i]= cell.m_aux_rod;
      // ... and add one success
      numSuccess++;
    }
  }

  /* Correct maximum displacements to generate ~ 50 % acceptance */
  DELTA_SPACE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);
  DELTA_ANGLE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);
}

void stepMontecarloNeighbourhood(AnnularCell& cell)
{
    /**
        'stepMontecarloNeighbourhood': Perform a Monte Carlo step on the 'cell'.
                                        I. e. for each grid neighborhood inside the 'cell', try
                                        to randomly make a (small) change in its
                                        'Rod's' position and orientation.
                                        The changes' amplitudes are dynamically adjusted
                                        so that the acceptation probability is ~ 50 %.
    **/

    // Auxiliary constant variables
    static const double MAX_RADIUS = OUTER_RADIUS-HALF_DIAGONAL;
    static const double INV_NUMBER_OF_RODS = 1.0d/NUMBER_OF_RODS;
    // Auxiliary variables:
    static double DELTA_SPACE = 0.1*LENGTH;   // Maximum variation of coordinates
    static double DELTA_ANGLE = 0.1*HALF_PI;  // Maximum variation of orientation

    // Auxiliary array of indexes of the 'Rod's
    std::vector<int> indexes;
    // Boolean used to exit loops and mark new position as valid
    bool validPosition;
    // Number of 'Rod's that have been successfully moved
    int  numSuccess = 0;

    static std::vector<int> gridBoxes(cell.m_grid.m_BOXES*cell.m_grid.m_BOXES);

    // Generate random index permutation
    std::iota(gridBoxes.begin(), gridBoxes.end(), 0);
    std::shuffle(gridBoxes.begin(), gridBoxes.end(), g);

   int gi, gj;
   for(int gbox : gridBoxes)
   {
    gi = gbox % cell.m_grid.m_BOXES;
    gj = gbox / cell.m_grid.m_BOXES;

    indexes = cell.m_grid.getBox(gi, gj);
    cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(gi, gj);
    for(int i : indexes)
    {
      // Get 'Rod' with index 'i'
      cell.m_aux_rod = cell.getRod(i);

      /* Add random displacement */
      // Coordinates (max displacement = DELTA_SPACE)
      cell.m_aux_rod.m_xPos  += dist(g)*DELTA_SPACE;
      cell.m_aux_rod.m_yPos  += dist(g)*DELTA_SPACE;
      // Orientation (max displacement = DELTA_ANGLE)
      cell.m_aux_rod.m_angle += dist(g)*DELTA_ANGLE;

      /* Adjust positions to fit inside the 'AnnularCell' */
      // Coordinates in interval [-OUTER_RADIUS, OUTER_RADIUS]
      if(cell.m_aux_rod.m_xPos > OUTER_RADIUS){ cell.m_aux_rod.m_xPos = MAX_RADIUS;}
      else if(cell.m_aux_rod.m_xPos < - OUTER_RADIUS){ cell.m_aux_rod.m_xPos = -MAX_RADIUS;}
      if(cell.m_aux_rod.m_yPos > OUTER_RADIUS){ cell.m_aux_rod.m_yPos = MAX_RADIUS;}
      else if(cell.m_aux_rod.m_yPos < - OUTER_RADIUS){ cell.m_aux_rod.m_yPos = -MAX_RADIUS;}
      // Orientation in [-HALF_PI, HALF_PI]
      while(cell.m_aux_rod.m_angle > HALF_PI)
      {
        cell.m_aux_rod.m_angle -= PI;
      }
      while(cell.m_aux_rod.m_angle < -HALF_PI)
      {
        cell.m_aux_rod.m_angle += PI;
      }

      // Assume the new position is valid
      validPosition = true;
      // If 'm_aux_rod' is touching inner or outer walls of 'cell',
      if(cell.rodIsTouchingInnerWall(cell.m_aux_rod)||cell.rodIsTouchingOuterWall(cell.m_aux_rod))
      {
        // ... new position is not valid
        validPosition = false;
      }else{
        // If not, for each neighbour
        for(int idx : cell.m_grid.m_neighbours)
        {
          if(idx != i)
          {
            // ... check if is touching any neighbour...
            if(cell.m_aux_rod.isTouchingRod(cell.getRod(idx)))
            {
              // ... and if it is, new position is not valid
              validPosition = false;
              break;
            }
          }
        }
      }

      // If position is valid
      if(validPosition)
      {
        // ... move index to new coordinates in the 'm_grid'
        cell.m_grid.moveIndex(i,cell.m_grid.getGridCoords(cell.m_bundle[i]), cell.m_grid.getGridCoords(cell.m_aux_rod));
        // ... save 'm_aux_rod' in the corresponding index
        cell.m_bundle[i]= cell.m_aux_rod;
        // ... and add one success
        numSuccess++;
      }
    }
  }

  /* Correct maximum displacements to generate ~ 50 % acceptance */
  DELTA_SPACE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);
  DELTA_ANGLE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);
}

void getOrderParameters(AnnularCell& cell, std::string& filename){
  /**
    'getOrderParameters': Compute the local nematic, tetratic and smectic 
                          order parameters of the liquid crystal of 'Rod's. 
                          Save the configuration and the order parameters
                          in a .txt file called 'filename.txt'.
  **/
  
  // Auxiliar constant variable
  static const double scale = 2.0d*PI/(1.2d*LENGTH);
  
  // File stream
  std::ofstream Qdat;
  Qdat.open(filename);

  /* Auxiliar angular variables */ 
  // Angles
  double tilt: // Local nematic director
  double zeta; // Double of relative orientation with respect to tilt 
  // Sines and cosines for order parameters
  double sum_cos2, sum_sin2; // Sum of cos(2*angle) and sin(2*angle)
  double sum_cos4, sum_sin4; // Sum of cos(4*angle) and sin(4*angle)
  double sum_cosS, sum_sinS; // Sum of cos(scale*h) and sin(scale*h), 
                             // were 'scale' is inverse distance between smectic layers 
                             // and 'h' is relative position with respect to the local layer

  for(int i=0; i<NUMBER_OF_RODS; i++){
    // Initialize sums to zero
    sum_cos2 = 0.0d;
    sum_sin2 = 0.0d;
    
    /* Compute local nematic director */
    // Get neighbours of 'Rod' of index 'i'
    cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_bundle[i]));
    for(int index : cell.m_grid.m_neighbours)
    {
      // NOTE: index INCLUDES current 'i' 'Rod'
      sum_cos2 += std::cos(2.0d*cell.m_bundle[index].m_angle);
      sum_sin2 += std::sin(2.0d*cell.m_bundle[index].m_angle);  
    }
    // Local nematic director
    tilt = 0.5d*std::atan2(sum_sin2,sum_cos2);

    // Initialize sums to zero
    sum_cos2 = 0.0d; sum_sin2 = 0.0d;
    sum_cos4 = 0.0d; sum_sin4 = 0.0d;
    sum_cosS = 0.0d; sum_sinS = 0.0d;

    /* Extract order parameters */
    for(int index : cell.m_grid.m_neighbours)
    {
      // Relative angle between 'Rod' orientation and local nematic director
      zeta = 2.0d*(cell.m_bundle[index].m_angle-tilt);
      
      // Sum for local nematic order parameter
      sum_cos2 += std::cos(zeta);
      sum_sin2 += std::sin(zeta);
      // Sum for local tetratic order parameter
      sum_cos4 += std::cos(2.0d*zeta);
      sum_sin4 += std::sin(2.0d*zeta);

      // Relative position between 'Rod' and local smectic layers
      zeta = scale*(std::cos(tilt)*(cell.m_bundle[index].m_xPos-cell.m_bundle[i].m_xPos) + std::sin(tilt)*(cell.m_bundle[index].m_yPos-cell.m_bundle[i].m_yPos));
      
      // Sum for local smectic order parameter
      sum_cosS += std::cos(zeta);
      sum_sinS += std::sin(zeta);
    }

    // Inverse number of neighbours
    zeta = 1.0d/cell.m_grid.m_neighbours.size();
    
    /* Save local information */
    Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " " // Index  and Position (x,y)
         << cell.m_bundle[i].m_angle << " " << tilt << " "                // Orientation and local tilt
         << std::sqrt(sum_cos2*sum_cos2 + sum_sin2*sum_sin2)*zeta << " "  // Nematic  order parameter Q1
         << std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)*zeta << " "  // Tetratic order parameter Q2
         << std::sqrt(sum_cosS*sum_cosS + sum_sinS*sum_sinS)*zeta <<      // Smectic  order parameter Qs
         '\r' << '\n'; // New line

  }
  
  // Close file stream
  Qdat.close();
}
