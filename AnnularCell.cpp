/****
  'AnnularCell.cpp':
  --------------------
  Implementation of 'AnnularCell' class.
  For declaration details, see 'AnnularCell.cpp'.

  Needs: 'Rod.h'
         'Grid.h'

  --------------------
  An 'AnnularCell' is the area contained within
  two concentric circles defined by the
    'm_INNER_RADIUS' and the
    'm_OUTER_RADIUS',
  and contains 'm_NUMBER_OF_RODS'
  in a 'm_bundle' of 'Rod's.
  The 'Rod's are numbered with an index
  saved in a 'm_grid' to know its
  approximate position inside the 'AnnularCell'.

  Note: the file 'parameters.cpp' is needed to specify the parameters
        'INNER_RADIUS' and 'OUTER_RADIUS',
        'ALPHA', 'HALF_PI' and 'PI'.

  --------------------
  Last modified: 2019-11-10
  By: M. E. Maza Cuello
****/

#include "AnnularCell.h"
#include "montecarloRoutines.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

/** PARAMETERS OBTAINED FROM 'parameters.cpp' **/
// Lengths
extern const double INNER_RADIUS;
extern const double OUTER_RADIUS;
// Angles
extern const double HALF_PI;
extern const double PI;
extern const double ALPHA;

/** 'AnnularCell' IMPLEMENTATION **/
/* Static variables */
// Length variables
const double AnnularCell::m_INNER_RADIUS = INNER_RADIUS;
const double AnnularCell::m_OUTER_RADIUS = OUTER_RADIUS;

/* Methods*/

// Random number generator
std::mt19937_64 gen(10082273);
std::uniform_real_distribution<double> rndmdist(-1.0,1.0);

void  AnnularCell::printRod(const int& index) const
{
  /**
  'printRod': Print information about 'Rod' to terminal.
              Included solely for debugging purposes.
  **/

  std::cout << "Rod #" << index << ": ( " << m_bundle[index].m_xPos << " , "
                                          << m_bundle[index].m_yPos << " , "
                                          << m_bundle[index].m_angle << " ) " << std::endl;
}

Rod  AnnularCell::getRod(const int& index) const
{
  /**
    'getRod': Given the 'index' of a 'Rod', get it from the 'm_bundle'.
  **/

  return m_bundle[index];
}

bool AnnularCell::rodIsTouchingInnerWall(const int& index)
{
  /**
    'rodIsTouchingInnerWall' : Check if 'Rod' referenced by 'index'
                               is overlapping inner boundary of the 'AnnularCell'.
  **/

  // If input 'index' is valid
  if(index < m_NUMBER_OF_RODS)
  {
    // Get Rod from 'index'
    return rodIsTouchingInnerWall(getRod(index));
  }
  else // if input 'index' is not valid
  {
    return true;
  }
}

bool AnnularCell::rodIsTouchingInnerWall(const Rod& rod)
{
  /**
    'rodIsTouchingInnerWall' : Check if 'rod' is overlapping inner boundary of the 'AnnularCell'.
  **/

  // Auxiliary constant variables
  static const double R_PLUS_HALF_L = INNER_RADIUS + Rod::m_HALF_LENGTH;
  static const double R_PLUS_HALF_W = INNER_RADIUS + Rod::m_HALF_WIDTH;
  static const double HALF_D_OVER_R = Rod::m_HALF_DIAGONAL/INNER_RADIUS;
  static const double PHI_ONE = std::atan2(Rod::m_HALF_WIDTH, R_PLUS_HALF_L);
  static const double PHI_TWO = std::atan2(R_PLUS_HALF_W, Rod::m_HALF_LENGTH);
  static const double INNER_MIN_DIST = m_INNER_RADIUS + Rod::m_HALF_WIDTH;
  static const double INNER_MAX_DIST = m_INNER_RADIUS + Rod::m_HALF_DIAGONAL;

  /* First criterion: Center of rod too far or too close of inner wall to overlap */

  // Cartesian distance between 'm_aux_rod' center and 'AnnularCell' center
  double distance = std::sqrt(rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos);

  // Check first criterion
  if(distance > INNER_MAX_DIST)
  {
    return false;
  }
  if(distance < INNER_MIN_DIST)
  {
    return true;
  }

  /* Second criterion: 'distance' less than minimum distance 'minDist',
                       obtained via an analytical expression which is function
                       of relatives angles 'theta' and 'phi' (defined below)   */

  // Angle between rod center and X (horizontal) axis: 'theta'
  double theta = std::atan2(rod.m_yPos, rod.m_xPos);

  // Relative angle between 'rod's orientation and 'theta': 'phi'
  double phi   = std::abs(rod.m_angle - theta);

  // 'phi' in interval [-pi/2, pi/2]
  if(phi > PI){ phi -= PI;}
  if(phi > HALF_PI){ phi = PI - phi;}

  // Computation of analytical 'minDist'
  double minDist;
  if(phi < PHI_ONE)
  {
    // First region: 'phi' in [-pi/2, PHI_ONE]
    minDist = R_PLUS_HALF_L/std::cos(phi);
  }
  else if(phi > PHI_TWO)
  {
    // Second region: 'phi' in [-pi/2, PHI_TWO]
    minDist = R_PLUS_HALF_W/std::sin(phi);
  }
  else
  {
    // Third region: 'phi' in [PHI_ONE, PHI_TWO]
    double lambda = std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
    if(phi < ALPHA)
    {
      minDist = (Rod::m_HALF_LENGTH+m_INNER_RADIUS*std::cos(phi-lambda))/std::cos(phi);
    }
    else
    {
      minDist = (Rod::m_HALF_WIDTH+m_INNER_RADIUS*std::sin(phi-lambda))/std::sin(phi);
    }
  }

  // Check second criterion
  return (distance < minDist);
}

bool AnnularCell::rodIsTouchingOuterWall(const int& index)
{
  /**
    'rodIsTouchingOuterWall' : Check if 'Rod' referenced by 'index'
                               is overlapping outer boundary of the 'AnnularCell'.
  **/

  // Auxiliary constant variables
  static const double OUTER_MIN_DIST = std::sqrt(m_OUTER_RADIUS*m_OUTER_RADIUS-Rod::m_HALF_LENGTH*Rod::m_HALF_LENGTH) - Rod::m_HALF_WIDTH;
  static const double OUTER_MAX_DIST = m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL;

  // If input 'index' is valid
  if(index < m_NUMBER_OF_RODS)
  {
    // Get 'm_aux_rod' from 'index'
    return rodIsTouchingOuterWall(getRod(index));
  }
  else // if input 'index' is not valid
  {
    return true;
  }
}

bool AnnularCell::rodIsTouchingOuterWall(const Rod& rod)
{
  /**
    'rodIsTouchingOuterWall' : Check if 'rod' is overlapping outer boundary of the 'AnnularCell'.
  **/

  // Auxiliary constant variables
  static const double OUTER_MIN_DIST = m_OUTER_RADIUS - Rod::m_HALF_WIDTH;
  static const double OUTER_MAX_DIST = m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL;

  /* First criterion: Center of rod too far or too close of outer wall to overlap */

  // Cartesian distance between 'rod' center and 'AnnularCell' center
  double distance = std::sqrt(rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos);

  // Check first criterion
  if(distance > OUTER_MIN_DIST)
  {
    return true;
  }
  if(distance < OUTER_MAX_DIST)
  {
    return false;
  }

  /* Second criterion: 'distance' less than minimum distance
                       obtained via an analytical expression which is function
                       of relatives angles 'theta' and 'phi' (defined below)   */

  // Angle between 'rod' center and X (horizontal) axis: 'theta'
  double theta = std::atan2(rod.m_yPos, rod.m_xPos);

  // Relative angle between 'rod's orientation and 'theta': 'phi'
  double phi = rod.m_angle - theta;

  // 'phi' in interval [-pi/2, pi/2]
  if(phi < -HALF_PI)
  {
      phi = phi + PI;
  }else if(phi > HALF_PI)
  {
      phi = phi - PI;
  }

  // 'phi' is not directly needed
  phi = std::cos(ALPHA-std::fabs(phi));

  // Check second criterion
  return (distance > (std::sqrt(m_OUTER_RADIUS*m_OUTER_RADIUS-Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL*(1.0-phi*phi))-Rod::m_HALF_DIAGONAL*phi ));
}

void AnnularCell::fillAnnularCell()
{
  /**
    'fillAnnularCell': Fill the 'AnnularCell' with 'Rod's at random positions with random orientation
                       that are saved in the 'm_bundle' and referenced in the 'm_grid'.
                       If after a given number of trials is unable to introduce one 'Rod',
                       it saves its index for reference in the 'm_missingRods' array.
  **/

  // Boolean used to exit loops and mark 'Rod' as missing
  bool rodsAreTouching;

  /* Try to fill 'AnnularCell' with 'm_NUMBER_OF_PARTICLES' 'Rod's */
  for(int i=0; i < m_NUMBER_OF_RODS; i++)
  {
    // Assume new 'Rod' will be touching some of the 'rod's already included
    rodsAreTouching = true;

    // For each index, try at most 1000 positions
    for(int trials = 0; trials < 1000; trials++)
    {
      // Obtain random position and orientation for new 'm_aux_rod'
      m_aux_rod.m_xPos  = rndmdist(gen)*m_OUTER_RADIUS;
      m_aux_rod.m_yPos  = rndmdist(gen)*m_OUTER_RADIUS;
      m_aux_rod.m_angle = rndmdist(gen)*HALF_PI;

      // Check if new position is touching some of the 'AnnularCell' walls
      if(rodIsTouchingInnerWall(m_aux_rod)||rodIsTouchingOuterWall(m_aux_rod))
      {
        // ... and, if so, try a new position
        continue;
      }
      else
      {
        // Assume new 'm_aux_rod' is not touching any of the 'rod's already included
        rodsAreTouching = false;

        // For all the 'Rod's already included
        for(int j = 0; j < i; j++)
        {
          // ... check if new 'm_aux_rod' is touching it
          rodsAreTouching = m_aux_rod.isTouchingRod(getRod(j));
          // ... and if so, break this loop and try again
          if(rodsAreTouching){ break;}
        }
      }

      // If the new position is valid,
      if(!rodsAreTouching)
      {
        // ...  save 'm_aux_rod' in the 'm_bundle' of 'Rod's
        m_bundle.emplace_back(m_aux_rod);
        // ... and break the trials for this index
        break;
      }
    }

    // If index could not be inserted
    if(rodsAreTouching)
    {
      // ... save index in array of 'm_missingRods'
      m_missingRods.push_back(i);
      // ... report to console that index could not be inserted
      std::cout << "Rod #" << i << " not included" << std::endl;
      // ... and put the 'Rod' temporarily outside the 'AnnularCell' in a corner
      m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
    }
  }

  /* Fill 'm_grid' with the indexes of the 'Rod's 'm_bundle' to know roughly were each rod is */
  m_grid.fillGrid(m_bundle);
}

void AnnularCell::fillMissingRods()
{
  /**
    'fillMissingRods': Attempts to introduce the 'm_missingRods' inside the 'AnnularCell', trying several
                       random orientations for each random position obtained.
                       Note: Assumes 'fillAnnularCell' has already been called.
  **/

  // Boolean used to exit loops and mark 'Rod' as missing
  bool rodsAreTouching;

  // Array to save the indexes of the 'Rod's that could not be introduced
  std::vector<int> newMissing;

  /* Attemp to  to introduce the 'm_missingRods' inside the 'AnnularCell' */

  // For each index in 'm_missingRods'
  for(int i : m_missingRods)
  {
    // Assume missing 'Rod' will be touching some of the 'rod's already included
    rodsAreTouching = true;

    // For each index, try at most 500 positions
    for(int trials = 0; trials < 500; trials++)
    {
      // Obtain random position for missing 'Rod'
      m_aux_rod.m_xPos  = rndmdist(gen)*m_OUTER_RADIUS;
      m_aux_rod.m_yPos  = rndmdist(gen)*m_OUTER_RADIUS;

      // For each position try at most 100 orientations
      for(int n = 0; n < 100; n++)
      {
        // Obtain random orientation
        m_aux_rod.m_angle = rndmdist(gen)*HALF_PI;

        // Check if new position is touching some of the 'AnnularCell' walls
        if(rodIsTouchingInnerWall(m_aux_rod)||rodIsTouchingOuterWall(m_aux_rod))
        {
          // ... and, if so, try a new orientation
          continue;
        }
        else
        {
          // Assume new 'm_aux_rod' is not touching any of the 'rod's already included
          rodsAreTouching = false;

          // For all the 'Rod's already included
          for(int j = 0; j < NUMBER_OF_RODS; j++)
          {
            if(j!=i)
            {
              // ... check if new 'm_aux_rod' is touching it
              rodsAreTouching = m_aux_rod.isTouchingRod(getRod(j));
              // ... and if so, break this loop and try again
              if(rodsAreTouching){ break;}
            }
          }
        }

        // If the new position is valid,
        if(!rodsAreTouching)
        {
          // ... move index to new coordinates in the 'm_grid'
          m_grid.moveIndex(i,m_grid.getGridCoords(m_bundle[i]),m_grid.getGridCoords(m_aux_rod));
          // ...  save 'm_aux_rod' in the 'm_bundle' of 'Rod's
          m_bundle[i] = m_aux_rod;
          // ... and break the orientation trials for this index
          break;
        }
      }

      // If the new position was accepted
      if(!rodsAreTouching){
        // ... break the position trials for this index
        break;
      }
    }

    // If index could not be inserted
    if(rodsAreTouching)
    {
      // ... save index in array of 'newMissing' 'Rod's
      newMissing.push_back(i);
    }
  }

  // Clear the content of the old 'm_missingRods'
  m_missingRods.clear();

  // Save array of 'newMissing' 'Rod's
  m_missingRods = newMissing;

  // Report to console the number of 'Rod's still missing
  std::cout << m_missingRods.size() << " rods missing" << std::endl;
}

void AnnularCell::fillAnnularCellFromFile(std::string filepath, const int& numRodsInFile)
{
  /**
    'fillAnnularCellFromFile': Fill the 'AnnularCell' from a file with full path 'filepath'
                               with 8 columns: index, X, Y, PHI, q1, q2, q3, q4
                               were q_ are order parameters of the liquid crystal.
                               Note: only the 'X', 'Y' and 'Phi' columns are used, the rest is discarded.
                               Note: it is assumed that 'numRodsInFile' <= 'm_NUMBER_OF_RODS'.
  **/

  // Auxiliar dummy variable
  double dummy;

  // Input stream of data from 'savedConfiguration'
  std::ifstream savedConfiguration;

  // Open 'savedConfiguration'
  savedConfiguration.open(filepath);
  // For each 'Rod' of data
  for(int i=0; i < numRodsInFile; i++)
  {
    // ... save data to 'm_aux_rod'
    savedConfiguration >> dummy >> m_aux_rod.m_xPos >> m_aux_rod.m_yPos >> m_aux_rod.m_angle
                       >> dummy >> dummy >> dummy >> dummy;
    // ... and put 'm_aux_rod' into 'm_bundle' of 'Rod's
    m_bundle.emplace_back(m_aux_rod);
  }

  // Close 'savedConfiguration'
  savedConfiguration.close();
  // If 'savedConfiguration' did not contain all the desired 'Rod's
  for(int i=numRodsInFile; i<NUMBER_OF_RODS; i++)
  {
    // ... save index as 'm_missingRods'
    m_missingRods.push_back(i);
    // ... and put the 'Rod' temporarily outside the 'AnnularCell' in a corner
    m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
  }

  // Fill 'm_grid' with the indexes of the 'Rod's 'm_bundle' to know roughly were each rod is
  m_grid.fillGrid(m_bundle);
}

void AnnularCell::saveAnnularCellToFile(std::string filename)
{
  /**
    'saveAnnularCellToFile': Save the 'AnnularCell' in CSV format in file 'filename'.
                             with 8 columns: index, X, Y, PHI, q1, q2, q3, q4
                             were q_ are order parameters of the liquid crystal.
  **/

  // Ouput stream
  std::ofstream configuration;
  configuration.open(filename);

  // Headers row
  configuration << "index" << "," << "x" << "," << "y" << "," << "angle" << '\r' << '\n';;

  // Save 'Rod's information
  for(int i=0; i < m_NUMBER_OF_PARTICLES; i++)
  {
    m_aux_rod = m_bundle[i];
      configuration << i << "," << m_aux_rod.m_xPos << "," << m_aux_rod.m_yPos << "," << m_aux_rod.m_angle << '\r' << '\n';;
  }

  // Close output stream
  configuration.close();
}
