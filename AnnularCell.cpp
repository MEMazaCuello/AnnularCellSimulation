/**
  * "AnnularCell.cpp":
  * --------------------
  * Implementation of the AnnularCell class.
  * For declaration details, see "AnnularCell.cpp".
  *
  * Needs: "parameters.cpp", "AnnularCell.h".
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

#include "AnnularCell.h"

// Parameters defined in "parameters.cpp" _____________________________
extern const double HALF_PI;
extern const double PI;
extern const double ALPHA;
extern const double INNER_RADIUS;
extern const double OUTER_RADIUS;

// Static variables ___________________________________________________
const double AnnularCell::m_INNER_RADIUS = INNER_RADIUS;
const double AnnularCell::m_OUTER_RADIUS = OUTER_RADIUS;

// Methods ____________________________________________________________

// Random number generator
std::mt19937_64 gen(99994825);
std::uniform_real_distribution<double> rndmdist(-1.0,1.0);

Rod  AnnularCell::getRod(const int& index) const
{
  return m_bundle[index];
}

bool AnnularCell::rodIsOutsideWalls(const int& index) const
{
  return rodIsOutsideWalls(getRod(index));
}

bool AnnularCell::rodIsOutsideWalls(const Rod& rod) const
{
  return (rodIsTouchingInnerWall(rod)|| rodIsTouchingOuterWall(rod));
}

bool AnnularCell::internalRodIsOutsideWalls() const
{
  return (internalRodIsTouchingInnerWall()|| internalRodIsTouchingOuterWall());
}

bool AnnularCell::rodIsTouchingInnerWall(const int& index) const
{
  return rodIsTouchingInnerWall(getRod(index));
}

bool AnnularCell::rodIsTouchingInnerWall(const Rod& rod) const
{
  // Auxiliary parameters needed for analytical constraints
  static const double R_PLUS_HALF_L = INNER_RADIUS + Rod::m_HALF_LENGTH;
  static const double R_PLUS_HALF_W = INNER_RADIUS + Rod::m_HALF_WIDTH;
  static const double HALF_D_OVER_R = Rod::m_HALF_DIAGONAL/INNER_RADIUS;
  static const double PHI_ONE = std::atan2(Rod::m_HALF_WIDTH, R_PLUS_HALF_L);
  static const double PHI_TWO = std::atan2(R_PLUS_HALF_W, Rod::m_HALF_LENGTH);
  static const double INNER_MIN_SQ_DIST = (m_INNER_RADIUS + Rod::m_HALF_WIDTH)*(m_INNER_RADIUS + Rod::m_HALF_WIDTH);
  static const double INNER_MAX_SQ_DIST = (m_INNER_RADIUS + Rod::m_HALF_DIAGONAL)*(m_INNER_RADIUS + Rod::m_HALF_DIAGONAL);

  double sqDist = rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos;

  if (sqDist < INNER_MIN_SQ_DIST)
  {
    return true;
  }
  else if (sqDist > INNER_MAX_SQ_DIST)
  {
    return false;
  }

  // Relative angles
  double theta = std::atan2(rod.m_yPos, rod.m_xPos);
  double phi   = std::abs(rod.m_angle - theta);

  // phi between 0 and PI/2
  if (phi > PI)
  {
    phi -= PI;
  }
  if (phi > HALF_PI)
  {
    phi = PI - phi;
  }

  // Analytically computed minimum distance
  if (phi < PHI_ONE)
  {
    double minDist = R_PLUS_HALF_L/std::cos(phi);
    return (sqDist < minDist*minDist);
  }
  else if (phi > PHI_TWO)
  {
    double minDist = R_PLUS_HALF_W/std::sin(phi);
    return (sqDist < minDist*minDist);
  }
  else
  {
    double minDist = phi-std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
    if (phi < ALPHA)
    {
      minDist = (Rod::m_HALF_LENGTH+m_INNER_RADIUS*std::cos(minDist))/std::cos(phi);
    }
    else
    {
      minDist = (Rod::m_HALF_WIDTH+m_INNER_RADIUS*std::sin(minDist))/std::sin(phi);
    }
    return (sqDist < minDist*minDist);
  }
}

bool AnnularCell::internalRodIsTouchingInnerWall() const
{
  // Auxiliary parameters needed for analytical constraints
  static const double R_PLUS_HALF_L = INNER_RADIUS + Rod::m_HALF_LENGTH;
  static const double R_PLUS_HALF_W = INNER_RADIUS + Rod::m_HALF_WIDTH;
  static const double HALF_D_OVER_R = Rod::m_HALF_DIAGONAL/INNER_RADIUS;
  static const double PHI_ONE = std::atan2(Rod::m_HALF_WIDTH, R_PLUS_HALF_L);
  static const double PHI_TWO = std::atan2(R_PLUS_HALF_W, Rod::m_HALF_LENGTH);
  static const double INNER_MIN_SQ_DIST = (m_INNER_RADIUS + Rod::m_HALF_WIDTH)*(m_INNER_RADIUS + Rod::m_HALF_WIDTH);
  static const double INNER_MAX_SQ_DIST = (m_INNER_RADIUS + Rod::m_HALF_DIAGONAL)*(m_INNER_RADIUS + Rod::m_HALF_DIAGONAL);

  double sqDist = m_aux_rod.m_xPos*m_aux_rod.m_xPos + m_aux_rod.m_yPos*m_aux_rod.m_yPos;

  if (sqDist < INNER_MIN_SQ_DIST)
  {
    return true;
  }
  else if (sqDist > INNER_MAX_SQ_DIST)
  {
    return false;
  }

  // Relative angles
  double theta = std::atan2(m_aux_rod.m_yPos, m_aux_rod.m_xPos);
  double phi   = std::abs(m_aux_rod.m_angle - theta);

  // phi between 0 and PI/2
  if (phi > PI)
  {
    phi -= PI;
  }
  if (phi > HALF_PI)
  {
    phi = PI - phi;
  }

  // Analytically computed minimum distance
  if (phi < PHI_ONE)
  {
    double minDist = R_PLUS_HALF_L/std::cos(phi);
    return (sqDist < minDist*minDist);
  }
  else if (phi > PHI_TWO)
  {
    double minDist = R_PLUS_HALF_W/std::sin(phi);
    return (sqDist < minDist*minDist);
  }
  else
  {
    double minDist = phi-std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
    if (phi < ALPHA)
    {
      minDist = (Rod::m_HALF_LENGTH+m_INNER_RADIUS*std::cos(minDist))/std::cos(phi);
    }
    else
    {
      minDist = (Rod::m_HALF_WIDTH+m_INNER_RADIUS*std::sin(minDist))/std::sin(phi);
    }
    return (sqDist < minDist*minDist);
  }
}

bool AnnularCell::rodIsTouchingOuterWall(const int& index) const
{
  return rodIsTouchingOuterWall(getRod(index));
}

bool AnnularCell::rodIsTouchingOuterWall(const Rod& rod) const
{
  static const double SQ_OUTER_RADIUS = OUTER_RADIUS*OUTER_RADIUS;
  static const double SQ_HALF_DIAGONAL = Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL;
  static const double OUTER_MIN_SQ_DIST = (m_OUTER_RADIUS - Rod::m_HALF_WIDTH)*(m_OUTER_RADIUS - Rod::m_HALF_WIDTH);
  static const double OUTER_MAX_SQ_DIST = (m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL)*(m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL);

  double sqDist = rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos;

  if (sqDist > OUTER_MIN_SQ_DIST)
  {
    return true;
  }
  else if (sqDist < OUTER_MAX_SQ_DIST)
  {
    return false;
  }

  // Analytically computed minimum distance
  double phi = rod.m_angle - std::atan2(rod.m_yPos, rod.m_xPos);

  if (phi < -HALF_PI)
  {
    phi = phi + PI;
  }
  else if (phi > HALF_PI)
  {
    phi = phi - PI;
  }

  phi = std::cos( ALPHA - std::fabs(phi) );
  phi = std::sqrt( SQ_OUTER_RADIUS-SQ_HALF_DIAGONAL*(1.0d-phi*phi) ) - Rod::m_HALF_DIAGONAL*phi;

  return (sqDist > phi*phi);
}

bool AnnularCell::internalRodIsTouchingOuterWall() const
{
  static const double SQ_OUTER_RADIUS = OUTER_RADIUS*OUTER_RADIUS;
  static const double SQ_HALF_DIAGONAL = Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL;
  static const double OUTER_MIN_SQ_DIST = (m_OUTER_RADIUS - Rod::m_HALF_WIDTH)*(m_OUTER_RADIUS - Rod::m_HALF_WIDTH);
  static const double OUTER_MAX_SQ_DIST = (m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL)*(m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL);

  double sqDist = m_aux_rod.m_xPos*m_aux_rod.m_xPos + m_aux_rod.m_yPos*m_aux_rod.m_yPos;

  if (sqDist > OUTER_MIN_SQ_DIST)
  {
    return true;
  }
  else if (sqDist < OUTER_MAX_SQ_DIST)
  {
    return false;
  }

  // Analytically computed minimum distance
  double phi = m_aux_rod.m_angle - std::atan2(m_aux_rod.m_yPos, m_aux_rod.m_xPos);

  if (phi < -HALF_PI)
  {
    phi = phi + PI;
  }
  else if (phi > HALF_PI)
  {
    phi = phi - PI;
  }

  phi = std::cos(ALPHA-std::fabs(phi));
  phi = std::sqrt(SQ_OUTER_RADIUS-SQ_HALF_DIAGONAL*(1.0d-phi*phi))-Rod::m_HALF_DIAGONAL*phi;

  return (sqDist > phi*phi);
}

bool AnnularCell::isInternalPositionValid(const int& index)
{
  if ( internalRodIsOutsideWalls() )
  {
    return false;
  }

  // New grid coordinates
  int newi = m_grid.getGridIdx(m_aux_rod.m_xPos);
  int newj = m_grid.getGridIdx(m_aux_rod.m_yPos);

  // Fills m_grid.m_neighbors to contain new neighbors of index
  m_grid.setNeighbors(newi,newj);

  for (int& idx : m_grid.m_neighbors)
  {
    if (idx != index)  // Rod is neighbor of itself
    {
      if (m_aux_rod.isTouchingRod(m_bundle[idx]))
      {
        return false;
      }
    }
  }

  // Update grid coordinates of rod
  m_grid.moveIndex(index,
                   m_grid.getGridIdx(m_bundle[index].m_xPos), m_grid.getGridIdx(m_bundle[index].m_yPos),
                   newi,newj);
  m_bundle[index] = m_aux_rod;

  return true;
}

void AnnularCell::fill()
{
  bool rodsAreTouching;
  for (int i = 0; i < m_NUMBER_OF_PARTICLES; i++)
  {
    rodsAreTouching = true;
    for (int trials = 0; trials < 1000; trials++) // 1000 trials per index
    {
      // Random position over square of size OUTER_RADIUS
      m_aux_rod.m_xPos  = rndmdist(gen)*m_OUTER_RADIUS;
      m_aux_rod.m_yPos  = rndmdist(gen)*m_OUTER_RADIUS;
      m_aux_rod.m_angle = rndmdist(gen)*HALF_PI;

      // Check it is inside the cell
      if ( rodIsTouchingInnerWall(m_aux_rod) || rodIsTouchingOuterWall(m_aux_rod) )
      {
        continue;
      }

      // Check overlapping
      rodsAreTouching = false;
      for (int j = 0; j < i; j++)
      {
        rodsAreTouching = m_aux_rod.isTouchingRod(getRod(j));
        if (rodsAreTouching)
        {
          break;
        }
      }

      // Save rod
      if (!rodsAreTouching)
      {
        m_bundle.emplace_back(m_aux_rod);
        break;
      }
    }

    // Create vector or the indexes that could not be inserted
    if (rodsAreTouching)
    {
      m_missingRods.push_back(i);
      std::cout << "Rod #" << i << " not included" << std::endl;
      m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
    }
  }

  // Fill grid
  m_grid.fill(m_bundle);

  std::cout << std::endl << "\t CELL FILLED" << std::endl;
}

void AnnularCell::fillMissingRods()
{
  bool rodsAreTouching;
  std::vector<int> newMissing;

  double c;
  double s;

  // Number of particles that were entered in the cell
  int insideParticles = m_NUMBER_OF_PARTICLES - m_missingRods.size();

  for (int i : m_missingRods)
  {
    rodsAreTouching = true;
    for (int p = 0; p < insideParticles; p++)
    {
      c = std::cos(m_bundle[p].m_angle);
      s = std::sin(m_bundle[p].m_angle);

      for (int q = 0; q < 4; q++)
      {
        m_aux_rod = m_bundle[p];

        // Try copy of rod displaced in four Cartesian directions
        switch (q)
        {
          case 0:
            m_aux_rod.m_xPos += 2.0*Rod::m_HALF_LENGTH*c + 0.0001d;
            m_aux_rod.m_yPos += 2.0*Rod::m_HALF_LENGTH*s + 0.0001d;
          case 1:
            m_aux_rod.m_xPos -= 2.0*Rod::m_HALF_LENGTH*c + 0.0001d;
            m_aux_rod.m_yPos -= 2.0*Rod::m_HALF_LENGTH*s + 0.0001d;
          case 2:
            m_aux_rod.m_xPos -= 2.0*Rod::m_HALF_WIDTH*c + 0.0001d;
            m_aux_rod.m_yPos += 2.0*Rod::m_HALF_WIDTH*s + 0.0001d;
          case 3:
            m_aux_rod.m_xPos -= 2.0*Rod::m_HALF_WIDTH*c + 0.0001d;
            m_aux_rod.m_yPos += 2.0*Rod::m_HALF_WIDTH*s + 0.0001d;
        }

        // Check if position is valid
        if (rodIsTouchingInnerWall(m_aux_rod)||rodIsTouchingOuterWall(m_aux_rod))
        {
          continue;
        }
        else
        {
          m_grid.setNeighbors(m_grid.getCoords(m_aux_rod));

          rodsAreTouching = false;
          for (int nb : m_grid.m_neighbors)
          {
            rodsAreTouching = m_aux_rod.isTouchingRod(getRod(nb));
            if (rodsAreTouching)
            {
              break;
            }
          }
        }

        // Save rod
        if (!rodsAreTouching)
        {
          // Update grid: rod has moved to another box of the grid
          m_grid.moveIndex(insideParticles,m_grid.getCoords(m_bundle[insideParticles]),m_grid.getCoords(m_aux_rod));
          // One more particle inside that needs to be checked for overlapping
          m_bundle[insideParticles] = m_aux_rod;
          insideParticles++;
          break;
        }
      }

      if (!rodsAreTouching)
      {
        break;
      }
    }

    // Those particles missing are saved again (for eventual recursion)
    if (rodsAreTouching)
    {
      newMissing.push_back(i);
    }
  }

  m_missingRods.clear();
  m_missingRods = newMissing;
  std::cout << m_missingRods.size() << " rods missing" << std::endl;
}

void AnnularCell::fillFromFile(std::string filepath, const int& numRodsInFile)
{
  // Empty cell
  m_bundle.clear();
  m_missingRods.clear();

  // Fill cell
  std::ifstream savedconfiguration;
  savedconfiguration.open(filepath);

  // Assumes numRodsInFile <= NUMBER_OF_RODS
  m_bundle.reserve(NUMBER_OF_RODS);

  double dummy;
  for (int i = 0; i < numRodsInFile; i++)
  {
    savedconfiguration >> dummy >> m_aux_rod.m_xPos >> m_aux_rod.m_yPos >> m_aux_rod.m_angle
                       >> dummy >> dummy >> dummy >> dummy;

    m_bundle.emplace_back(m_aux_rod);
  }
  savedconfiguration.close();

  // Missing particles (if numRodsInFile < NUMBER_OF_RODS)
  for (int i = numRodsInFile; i < NUMBER_OF_RODS; i++)
  {
    m_missingRods.push_back(i);
    m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
  }

  // Fill grid
  m_grid.fill(m_bundle);
}

void AnnularCell::saveToCSV(std::string filename)
{
  std::ofstream configurationToSave;
  configurationToSave.open(filename);
  configurationToSave << "index" << "," << "x" << "," << "y" << "," << "angle" << std::endl;
  for (int i = 0; i < m_NUMBER_OF_PARTICLES; i++)
  {
    m_aux_rod = m_bundle[i];
    configurationToSave << i << "," << m_aux_rod.m_xPos << "," << m_aux_rod.m_yPos << ","
                        << m_aux_rod.m_angle << std::endl;
  }
  configurationToSave.close();
}
