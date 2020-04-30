/**
  * "AnnularCell.h":
  * --------------------
  * Declaration of the AnnularCell class.
  * For implementation details, see "AnnularCell.cpp".
  *
  * Needs: "Rod.h", "Grid.h".
  *
  * --------------------
  * An "AnnularCell" is an the region between
  * two concentric circles of radii
  *   "m_INNER_RADIUS" and
  *   "m_OUTER_RADIUS".
  * It contains a "m_bundle" of "m_NUMBER_OF_PARTICLES" rods,
  * complemented with a "m_grid" to point to their locations.
  * An "AnnularCell" can be filled with a pseudorandom algorithm,
  * or from a specific type of file (see below).
  * The rods on the "m_bundle" can be saved to a CSV file.
  *
  * --------------------
  * Last modified: 2020-04-29
  * By: M. E. Maza-Cuello
  */

#ifndef ANNULARCELL_H_INCLUDED
#define ANNULARCELL_H_INCLUDED

#include <vector>
#include <string>

#include "Rod.h"
#include "Grid.h"

extern const int NUMBER_OF_RODS;

class AnnularCell
{
public:
  static const double m_INNER_RADIUS;
  static const double m_OUTER_RADIUS;
  const  int   m_NUMBER_OF_PARTICLES;
  Rod          m_aux_rod;
  std::vector<Rod> m_bundle;
  Grid         m_grid;
  std::vector<int> m_missingRods;

  AnnularCell()
    : m_NUMBER_OF_PARTICLES(NUMBER_OF_RODS)
  {
    m_aux_rod = Rod();
    m_bundle.reserve(NUMBER_OF_RODS);
  }

  Rod  getRod(const int& index)   const;

  bool rodIsOutsideWalls(const int& index);
  bool rodIsOutsideWalls(const Rod& rod);

  bool rodIsTouchingInnerWall(const int& index);
  bool rodIsTouchingInnerWall(const Rod& rod);
  bool rodIsTouchingOuterWall(const int& index);
  bool rodIsTouchingOuterWall(const Rod& rod);

  /**
    * Fill annular cell with Rods at random positions.
    * After a maximum number of attempts, those rods not included
    * are saved in "m_missingRods" and a message is displayed.
    * Internally, these Rods do exist, they are located outside the cell.
    */
  void fill();

  /**
    * Attempts to move missing rods inside the cell.
    * Upon exit, displays how many rods are still in "m_missingRods".
    */
  void fillMissingRods();

  /**
    * Fill with rods from file with eight columns:
    * index, X, Y, PHI, q1, q2, q3, q4
    * were q_ are order parameters of the bundle of rods.
    * Only X, Y, and PHI are used to create a Rod(X,Y,PHI).
    * Note: erases all previous rods on the cell.
    */
  void fillFromFile(std::string filename, const int& numRodsInFile);

  /**
    * Saves bundle of rods to CSV with four columns:
    * index,x,y,angle.
    * First row: header names
    */
  void saveToCSV(std::string filename);
};

#endif // ANNULARCELL_H_INCLUDED
