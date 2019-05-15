/****
  'AnnularCell.h':
  --------------------
  Declaration of 'AnnularCell' class.
  For implementation details, see 'AnnularCell.cpp'.
  
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
  
  --------------------
  Last modified: 2019-05-15 
  By: M. E. Maza Cuello
****/

#ifndef ANNULARCELL_H_INCLUDED
#define ANNULARCELL_H_INCLUDED

#include <vector>
#include <iostream>
#include <string>
#include "Rod.h"
#include "Grid.h"

/** PARAMETERS OBTAINED FROM 'parameters.cpp' **/
// Number of 'Rod's
extern const int NUMBER_OF_RODS;

/** 'AnnularCell' DECLARATION **/
class AnnularCell
{
public:
  /* Static variables */
  // Number of 'Rod's   
  const int m_NUMBER_OF_RODS;
  // Lengths
  static const double m_INNER_RADIUS;
  static const double m_OUTER_RADIUS;
  
  /* Public variables */
  // 'Rod's contained in the 'Annular Cell'
  std::vector<Rod> m_bundle;
  // Auxiliar 'Rod'
  Rod m_aux_rod;
  // 'Grid' of indexes to know approximate position of any 'Rod' 
  Grid m_grid;
  // Indexes references to 'Rod's that could not be inserted 
  // with the filling algorithm of the 'AnnularCell'
  std::vector<int> m_missingRods;

  /* Constructors */
  // Default
  AnnularCell()
      : m_NUMBER_OF_RODS(NUMBER_OF_RODS)
  {
      m_aux_rod = Rod();
      m_bundle.reserve(NUMBER_OF_RODS);
  }

  /* Methods */

  void printRod(const int& index) const;
  /**
    'printRod': Print information about 'Rod' to terminal. 
                Included solely for debugging purposes.
  **/  
  
  Rod  getRod(const int& index)   const;
  /**
    'getRod': Given the 'index' of a 'Rod', get it from the 'm_bundle'. 
  **/

  bool rodIsTouchingInnerWall(const int& index);
  /**
    'rodIsTouchingInnerWall' : Check if 'Rod' referenced by 'index' 
                               is overlapping inner boundary of the 'AnnularCell'.
  **/
  
  bool rodIsTouchingInnerWall(const Rod& rod);
  /**
    'rodIsTouchingInnerWall' : Check if 'rod' is overlapping inner boundary of the 'AnnularCell'.
  **/
  
  bool rodIsTouchingOuterWall(const int& index);
  /**
    'rodIsTouchingOuterWall' : Check if 'Rod' referenced by 'index' 
                               is overlapping outer boundary of the 'AnnularCell'.
  **/
  
  bool rodIsTouchingOuterWall(const Rod& rod);  
  /**
    'rodIsTouchingOuterWall' : Check if 'rod' is overlapping outer boundary of the 'AnnularCell'.
  **/
  
  void fillAnnularCell();
  /**
    'fillAnnularCell': Fill the 'AnnularCell' with 'Rod's at random positions with random orientation
                       that are saved in the 'm_bundle' and referenced in the 'm_grid'.
                       If after a given number of trials is unable to introduce one 'Rod',
                       it saves its index for reference in the 'm_missingRods' array.
  **/
  
  void fillMissingRods();
  /**
    'fillMissingRods': Attempts to introduce the 'm_missingRods' inside the 'AnnularCell', trying several
                       random orientations for each random position obtained.
  **/
  
  void fillAnnularCellFromFile(std::string filepath, const int& numRodsInFile);
  /**
    'fillAnnularCellFromFile': Fill the 'AnnularCell' from a file with full path 'filepath' 
                               with 8 columns: index, X, Y, PHI, q1, q2, q3, q4
                               were q_ are order parameters of the liquid crystal.
                               Note: only the 'X', 'Y' and 'Phi' are used, the res tis discarded.
                               Note: it is assumed that 'numRodsInFile' <= 'm_NUMBER_OF_RODS'
  **/
};

#endif // ANNULARCELL_H_INCLUDED
