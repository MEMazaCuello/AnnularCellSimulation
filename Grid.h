/****
  'Grid.h':
  --------------------
  Declaration of 'Grid' class.
  For implementation details, see 'Grid.cpp'.
  
  Needs: 'Rods.h'
  
  --------------------
  A 'Grid' is a 'm_map' of 'Rod's indexes, 
  divided in square 'm_BOXES' with 
    side length = 'm_BOX_LENGTH'.
        
  --------------------
  Last modified: 2019-05-14 
  By: M. E. Maza Cuello
****/

#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Rod.h"
#include <vector>

/** 'Grid' DECLARATION **/
class Grid
{
public:
  /* Static variables */
  // Number of boxes
  static const int    m_BOXES;
  // Length of boxes
  static const double m_BOX_LENGTH;
  static const double m_INV_BOX_LENGTH;
  static const double m_HALF_LENGTH;
  
  /* Public variables */
  // Basic map 
  std::vector< std::vector< std::vector<int> > > m_map;
  // Coordinates of neightbouring boxes of each map box
  std::vector< std::vector< std::vector<int> > > m_NB_XCOORD;
  std::vector< std::vector< std::vector<int> > > m_NB_YCOORD;
  // Vector of indexes saved at a box and in first neighbouring boxes
  std::vector<int>    m_neighbours;

  /* Constructors */
  // Default
  Grid();

  /* Methods */
  
  int getGridX(const Rod& rod) const;
  /**
    'getGridX': Given a 'rod', compute horizontal coordinate in grid,
                to be used as first index for 'm_map'. 
  **/
  
  int getGridY(const Rod& rod) const;
  /**
    'getGridY': Given a 'rod', compute vertical coordinate in grid,
                to be used as second index for 'm_map'. 
  **/
  
  std::vector<int> getGridCoords(const Rod& rod) const;
  /**
    'getGridCoords': Given a 'rod', compute pair of coordinates in grid,
                     to be used as indexes for 'm_map'. 
  **/
  
  void fillGrid(std::vector<Rod>& bundle);
  /** 
    'fillGrid': Given a bundle of 'Rod's, compute coordinates in the grid of each 'Rod'
                and save 'Rod's index on corresponding box on 'm_map'.
  **/
  
  int  getNumberOfNeighbours(const int& coordx, const int& coordy) const;
  /**
    'getNumberOfNeighbours': Given the coordinates ('coordx', 'coordy') of a box in the grid,
                             obtain number of indexes saved on 'm_map'['coordx']['coordy'] and in
                             the (generally 8) first-neighbouring boxes around it.  
  **/  
  
  std::vector<int> getNeighbours(const int& coordx, const int& coordy);
  /**
    'getNeighbours': Given the coordinates ('coordx', 'coordy') of a box in the grid,
                     obtain indexes saved on 'm_map'['coordx']['coordy'] and in
                     the (generally 8) first-neighbouring boxes around it.  
  **/
  
  std::vector<int> getNeighbours(const std::vector<int>& coords);
  /**
    'getNeighbours': Given the pair of coordinates 'coords' of a box in the grid,
                     obtain indexes saved on 'm_map'['coords[0]']['coords[1]'] and in
                     the (generally 8) first-neighbouring boxes around it.  
  **/
  
  void moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy);
  /**
    'moveIndex': Update the 'm_map' by moving the index 'movingIdx'
                 from the box at ('oldx','oldy')
                 to the box at ('newx','newy').  
  **/
  
  void moveIndex(const int& movingIdx, const std::vector<int>& oldCoords, const std::vector<int>& newCoords);
  /**
    'moveIndex': Update the 'm_map' by moving the index 'movingIdx'
                 from the box at ('oldCoords[0]','oldCoords[1]')
                 to the box at ('newCoords[0]','newCoords[1]').
  **/
};

#endif // GRID_H_INCLUDED
