/****
  'Grid.cpp':
  --------------------
  Implementation of 'Grid' class.
  For declaration details, see 'Grid.cpp'.

  Needs: 'Rods.h'

  --------------------
  A 'Grid' is a 'm_map' of 'Rod's indexes, 
  divided in square 'm_BOXES' with 
    side length = 'm_BOX_LENGTH'.
  It is used to optimise the checking of 
  overlapping of a 'Rod' with the rest of
  the 'Rod's, taking into account only those
  standing in neighbouring boxes.

  Note: the file 'parameters.cpp' is needed to specify the parameters
        'DIAGONAL' and 'OUTER_RADIUS'.

  Note: each 'm_BOXES' is neighbour of itself. 
  --------------------
  Last modified: 2019-05-14 
  By: M. E. Maza Cuello
****/

#include "Grid.h"
#include <cmath>
#include <iostream>

/** PARAMETERS OBTAINED FROM 'parameters.cpp' **/
// Lengths
extern const double DIAGONAL;
extern const double OUTER_RADIUS;

/** 'Grid' IMPLEMENTATION **/
/* Static variables */
// Number of 'm_BOXES' per side
const int    Grid::m_BOXES = int(std::ceil(2.0d*OUTER_RADIUS/Grid::m_BOX_LENGTH)) + 1;
// Length variables
const double Grid::m_BOX_LENGTH     = DIAGONAL;
const double Grid::m_HALF_LENGTH    = 0.5d*Grid::m_BOXES*Grid::m_BOX_LENGTH;
const double Grid::m_INV_BOX_LENGTH = 1.0d/Grid::m_BOX_LENGTH;

/* Constructors */
// Default
Grid::Grid()
{
  /* Create 'm_map' as an empty 2D matrix of dimension 'm_BOXES' * 'm_BOXES' */
  m_map.resize(m_BOXES);
  for(auto& v : m_map){ v.resize(m_BOXES);}

  /* Create 'm_NB_XCOORD' and 'm_NB_YCOORD' to save the coordinates
     of the adjacent boxes (horizontal, vertical and diagonal directions) 
     for each 'm_BOXES', to optimise the extraction of neighbouring 'm_BOXES'.
     'm_NB_XCOORD' will contain indexes in horizontal direction, and
     'm_NB_YCOORD' will contain indexes in vertical direction. 
     
     Note: each 'm_BOXES' is neighbour of itself. */
  
  // Declare 'm_NB_XCOORD' and 'm_NB_YCOORD' with the same size as 'm_map'
  m_NB_XCOORD = m_map;
  m_NB_YCOORD = m_map;

  // Fill each grid point with coordinates of adjacent 'm_boxes'.
  int grid_boxes_minus_one = m_BOXES - 1;
  for(int i = 0; i < m_BOXES; i++) // For each column in 'm_map'
  {
    if(i==0) // Leftmost column of 'm_map'
    {
        // Top left corner of 'm_map'
        m_NB_XCOORD[i][0] = {i,i+1};
        m_NB_YCOORD[i][0] = {0,1};

        // General row
        for(int j = 1; j < grid_boxes_minus_one; j++)
        {
            m_NB_XCOORD[i][j] = {i,i+1};
            m_NB_YCOORD[i][j] = {j-1,j,j+1};
        }

        // Bottom left corner of 'm_map'
        m_NB_XCOORD[i][grid_boxes_minus_one] = {i,i+1};
        m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
    else if(i==grid_boxes_minus_one) // Rightmost column of 'm_map'
    {
        // Top right corner of 'm_map'
        m_NB_XCOORD[i][0] = {i-1,i};
        m_NB_YCOORD[i][0] = {0,1};

        // General row of 'm_map'
        for(int j = 1; j < grid_boxes_minus_one; j++)
        {
            m_NB_XCOORD[i][j] = {i-1,i};
            m_NB_YCOORD[i][j] = {j-1,j,j+1};
        }

        // Bottom right corner of 'm_map'
        m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i};
        m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
    else  // General column of 'm_map'
    {
        // Top row of 'm_map'
        m_NB_XCOORD[i][0] = {i-1,i,i+1};
        m_NB_YCOORD[i][0] = {0,1};

        // General row of 'm_map'
        for(int j = 1; j < grid_boxes_minus_one; j++)
        {
            m_NB_XCOORD[i][j] = {i-1,i,i+1};
            m_NB_YCOORD[i][j] = {j-1,j,j+1};
        }

        // Bottom row of 'm_map'
        m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i,i+1};
        m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
  }
}

/* Methods */

inline int Grid::getGridX(const Rod& rod) const
{
  /**
    'getGridX': Given a 'rod', compute X (horizontal) coordinate in grid,
                to be used as first index for 'm_map'. 
  **/
  
  // Return X coordinate
  return int(std::floor((rod.m_xPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

inline int Grid::getGridY(const Rod& rod) const
{
  /**
    'getGridY': Given a 'rod', compute Y (vertical) coordinate in grid,
                to be used as second index for 'm_map'. 
  **/
  
  // Return X coordinate
  return int(std::floor((rod.m_yPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

std::vector<int> Grid::getGridCoords(const Rod& rod) const
{
  /**
    'getGridCoords': Given a 'rod', compute pair (X,Y) of coordinates in grid,
                     to be used as indexes for 'm_map'. 
  **/
  
  // Return (X,Y) coordinates
  return std::vector<int> {int(std::floor((rod.m_xPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH)),
                           int(std::floor((rod.m_yPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH))};
}

void Grid::fillGrid(std::vector<Rod>& bundle)
{
  /** 
    'fillGrid': Given a 'bundle' of 'Rod's, compute coordinates in the grid of each 'Rod'
                and save 'Rod's index on corresponding box on 'm_map'.
  **/
  
  // 'Rod's index
  int index = 0;
  for(const Rod& rod : bundle)
  {
      // Push the index of the 'rod' in the correspnding 'm_BOXES' of the 'm_map' 
      m_map[getGridX(rod)][getGridY(rod)].push_back(index);
      // Next index
      index++;
  }
}

int Grid::getNumberOfNeighbours(const int& coordx, const int& coordy) const
{
  /**
    'getNumberOfNeighbours': Given the coordinates ('coordx', 'coordy') of a box in the grid,
                             obtain number of indexes saved on 'm_map'['coordx']['coordy'] and in
                             the (generally 8) first-neighbouring boxes around it.  
  **/
  
  // Number of neighbours
  int number = 0;

  // Get neighbour's X coordinate 
  for(const int& idx : m_NB_XCOORD[coordx][coordy])
  {
    // Get neighbour's Y coordinate
    for(const int& idy : m_NB_YCOORD[coordx][coordy])
    {
      // Sum number of indexes in neighbouring 'm_BOXES'
      number += m_map[idx][idy].size();
    }
  }

  // Return number of neighbours
  return number;
}

std::vector<int> Grid::getNeighbours(const int& coordx, const int& coordy)
{
  /**
    'getNeighbours': Given the coordinates ('coordx', 'coordy') of a box in the grid,
                     obtain indexes saved on 'm_map'['coordx']['coordy'] and in
                     the (generally 8) first-neighbouring boxes around it.  
  **/
  
  // Clear previous vector of indexes
  m_neighbours.clear();
  
  // Reserve space for indexes in neighbours
  m_neighbours.reserve(getNumberOfNeighbours(coordx,coordy));

  // Get neighbour's X coordinate 
  for(const int& idx : m_NB_XCOORD[coordx][coordy])
  {
    // Get neighbour's Y coordinate 
    for(const int& idy : m_NB_YCOORD[coordx][coordy])
    {
      // Get indexes in neighbour 'm_BOXES'
      for(const int& idz : m_map[idx][idy])
      {
        // Save indexes
        m_neighbours.emplace_back(idz);
      }
    }
  }

  // Return vector of indexes
  return m_neighbours;
}

std::vector<int> Grid::getNeighbours(const std::vector<int>& coords)
{
  /**
    'getNeighbours': Given the pair of coordinates 'coords' of a box in the grid,
                     obtain indexes saved on 'm_map'['coords[0]']['coords[1]'] and in
                     the (generally 8) first-neighbouring boxes around it.  
  **/

  // Clear previous vector of indexes
  m_neighbours.clear();
  
  // Reserve space for indexes in neighbours
  m_neighbours.reserve(getNumberOfNeighbours(coords[0],coords[1]));

  // Get neighbour's X coordinate 
  for(const int& idx : m_NB_XCOORD[coords[0]][coords[1]])
  {
    // Get neighbour's Y coordinate
    for(const int& idy : m_NB_YCOORD[coords[0]][coords[1]])
    {
      // Get indexes in neighbour 'm_BOXES'
      for(const int& idz : m_map[idx][idy])
      {
        // Save indexes
        m_neighbours.emplace_back(idz);
      }
    }
  }

  // Return vector of indexes
  return m_neighbours;
}

void Grid::moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy)
{
  /**
    'moveIndex': Update the 'm_map' by moving the index 'movingIdx'
                 from the box at ('oldx','oldy')
                 to the box at ('newx','newy').  
  **/
  
  // If the 'm_map' needs to be updated
  if(oldx != newx || oldy != newy)
  {
    // Error boolean
    bool indexError = true;
    
    // Get indexes in neighbour 'm_BOXES'
    for(auto& index : m_map[oldx][oldy])
    {
      // When 'movingIdx' is found
      if(index == movingIdx)
      {
        // Move 'movingIdx' at the end of the array...
        index = m_map[oldx][oldy].back();
        // ... and erase it
        m_map[oldx][oldy].pop_back();
        
        // Save 'movingIdx' into its new 'm_BOXES'
        m_map[newx][newy].push_back(movingIdx);
        
        // No error
        indexError = false;
        
        // Break for-loop
        break;
      }
    }

    // If index not found at expected location, prompt a message in the terminal, 
    // then proceed normally
    if(indexError)
    {
      std::cout << "ERROR: Rod #" << movingIdx
                << " not found when moving from (" << oldx << ", " << oldy
                << ") to ("  << newx << ", " << newy << ")"
                << std::endl;
    }
  }

}

void Grid::moveIndex(const int& movingIdx, const std::vector<int>& oldCoords, const std::vector<int>& newCoords)
{
  /**
    'moveIndex': Update the 'm_map' by moving the index 'movingIdx'
                 from the box at ('oldCoords[0]','oldCoords[1]')
                 to the box at ('newCoords[0]','newCoords[1]').
  **/
  
  // If the 'm_map' needs to be updated
  if(oldCoords[0] != newCoords[0] || oldCoords[1] != newCoords[1])
  {
    // Error boolean
    bool indexError = true;
    
    // Get indexes in neighbour 'm_BOXES'
    for(auto& index : m_map[oldCoords[0]][oldCoords[1]])
    {
      // When 'movingIdx' is found
      if(index == movingIdx)
      {
        // Move 'movingIdx' at the end of the array...
        index = m_map[oldCoords[0]][oldCoords[1]].back();
        // ... and erase it
        m_map[oldCoords[0]][oldCoords[1]].pop_back();
        
        // Save 'movingIdx' into its new 'm_BOXES'
        m_map[newCoords[0]][newCoords[1]].push_back(movingIdx);
        
        // No error
        indexError = false;
        
        // Break for-loop
        break;
      }
    }

    // If index not found at expected location, prompt a message in the terminal, 
    // then proceed normally
    if(indexError)
    {
      std::cout << "ERROR: Rod #" << movingIdx
                << " not found when moving from (" << oldCoords[0] << ", " << oldCoords[1]
                << ") to ("  << newCoords[0] << ", " << newCoords[1] << ")"
                << std::endl;
    }
  }
}
