/**
  * "Grid.cpp":
  * --------------------
  * Implementation of the Grid class.
  * For declaration details, see "Grid.cpp".
  *
  * Needs: "Grid.h".
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#include <cmath>

#include "Grid.h"

// Parameters defined in "parameters.cpp" _____________________________
extern const double DIAGONAL;
extern const double OUTER_RADIUS;

// Static variables ___________________________________________________
const double Grid::m_BOX_LENGTH  = DIAGONAL;
const int    Grid::m_BOXES       = int(std::ceil(2.0d*OUTER_RADIUS/Grid::m_BOX_LENGTH)) + 1;
const double Grid::m_HALF_LENGTH = 0.5d*Grid::m_BOXES*Grid::m_BOX_LENGTH;
const double Grid::m_INV_BOX_LENGTH = 1.0d/Grid::m_BOX_LENGTH;

// Methods ____________________________________________________________
Grid::Grid()
{
  m_map.resize(m_BOXES);
  for(auto& v : m_map){ v.resize(m_BOXES);}

  m_NB_XCOORD = m_map;
  m_NB_YCOORD = m_map;

  int grid_boxes_minus_one = m_BOXES - 1;
  for (int i = 0; i < m_BOXES; i++)
  {
    if (i == 0) // Left-most column
    {
      for (int j = 1; j < grid_boxes_minus_one; j++)
      {
        m_NB_XCOORD[i][j] = {i,i+1};
        m_NB_YCOORD[i][j] = {j-1,j,j+1};
      }

      // Top left corner
      m_NB_XCOORD[i][0] = {i,i+1};
      m_NB_YCOORD[i][0] = {0,1};

      // Bottom left corner
      m_NB_XCOORD[i][grid_boxes_minus_one] = {i,i+1};
      m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
    else if (i == grid_boxes_minus_one) // Right-most column
    {
      for (int j = 1; j < grid_boxes_minus_one; j++)
      {
        m_NB_XCOORD[i][j] = {i-1,i};
        m_NB_YCOORD[i][j] = {j-1,j,j+1};
      }

      // Top right corner
      m_NB_XCOORD[i][0] = {i-1,i};
      m_NB_YCOORD[i][0] = {0,1};

      // Bottom right corner
      m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i};
      m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
    else
    {
      for (int j = 1; j < grid_boxes_minus_one; j++)
      {
        m_NB_XCOORD[i][j] = {i-1,i,i+1};
        m_NB_YCOORD[i][j] = {j-1,j,j+1};
      }

      // Top row
      m_NB_XCOORD[i][0] = {i-1,i,i+1};
      m_NB_YCOORD[i][0] = {0,1};

      // Bottom row
      m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i,i+1};
      m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};
    }
  }
}

int Grid::getGridIdx(const double& val) const
{
  return int(std::floor((val+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

int Grid::getGridX(const Rod& rod) const
{
  return int(std::floor((rod.m_xPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

int Grid::getGridY(const Rod& rod) const
{
  return int(std::floor((rod.m_yPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

std::vector<int> Grid::getCoords(const Rod& rod) const
{
  return std::vector<int> { getGridX(rod) , getGridY(rod) };
}

void Grid::fill(std::vector<Rod>& bundle)
{
  int index = 0;
  for (const Rod& rod : bundle)
  {
    m_map[getGridX(rod)][getGridY(rod)].emplace_front(index);
    index++;
  }
}

std::vector<int> Grid::getBox(const int& coordx, const int& coordy)
{
  std::vector<int> indexes;

  for (const int& idz : m_map[coordx][coordy])
  {
    indexes.push_back(idz);
  }

  return indexes;
}

std::vector<int> Grid::getNeighbors(const int& coordx, const int& coordy)
{
  // Note: index of given rod IS included (is neighbor of itself)

  std::vector<int> neighbors;

  for (const int& idx : m_NB_XCOORD[coordx][coordy])
  {
    for (const int& idy : m_NB_YCOORD[coordx][coordy])
    {
      for (const int& idz : m_map[idx][idy])
      {
        neighbors.push_back(idz);
      }
    }
  }

  return neighbors;
}

std::vector<int> Grid::getNeighbors(const std::vector<int>& coords)
{
  return getNeighbors(coords[0],coords[1]);
}

void Grid::setNeighbors(const int& coordx, const int& coordy)
{
  // Note: m_neighbors is modified.
  // Note: index of given rod IS included (is neighbor of itself)

  m_neighbors.clear();

  for (const int& idx : m_NB_XCOORD[coordx][coordy])
  {
    for (const int& idy : m_NB_YCOORD[coordx][coordy])
    {
      for (const int& idz : m_map[idx][idy])
      {
        m_neighbors.push_back(idz);
      }
    }
  }
}

void Grid::setNeighbors(const std::vector<int>& coords)
{
  setNeighbors(coords[0],coords[1]);
}

void Grid::moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy)
{
  if (oldx != newx || oldy != newy)
  {
    m_map[oldx][oldy].remove(movingIdx);
    m_map[newx][newy].emplace_front(movingIdx);
  }
}

void Grid::moveIndex(const int& movingIdx, const std::vector<int>& oldCoords, const std::vector<int>& newCoords)
{
  moveIndex(movingIdx,oldCoords[0],oldCoords[1],newCoords[0],newCoords[1]);
}
