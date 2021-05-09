/**
  * "Grid.h":
  * --------------------
  * Declaration of the Grid class.
  * For implementation details, see "Grid.cpp".
  *
  * Needs: "Rod.h".
  *
  * --------------------
  * A "Grid" is a division of the plane in a "m_map"
  * consisting of a square lattice of
  *   ("m_BOXES" x "m_BOXES")
  * squares of side
  *   "m_BOX_LENGTH".
  * Each box contains the indexes of those rods
  * whose centers are located within the box.
  * The "m_map" is complemented with "neighboring boxes" maps
  * both for the vertical and horizontal direction, used for
  * retrieving the indexes on a given box and its adjacent boxes.
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <forward_list>

#include "Rod.h"

class Grid
{
public:

  static const int    m_BOXES;
  static const double m_BOX_LENGTH;
  static const double m_INV_BOX_LENGTH;
  static const double m_HALF_LENGTH;

  std::vector< std::vector< std::forward_list<int> > > m_map;

  /**
    * Coordinates in the map of neighboring boxes a given box.
    */
  std::vector< std::vector< std::forward_list<int> > > m_NB_XCOORD;
  std::vector< std::vector< std::forward_list<int> > > m_NB_YCOORD;

  /**
    * Vector of all the indexes inside a box and its 8 neighboring boxes.
    */
  std::vector<int> m_neighbors;

  Grid();

  void fill(std::vector<Rod>& bundle);

  int getGridIdx(const double& val) const;
  int getGridX(const Rod& rod) const;
  int getGridY(const Rod& rod) const;

  std::vector<int> getCoords(const Rod& rod) const;
  std::vector<int> getBox(const int& coordx, const int& coordy);

  /**
    * Gets vector of all the indexes inside a box and its 8 neighboring boxes.
    * Note: each rod is a neighbor of itself.
    */
  std::vector<int> getNeighbors(const int& coordx, const int& coordy);
  std::vector<int> getNeighbors(const std::vector<int>& coords);

  /**
    * Set m_neighbors to have all the indexes inside a box and its 8 neighboring boxes.
    * Note: m_neighbors is modified.
    * Note: each rod is a neighbor of itself.
    */
  void setNeighbors(const int& coordx, const int& coordy);
  void setNeighbors(const std::vector<int>& coords);

  /**
    * Update the box to which a given index belongs to.
    */
  void moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy);
  void moveIndex(const int& movingIdx, const std::vector<int>& oldCoords, const std::vector<int>& newCoords);
};

#endif // GRID_H_INCLUDED
