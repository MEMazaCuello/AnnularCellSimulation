#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Rod.h"
#include <vector>

class Grid
{
public:
    static const double m_BOX_LENGTH;
    static const double m_INV_BOX_LENGTH;
    static const int    m_BOXES;
    static const double m_HALF_LENGTH;
    std::vector< std::vector< std::vector<int> > > m_map;
    std::vector< std::vector< std::vector<int> > > m_NB_XCOORD;
    std::vector< std::vector< std::vector<int> > > m_NB_YCOORD;
    std::vector<int>    m_neighbours;

    Grid();

    void fillGrid(std::vector<Rod>& bundle);
    int getGridX(const Rod& rod) const;
    int getGridY(const Rod& rod) const;
    std::vector<int> getGridCoords(const Rod& rod) const;
    int  getNumberOfNeighbours(const int& coordx, const int& coordy) const;
    std::vector<int> getNeighbours(const int& coordx, const int& coordy);
    std::vector<int> getNeighbours(const std::vector<int>& coords);
    void moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy);
    void moveIndex(const int& movingIdx, const std::vector<int>& oldCoords, const std::vector<int>& newCoords);
};

#endif // GRID_H_INCLUDED
