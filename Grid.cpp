#include "Grid.h"
#include <cmath>
#include <iostream>
extern const double DIAGONAL;
extern const double OUTER_RADIUS;
const double Grid::m_BOX_LENGTH  = DIAGONAL;
const int    Grid::m_BOXES       = int(std::ceil(2.0d*OUTER_RADIUS/Grid::m_BOX_LENGTH)) + 1;
const double Grid::m_HALF_LENGTH = 0.5d*Grid::m_BOXES*Grid::m_BOX_LENGTH;
const double Grid::m_INV_BOX_LENGTH = 1.0d/Grid::m_BOX_LENGTH;
Grid::Grid()
{
    m_map.resize(m_BOXES);
    for(auto& v : m_map){ v.resize(m_BOXES);}//, std::vector<int> {0});}

    m_NB_XCOORD = m_map;
    m_NB_YCOORD = m_map;
    int grid_boxes_minus_one = m_BOXES - 1;
    for(int i = 0; i < m_BOXES; i++)
    {
        if(i==0)
        {
            for(int j = 1; j < grid_boxes_minus_one; j++)
            {
                m_NB_XCOORD[i][j] = {i,i+1};//{grid_boxes_minus_one,i,i+1};
                m_NB_YCOORD[i][j] = {j-1,j,j+1};
            }

            // Top left corner
            m_NB_XCOORD[i][0] = {i,i+1};//{grid_boxes_minus_one,i,i+1};
            m_NB_YCOORD[i][0] = {0,1};//{grid_boxes_minus_one,0,1};

            // Bottom left corner
            m_NB_XCOORD[i][grid_boxes_minus_one] = {i,i+1};//{grid_boxes_minus_one,i,i+1};
            m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};//{grid_boxes_minus_one-1,grid_boxes_minus_one,0};

        }
        else if(i==grid_boxes_minus_one)
        {
            for(int j = 1; j < grid_boxes_minus_one; j++)
            {
                m_NB_XCOORD[i][j] = {i-1,i};//{i-1,i,0};
                m_NB_YCOORD[i][j] = {j-1,j,j+1};
            }

            // Top right corner
            m_NB_XCOORD[i][0] = {i-1,i};//{i-1,i,0};
            m_NB_YCOORD[i][0] = {0,1};//{grid_boxes_minus_one,0,1};

            // Bottom right corner
            m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i};//{i-1,i,0};
            m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};//{grid_boxes_minus_one-1,grid_boxes_minus_one,0};
        }else
        {
            for(int j = 1; j < grid_boxes_minus_one; j++)
            {
                m_NB_XCOORD[i][j] = {i-1,i,i+1};
                m_NB_YCOORD[i][j] = {j-1,j,j+1};
            }

            // Top row
            m_NB_XCOORD[i][0] = {i-1,i,i+1};
            m_NB_YCOORD[i][0] = {0,1};//{grid_boxes_minus_one,0,1};

            // Bottom row
            m_NB_XCOORD[i][grid_boxes_minus_one] = {i-1,i,i+1};
            m_NB_YCOORD[i][grid_boxes_minus_one] = {grid_boxes_minus_one-1,grid_boxes_minus_one};//{grid_boxes_minus_one-1,grid_boxes_minus_one,0};
        }
    }
}

inline int Grid::getGridX(const Rod& rod) const
{
    return int(std::floor((rod.m_xPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

inline int Grid::getGridY(const Rod& rod) const
{
    return int(std::floor((rod.m_yPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH));
}

std::vector<int> Grid::getGridCoords(const Rod& rod) const
{
    return std::vector<int> {int(std::floor((rod.m_xPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH)),
                             int(std::floor((rod.m_yPos+Grid::m_HALF_LENGTH)*Grid::m_INV_BOX_LENGTH))};
}

void Grid::fillGrid(std::vector<Rod>& bundle)
{
    int index = 0;
    for(const Rod& rod : bundle)
    {
        m_map[getGridX(rod)][getGridY(rod)].push_back(index);
        index++;
    }
}

int Grid::getNumberOfNeighbours(const int& coordx, const int& coordy) const
{
    int number = 0;

    for(const int& idx : m_NB_XCOORD[coordx][coordy])
    {
        for(const int& idy : m_NB_YCOORD[coordx][coordy])
        {
            number += m_map[idx][idy].size();
        }
    }

    return number;
}

std::vector<int> Grid::getNeighbours(const int& coordx, const int& coordy)
{
    m_neighbours.clear();
    m_neighbours.reserve(getNumberOfNeighbours(coordx,coordy));

    for(const int& idx : m_NB_XCOORD[coordx][coordy])
    {
        for(const int& idy : m_NB_YCOORD[coordx][coordy])
        {
            //m_neighbours.insert(m_neighbours.end(),m_map[idx][idy].begin(),m_map[idx][idy].end());
            for(const int& idz : m_map[idx][idy])
            {
                m_neighbours.emplace_back(idz);
            }
        }
    }

    return m_neighbours;
}

std::vector<int> Grid::getNeighbours(const std::vector<int>& coords)
{
    // Note: index of given rod IS included (is neighbour of itself)

    m_neighbours.clear();
    m_neighbours.reserve(getNumberOfNeighbours(coords[0],coords[1]));

    for(const int& idx : m_NB_XCOORD[coords[0]][coords[1]])
    {
        for(const int& idy : m_NB_YCOORD[coords[0]][coords[1]])
        {
            //m_neighbours.insert(m_neighbours.end(),m_map[idx][idy].begin(),m_map[idx][idy].end());
            for(const int& idz : m_map[idx][idy])
            {
                m_neighbours.emplace_back(idz);
            }
        }
    }

    return m_neighbours;
}

void Grid::moveIndex(const int& movingIdx, const int& oldx, const int& oldy, const int& newx, const int& newy)
{
    if(oldx != newx || oldy != newy)
    {
        bool indexError = true;
        for(auto& index : m_map[oldx][oldy]) // This can be improved using iterators
        {
            if(index == movingIdx)
            {
                index = m_map[oldx][oldy].back();
                m_map[oldx][oldy].pop_back();
                m_map[newx][newy].push_back(movingIdx);
                indexError = false;
                break;
            }
        }

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
    if(oldCoords[0] != newCoords[0] || oldCoords[1] != newCoords[1])
    {
        bool indexError = true;
        for(auto& index : m_map[oldCoords[0]][oldCoords[1]]) // This can be improved using iterators
        {
            if(index == movingIdx)
            {
                index = m_map[oldCoords[0]][oldCoords[1]].back();
                m_map[oldCoords[0]][oldCoords[1]].pop_back();
                m_map[newCoords[0]][newCoords[1]].push_back(movingIdx);
                indexError = false;
                break;
            }
        }

        if(indexError)
        {
            std::cout << "ERROR: Rod #" << movingIdx
                      << " not found when moving from (" << oldCoords[0] << ", " << oldCoords[1]
                      << ") to ("  << newCoords[0] << ", " << newCoords[1] << ")"
                      << std::endl;
        }
    }
}
