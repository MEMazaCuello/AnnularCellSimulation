#ifndef ANNULARCELL_H_INCLUDED
#define ANNULARCELL_H_INCLUDED

#include <vector>
#include <iostream>
#include <string>
#include "Rod.h"
#include "Grid.h"

extern const int NUMBER_OF_RODS;

class AnnularCell
{
private:


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

    Rod  getRod(const int& index)   const; // Private ?
    void printRod(const int& index) const;

    void fillAnnularCell(); // Private ?
    void fillMissingRods();
    void fillAnnularCellFromFile(std::string filename, const int& numRodsInFile);
    //void fillFromFile(std::ifstream& filepath);
    bool rodIsTouchingInnerWall(const int& index);
    bool rodIsTouchingInnerWall(const Rod& rod);
    bool rodIsTouchingOuterWall(const int& index);
    bool rodIsTouchingOuterWall(const Rod& rod);
};

#endif // ANNULARCELL_H_INCLUDED
