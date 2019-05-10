#include "AnnularCell.h"
#include "montecarloRoutines.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

// Auxiliary parameters from "parameters.cpp"
extern const double HALF_PI;
extern const double PI;
//extern const double QUARTER_PI;
extern const double ALPHA;

// Random number generator
std::mt19937_64 gen(10082273);
std::uniform_real_distribution<double> rndmdist(-1.0,1.0);

// Load dimensions from "parameters.cpp"
extern const double INNER_RADIUS;
extern const double OUTER_RADIUS;
const double AnnularCell::m_INNER_RADIUS = INNER_RADIUS;
const double AnnularCell::m_OUTER_RADIUS = OUTER_RADIUS;

/* Get rod from bundle */
Rod  AnnularCell::getRod(const int& index) const  // Private ?
{
    return m_bundle[index];
}

/* Print information about rod to terminal */
void  AnnularCell::printRod(const int& index) const // Private ?
{
    std::cout << "Rod #" << index << ": ( " << m_bundle[index].m_xPos << " , "
                                            << m_bundle[index].m_yPos << " , "
                                            << m_bundle[index].m_angle << " ) " << std::endl;
}


/* Is rectangle touching the circle in center of the cell? */
// Auxiliary global parameters
const double R_PLUS_HALF_L = INNER_RADIUS + Rod::m_HALF_LENGTH;
const double R_PLUS_HALF_W = INNER_RADIUS + Rod::m_HALF_WIDTH;
const double HALF_D_OVER_R = Rod::m_HALF_DIAGONAL/INNER_RADIUS;
const double PHI_ONE = std::atan2(Rod::m_HALF_WIDTH, R_PLUS_HALF_L);
const double PHI_TWO = std::atan2(R_PLUS_HALF_W, Rod::m_HALF_LENGTH);

// Index version
bool AnnularCell::rodIsTouchingInnerWall(const int& index)
{
    static const double INNER_MIN_DIST = m_INNER_RADIUS + Rod::m_HALF_WIDTH;
    static const double INNER_MAX_DIST = m_INNER_RADIUS + Rod::m_HALF_DIAGONAL;

    if(index < m_NUMBER_OF_PARTICLES)
    {
        // Get rod
        m_aux_rod = getRod(index);

        /* Center too far or too close of inner wall*/
        double distance = std::sqrt(m_aux_rod.m_xPos*m_aux_rod.m_xPos + m_aux_rod.m_yPos*m_aux_rod.m_yPos);
        if(distance > INNER_MAX_DIST){ return false; }
        if(distance < INNER_MIN_DIST){ return  true; }

        /* Intermediate region */
        double theta = std::atan2(m_aux_rod.m_yPos, m_aux_rod.m_xPos);
        double phi   = std::abs(m_aux_rod.m_angle - theta);

        // phi between 0 and PI/2
        if(phi > PI){ phi -= PI;}
        if(phi > HALF_PI){ phi = PI - phi;}

        double minDist;
        if(phi < PHI_ONE)
        {
            minDist = R_PLUS_HALF_L/std::cos(phi);
            if(distance < minDist){return true;}else{return false;}
        }
        else if(phi > PHI_TWO)
        {
            minDist = R_PLUS_HALF_W/std::sin(phi);
            if(distance < minDist){return true;}else{return false;}
        }
        else
        {
            double lambda = std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
            if(phi < ALPHA)
            {
                minDist = (Rod::m_HALF_LENGTH+m_INNER_RADIUS*std::cos(phi-lambda))/std::cos(phi);
                if(distance < minDist){return true;}else{return false;}
            }
            else
            {
                minDist = (Rod::m_HALF_WIDTH+m_INNER_RADIUS*std::sin(phi-lambda))/std::sin(phi);
                if(distance < minDist){return true;}else{return false;}
            }
        }
    }
    else
    {
        return true;
    }
}

// Rod version
bool AnnularCell::rodIsTouchingInnerWall(const Rod& rod)
{
    static const double INNER_MIN_DIST = m_INNER_RADIUS + Rod::m_HALF_WIDTH;
    static const double INNER_MAX_DIST = m_INNER_RADIUS + Rod::m_HALF_DIAGONAL;

    /* Center too far or too close of inner wall*/
    double distance = std::sqrt(rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos);
    if(distance > INNER_MAX_DIST){ return false; }
    if(distance < INNER_MIN_DIST){ return  true; }

    /* Intermediate region */
    double theta = std::atan2(rod.m_yPos, rod.m_xPos);
    double phi   = std::abs(rod.m_angle - theta);

    // phi between 0 and PI/2
    if(phi > PI){ phi -= PI;}
    if(phi > HALF_PI){ phi = PI - phi;}

    double minDist;
    if(phi < PHI_ONE)
    {
        minDist = R_PLUS_HALF_L/std::cos(phi);
        if(distance < minDist){return true;}else{return false;}
    }
    else if(phi > PHI_TWO)
    {
        minDist = R_PLUS_HALF_W/std::sin(phi);
        if(distance < minDist){return true;}else{return false;}
    }
    else
    {
        double lambda = std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
        if(phi < ALPHA)
        {
            minDist = (Rod::m_HALF_LENGTH+m_INNER_RADIUS*std::cos(phi-lambda))/std::cos(phi);
            if(distance < minDist){return true;}else{return false;}
        }
        else
        {
            minDist = (Rod::m_HALF_WIDTH+m_INNER_RADIUS*std::sin(phi-lambda))/std::sin(phi);
            if(distance < minDist){return true;}else{return false;}
        }
    }
}

/* Is rectangle touching the external wall of the cell? */
// Index version
bool AnnularCell::rodIsTouchingOuterWall(const int& index)
{
    // OUTER RADIUS GREATER THAN HALF_DIAGONAL!!
    // UNKNOWN BEHAVIOUR OTHERWISE
    static const double OUTER_MIN_DIST = std::sqrt(m_OUTER_RADIUS*m_OUTER_RADIUS-Rod::m_HALF_LENGTH*Rod::m_HALF_LENGTH) - Rod::m_HALF_WIDTH;
    static const double OUTER_MAX_DIST = m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL;

    if(index < m_NUMBER_OF_PARTICLES)
    {
        m_aux_rod = getRod(index);

        /* Center too far / too close of inner wall*/
        double distance = std::sqrt(m_aux_rod.m_xPos*m_aux_rod.m_xPos + m_aux_rod.m_yPos*m_aux_rod.m_yPos);
        if(distance > OUTER_MIN_DIST){ return true; }
        if(distance < OUTER_MAX_DIST){ return false; }

        /* Intermediate region */
        double theta = std::atan2(m_aux_rod.m_yPos, m_aux_rod.m_xPos);

        double phi = m_aux_rod.m_angle - theta;
        if(phi < -HALF_PI)
        {
            phi = phi + PI;
        }else if(phi > HALF_PI)
        {
            phi = phi - PI;
        }

        phi = std::cos(ALPHA-std::fabs(phi));

        if(distance > (std::sqrt(m_OUTER_RADIUS*m_OUTER_RADIUS-Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL*(1.0-phi*phi))-Rod::m_HALF_DIAGONAL*phi ))
        {
            return true;
        }else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}

// Rod version
bool AnnularCell::rodIsTouchingOuterWall(const Rod& rod)
{

    static const double OUTER_MIN_DIST = m_OUTER_RADIUS - Rod::m_HALF_WIDTH;
    static const double OUTER_MAX_DIST = m_OUTER_RADIUS - Rod::m_HALF_DIAGONAL;

    /* Center too far / too close of inner wall*/
    double distance = std::sqrt(rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos);
    if(distance > OUTER_MIN_DIST){ return true; }
    if(distance < OUTER_MAX_DIST){ return false; }

    /* Intermediate region */
    double theta = std::atan2(rod.m_yPos, rod.m_xPos);

    double phi = rod.m_angle - theta;
    if(phi < -HALF_PI)
    {
        phi = phi + PI;
    }else if(phi > HALF_PI)
    {
        phi = phi - PI;
    }

    phi = std::cos(ALPHA-std::fabs(phi));

    if(distance > (std::sqrt(m_OUTER_RADIUS*m_OUTER_RADIUS-Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL*(1.0-phi*phi))-Rod::m_HALF_DIAGONAL*phi ))
    {
        return true;
    }else
    {
        return false;
    }
}

/* Fill annular cell with rectangles at random positions */
void AnnularCell::fillAnnularCell() // Private ?
{
    bool rodsAreTouching;
    for(int i=0; i < m_NUMBER_OF_PARTICLES; i++)
    {
        rodsAreTouching = true;
        for(int trials = 0; trials < 1000; trials++)
        {
            /* Global approach */
            m_aux_rod.m_xPos  = rndmdist(gen)*m_OUTER_RADIUS;
            m_aux_rod.m_yPos  = rndmdist(gen)*m_OUTER_RADIUS;
            m_aux_rod.m_angle = rndmdist(gen)*HALF_PI;

            /* Check if position is valid */
            if(rodIsTouchingInnerWall(m_aux_rod)||rodIsTouchingOuterWall(m_aux_rod))
            {
                continue;
            }
            else
            {
                rodsAreTouching = false;
                for(int j = 0; j < i; j++)
                {
                    rodsAreTouching = m_aux_rod.isTouchingRod(getRod(j));
                    if(rodsAreTouching){ break;}
                }
            }

            /* Save rod */
            if(!rodsAreTouching)
            {
                m_bundle.emplace_back(m_aux_rod);
                break;
            }
        }

        /* Create vector or the indexes that could not be inserted */
        if(rodsAreTouching)
        {
            m_missingRods.push_back(i);
            std::cout << "Rod #" << i << " not included" << std::endl;
            m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
        }
    }

    /* Fill grid of indexes to know roughly were each rod is*/
    m_grid.fillGrid(m_bundle);

    std::cout << std::endl << "\t CELL FILLED" << std::endl;
}

/* Try to fill rods that were missing */
void AnnularCell::fillMissingRods()
{
    bool rodsAreTouching;
    std::vector<int> newMissing;
    for(int i : m_missingRods)
    {
        rodsAreTouching = true;
        for(int trials = 0; trials < 500; trials++)
        {
            m_aux_rod.m_xPos  = rndmdist(gen)*m_OUTER_RADIUS;
            m_aux_rod.m_yPos  = rndmdist(gen)*m_OUTER_RADIUS;

            /* For each position try several different angles*/
            for(int n = 0; n < 100; n++)
            {
                m_aux_rod.m_angle = rndmdist(gen)*HALF_PI;

                /* Check if position is valid */
                if(rodIsTouchingInnerWall(m_aux_rod)||rodIsTouchingOuterWall(m_aux_rod))
                {
                    continue;
                }
                else
                {
                    rodsAreTouching = false;
                    for(int j = 0; j < NUMBER_OF_RODS; j++)
                    {
                        if(j!=i)
                        {
                            rodsAreTouching = m_aux_rod.isTouchingRod(getRod(j));
                            if(rodsAreTouching){ break;}
                        }
                    }
                }

                /* Save rod */
                if(!rodsAreTouching)
                {
                    m_grid.moveIndex(i,m_grid.getGridCoords(m_bundle[i]),m_grid.getGridCoords(m_aux_rod));
                    m_bundle[i] = m_aux_rod;
                    break;
                }
            }
            if(!rodsAreTouching){break;}

        }

        if(rodsAreTouching)
        {
            newMissing.push_back(i);
        }
    }

    m_missingRods.clear();
    m_missingRods = newMissing;
    std::cout << m_missingRods.size() << " rods missing" << std::endl;
}

/*
    Fill rods from file with cols: index, X, Y, PHI, q1, q2, q3, q4
    were q_ are order parameters of the liquid crystal
*/
void AnnularCell::fillAnnularCellFromFile(std::string filepath, const int& numRodsInFile)
{
    /** ASSUMES numRodsInFile <= NUMBER_OF_RODS **/
    double dummy;
    std::ifstream savedconfiguration;
    savedconfiguration.open(filepath);
    for(int i=0; i < numRodsInFile; i++)
    {
        savedconfiguration >> dummy >> m_aux_rod.m_xPos >> m_aux_rod.m_yPos >> m_aux_rod.m_angle
                           >> dummy >> dummy >> dummy >> dummy; // From order parameters files
        m_bundle.emplace_back(m_aux_rod);
    }
    savedconfiguration.close();

    for(int i=numRodsInFile; i<NUMBER_OF_RODS; i++)
    {
        m_missingRods.push_back(i);
        m_bundle.emplace_back(m_OUTER_RADIUS,m_OUTER_RADIUS,-ALPHA);
    }

    m_grid.fillGrid(m_bundle);
}
