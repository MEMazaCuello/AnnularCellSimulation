#include "montecarloRoutines.h"
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

std::mt19937_64 g(98710536);
std::uniform_real_distribution<double> dist(-1.0,1.0);

/* External parameters from 'parameters.cpp' */
extern const int    NUMBER_OF_RODS;
extern const double LENGTH;
extern const double HALF_PI;
extern const double PI;
extern const double OUTER_RADIUS;
extern const double HALF_DIAGONAL;

void stepMontecarlo(AnnularCell& cell)
{
    static double DELTA_SPACE = 0.1*LENGTH;
    static double DELTA_ANGLE = 0.1*HALF_PI;
    static const double MAX_RADIUS = OUTER_RADIUS-HALF_DIAGONAL;
    static const double INV_NUMBER_OF_RODS = 1.0d/NUMBER_OF_RODS;
    static std::vector<int> indexes(NUMBER_OF_RODS);

    bool validPosition;
    int  numSuccess = 0;

    // Generate random index permutation
    std::iota(indexes.begin(), indexes.end(), 0);
    std::shuffle(indexes.begin(), indexes.end(), g);
    //for(int i = 0; i < NUMBER_OF_RODS; i++)
    for(int i : indexes)
    {
        cell.m_aux_rod = cell.getRod(i);

        cell.m_aux_rod.m_xPos  += dist(g)*DELTA_SPACE;
        if(cell.m_aux_rod.m_xPos > OUTER_RADIUS){ cell.m_aux_rod.m_xPos = MAX_RADIUS;}
        else if(cell.m_aux_rod.m_xPos < - OUTER_RADIUS){ cell.m_aux_rod.m_xPos = -MAX_RADIUS;}
        cell.m_aux_rod.m_yPos  += dist(g)*DELTA_SPACE;
        if(cell.m_aux_rod.m_yPos > OUTER_RADIUS){ cell.m_aux_rod.m_yPos = MAX_RADIUS;}
        else if(cell.m_aux_rod.m_yPos < - OUTER_RADIUS){ cell.m_aux_rod.m_yPos = -MAX_RADIUS;}
        cell.m_aux_rod.m_angle += dist(g)*DELTA_ANGLE;
        while(cell.m_aux_rod.m_angle > HALF_PI)
        {
            cell.m_aux_rod.m_angle -= PI;
        }
        while(cell.m_aux_rod.m_angle < -HALF_PI)
        {
            cell.m_aux_rod.m_angle += PI;
        }

        cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_aux_rod));

        validPosition = true;
        if(cell.rodIsTouchingInnerWall(cell.m_aux_rod)||cell.rodIsTouchingOuterWall(cell.m_aux_rod))
        {
            validPosition = false;
        }else{
            for(int idx : cell.m_grid.m_neighbours)
            {
                if(idx != i)
                {
                    if(cell.m_aux_rod.isTouchingRod(cell.getRod(idx)))
                    {
                        validPosition = false;
                        break;
                    }
                }
            }
        }

        if(validPosition)
        {
            cell.m_grid.moveIndex(i,cell.m_grid.getGridCoords(cell.m_bundle[i]), cell.m_grid.getGridCoords(cell.m_aux_rod));
            cell.m_bundle[i]= cell.m_aux_rod;
            numSuccess++;
        }
    }

    /* Correction of deltas to generate 50% acceptance */
    DELTA_SPACE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);
    DELTA_ANGLE *= (0.5d+double(numSuccess)*INV_NUMBER_OF_RODS);

}

void getOrderParameters(AnnularCell& cell){

    static double scale = 2.0d*PI/(1.2d*LENGTH);
    static std::ofstream Qdat;

    //int sec = (int)std::floor(steps/1000000);
    std::ostringstream filename;
    filename << "doc";
    //filename << std::setfill( '0' ) << std::setw( 3 ) << sec;
    filename << ".txt";

    Qdat.open( filename.str().c_str() );

    double tilt, zeta;
    double sum_cos2, sum_sin2;
    double sum_cos4, sum_sin4;
    double sum_cosS, sum_sinS;

    for(int i=0; i<NUMBER_OF_RODS; i++){

        sum_cos2 = 0.0d;
        sum_sin2 = 0.0d;
        cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_bundle[i]));
        for(int index : cell.m_grid.m_neighbours)
        {
            sum_cos2 += std::cos(2.0d*cell.m_bundle[index].m_angle);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin2 += std::sin(2.0d*cell.m_bundle[index].m_angle);      // index INCLUDES current 'i' rod
        }
        tilt = 0.5d*std::atan2(sum_sin2,sum_cos2);

        sum_cos2 = 0.0d; sum_sin2 = 0.0d;
        sum_cos4 = 0.0d; sum_sin4 = 0.0d;
        sum_cosS = 0.0d; sum_sinS = 0.0d;

        //cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_aux_rod));
        for(int index : cell.m_grid.m_neighbours)
        {
            zeta = 2.0d*(cell.m_bundle[index].m_angle-tilt);
            sum_cos2 += std::cos(zeta);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin2 += std::sin(zeta);      // index INCLUDES current 'i' rod

            sum_cos4 += std::cos(2.0d*zeta);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin4 += std::sin(2.0d*zeta);      // index INCLUDES current 'i' rod

            zeta = scale*(std::cos(tilt)*(cell.m_bundle[index].m_xPos-cell.m_bundle[i].m_xPos) + std::sin(tilt)*(cell.m_bundle[index].m_yPos-cell.m_bundle[i].m_yPos));
            sum_cosS += std::cos(zeta); // cos( 2pi/(1.2L) * zeta)
            sum_sinS += std::sin(zeta);
        }

        zeta = 1.0d/cell.m_grid.m_neighbours.size();
        Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " "
             << cell.m_bundle[i].m_angle << " " << tilt << " "                           // Save configuration (and tilt (eigen)angle)
             << std::sqrt(sum_cos2*sum_cos2 + sum_sin2*sum_sin2)*zeta << " "             // Nematic order parameter Q1
             << std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)*zeta << " "             // Tetratic order parameter Q2
             << std::sqrt(sum_cosS*sum_cosS + sum_sinS*sum_sinS)*zeta <<                 // Smectic order parameter Qs
             '\r' << '\n';

    }
    Qdat.close();

}

void getOrderParameters(AnnularCell& cell, std::string& filename){

    static double scale = 2.0d*PI/(1.2d*LENGTH);
    static std::ofstream Qdat;

//    //int sec = (int)std::floor(steps/1000000);
//    std::ostringstream filename;
//    filename << "doc";
//    //filename << std::setfill( '0' ) << std::setw( 3 ) << sec;
//    filename << ".txt";

    Qdat.open(filename);

    double tilt, zeta;
    double sum_cos2, sum_sin2;
    double sum_cos4, sum_sin4;
    double sum_cosS, sum_sinS;

    for(int i=0; i<NUMBER_OF_RODS; i++){

        sum_cos2 = 0.0d;
        sum_sin2 = 0.0d;
        cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_bundle[i]));
        for(int index : cell.m_grid.m_neighbours)
        {
            sum_cos2 += std::cos(2.0d*cell.m_bundle[index].m_angle);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin2 += std::sin(2.0d*cell.m_bundle[index].m_angle);      // index INCLUDES current 'i' rod
        }
        tilt = 0.5d*std::atan2(sum_sin2,sum_cos2);

        sum_cos2 = 0.0d; sum_sin2 = 0.0d;
        sum_cos4 = 0.0d; sum_sin4 = 0.0d;
        sum_cosS = 0.0d; sum_sinS = 0.0d;

        //cell.m_grid.m_neighbours = cell.m_grid.getNeighbours(cell.m_grid.getGridCoords(cell.m_aux_rod));
        for(int index : cell.m_grid.m_neighbours)
        {
            zeta = 2.0d*(cell.m_bundle[index].m_angle-tilt);
            sum_cos2 += std::cos(zeta);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin2 += std::sin(zeta);      // index INCLUDES current 'i' rod

            sum_cos4 += std::cos(2.0d*zeta);      // Mean value of cos(2phi) and sin(2phi) to get the nematic (eigen)direction
            sum_sin4 += std::sin(2.0d*zeta);      // index INCLUDES current 'i' rod

            zeta = scale*(std::cos(tilt)*(cell.m_bundle[index].m_xPos-cell.m_bundle[i].m_xPos) + std::sin(tilt)*(cell.m_bundle[index].m_yPos-cell.m_bundle[i].m_yPos));
            sum_cosS += std::cos(zeta); // cos( 2pi/(1.2L) * zeta)
            sum_sinS += std::sin(zeta);
        }

        zeta = 1.0d/cell.m_grid.m_neighbours.size();
        Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " "
             << cell.m_bundle[i].m_angle << " " << tilt << " "                           // Save configuration (and tilt (eigen)angle)
             << std::sqrt(sum_cos2*sum_cos2 + sum_sin2*sum_sin2)*zeta << " "             // Nematic order parameter Q1
             << std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)*zeta << " "             // Tetratic order parameter Q2
             << std::sqrt(sum_cosS*sum_cosS + sum_sinS*sum_sinS)*zeta <<                 // Smectic order parameter Qs
             '\r' << '\n';

    }
    Qdat.close();

}

