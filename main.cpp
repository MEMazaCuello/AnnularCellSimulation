#include <fstream>
#include <iostream>
#include "AnnularCell.h"
#include "montecarloRoutines.h"
//#include "Grid.h"


extern const double PACKING_FRACTION;

//#include <iostream>
//#include "Rod.h"
//#include <random>

//extern const double QUARTER_PI;
//extern const double HALF_PI;
//extern const double PI;
//extern const double ALPHA;

//inline double excludedArea(const Rod& rod)
//{
//    return 4.0*(2.0*Rod::m_HALF_LENGTH*Rod::m_HALF_WIDTH*(1.0+std::fabs(std::cos(rod.m_angle)))
//                + Rod::m_HALF_DIAGONAL*Rod::m_HALF_DIAGONAL*fabs(std::sin(rod.m_angle)));
//}

//extern const double INNER_RADIUS;
//extern const double OUTER_RADIUS;
//const double R_PLUS_HALF_L = INNER_RADIUS + Rod::m_HALF_LENGTH;
//const double R_PLUS_HALF_W = INNER_RADIUS + Rod::m_HALF_WIDTH;
//const double HALF_D_OVER_R = Rod::m_HALF_DIAGONAL/INNER_RADIUS;
//const double PHI_ONE = std::atan2(Rod::m_HALF_WIDTH, R_PLUS_HALF_L);
//const double PHI_TWO = std::atan2(R_PLUS_HALF_W, Rod::m_HALF_LENGTH);
//bool rodIsTouchingInnerWall(const Rod& rod)
//{
//    static const double INNER_MIN_DIST = INNER_RADIUS + Rod::m_HALF_WIDTH;
//    static const double INNER_MAX_DIST = INNER_RADIUS + Rod::m_HALF_DIAGONAL;
//
//    /* Center too far or too close of inner wall*/
//    double distance = std::sqrt(rod.m_xPos*rod.m_xPos + rod.m_yPos*rod.m_yPos);
//    if(distance > INNER_MAX_DIST){ return false; }
//    else if(distance < INNER_MIN_DIST){ return true; }
//
//    /* Intermediate region */
//    double theta = std::atan2(rod.m_yPos, rod.m_xPos);
//    double phi   = std::abs(rod.m_angle - theta);
//
//    // phi between 0 and PI/2
//    if(phi > PI){ phi -= PI;}
//    if(phi > HALF_PI){ phi = PI - phi;}
//
//    double minDist;
//    if(phi < PHI_ONE)
//    {
//        minDist = R_PLUS_HALF_L/std::cos(phi);
//        if(distance < minDist){return true;}else{return false;}
//    }
//    else if(phi > PHI_TWO)
//    {
//        minDist = R_PLUS_HALF_W/std::sin(phi);
//        if(distance < minDist){return true;}else{return false;}
//    }
//    else
//    {
//        double lambda = std::asin( HALF_D_OVER_R*std::sin(ALPHA-phi) );
//        if(phi < ALPHA)
//        {
//            minDist = (Rod::m_HALF_LENGTH+INNER_RADIUS*std::cos(phi-lambda))/std::cos(phi);
//            if(distance < minDist){return true;}else{return false;}
//        }
//        else
//        {
//            minDist = (Rod::m_HALF_WIDTH+INNER_RADIUS*std::sin(phi-lambda))/std::sin(phi);
//            if(distance < minDist){return true;}else{return false;}
//        }
//    }
//}

//#include <random>
//extern const int NUMBER_OF_RODS;
//std::mt19937_64 g(12345);
//std::uniform_real_distribution<double> dist(-1.0,1.0);


int main()
{

/**_____________________________________________________________________**/

/** Check overlaping rod-rod **/

//    Rod rodRef = Rod(6.40671,-9.12292,-0.78476);
////
////    std::mt19937_64 gen(9547110);
////    std::uniform_real_distribution<double> rndmdist(0.0,1.0);
////
//    Rod rod = Rod(7.51301,-9.83874,1.48417);
////
//////    rod.m_xPos  = -0.4-std::sqrt(0.5);
//////    rod.m_yPos  = -0.4-std::sqrt(0.5);
//////    rod.m_angle = 0.5*HALF_PI;
//////
//    if(rod.isTouchingRod(rodRef))
//        std::cout << "Yes" << std::endl;
//    else
//        std::cout << "No" << std::endl;
//
//    std::ofstream datafile;
//    datafile.open("areas.txt");
//    bool rodsAreTouching;
//    int    trials  = 5e6;
//    double counter;
//    double aux = 5.0*rod.m_HALF_DIAGONAL;
//    rodRef.m_angle = -1.1*HALF_PI;
//    do
//    {
//        rodRef.m_angle += 0.1*HALF_PI;
//        counter = 0.0;
//        for(int i=0; i<trials; i++)
//        {
//            rodRef.m_xPos = (2.0*rndmdist(gen)-1.0)*aux;
//            rodRef.m_yPos = rndmdist(gen)*aux;//(2.0*rndmdist(gen)-1.0)*aux;
//            rodsAreTouching = rodRef.isTouchingRod(rod);
//            if(rodsAreTouching){counter += 1.0;}
//        }
//
//        std::cout << "> Phi = " << rodRef.m_angle/HALF_PI << "*HALF_PI : "
//                  << counter*(4*aux*aux)/trials << " vs " << excludedArea(rodRef) << "\t"
//                  << 50*(counter*(4*aux*aux)/trials - excludedArea(rodRef))/(0.5*excludedArea(rodRef)) << "%" << std::endl;
//
//        datafile << rodRef.m_angle/HALF_PI << " " << counter*(4*aux*aux)/trials << " " << excludedArea(rodRef) << std::endl;
//    }while(rodRef.m_angle < HALF_PI);
//    datafile.close();

/*_____________________________________________________________________*/

/** Check overlaping rod-inner circle **/
//
//    Rod rod = Rod(0.0,0.0,0.0);
//    bool rodsAreTouching;
//    int    trials  = 5e7;
//    double counter = 0;
//    double aux = 5.0*rod.m_HALF_DIAGONAL;
//    std::mt19937_64 gen(495118);
//    std::uniform_real_distribution<double> rndmdist(0.0,1.0);
//
//        for(int i=0; i<trials; i++)
//        {
//            rod.m_xPos = (2.0*rndmdist(gen)-1.0)*aux;
//            rod.m_yPos = (2.0*rndmdist(gen)-1.0)*aux;
//            rodsAreTouching = rodIsTouchingInnerWall(rod);
//            if(rodsAreTouching){counter += 1.0;}
//        }
//
//        std::cout << "Fraction area = " << counter*(4.0*aux*aux)/trials
//                  << " vs " << 4.0*(Rod::m_HALF_LENGTH*Rod::m_HALF_WIDTH+QUARTER_PI*INNER_RADIUS*INNER_RADIUS+(Rod::m_HALF_LENGTH+Rod::m_HALF_WIDTH)*INNER_RADIUS) << "\t" << std::endl;

/*_____________________________________________________________________*/


//    AnnularCell cell = AnnularCell();
//
//    cell.fillAnnularCell();
//
////    for(int i = 0; i<cell.m_NUMBER_OF_PARTICLES; i++)
////    {
////        cell.printRod(i);
////    }
//
//    std::cout << "Packing fraction = " << PACKING_FRACTION << std::endl;
//
//    Rod rod;
//    std::ofstream datafile;
//    datafile.open("final.txt");
//    for(int i=0; i<cell.m_NUMBER_OF_PARTICLES; i++){
//        rod = cell.getRod(i);
//        datafile << i+1 << " " << rod.m_xPos << " " << rod.m_yPos << " " << rod.m_angle << '\r' << '\n';
//    }
//    datafile.close();
//
//    do
//    {
//
//        for(int steps = 0; steps < 1000; steps++)
//        {
//            stepMontecarlo(cell);
//        }
//
//        cell.fillMissingRods();
//    }while(cell.m_missingRods.size()>0);
//
//
//    datafile.open("final2.txt");
//    for(int i=0; i<cell.m_NUMBER_OF_PARTICLES; i++){
//        rod = cell.getRod(i);
//        datafile << i+1 << " " << rod.m_xPos << " " << rod.m_yPos << " " << rod.m_angle << '\r' << '\n';
//
//        if(rod.m_xPos != cell.m_OUTER_RADIUS)
//        {
//            for(int j=0; j<i; j++)
//            {
//                if(rod.isTouchingRod(cell.getRod(j)))
//                {
//                    std::cout << "Error with " << i+1 << " and " << j+1;
//                }
//            }
//        }
//    }
//    datafile.close();
//
//    getOrderParameters(cell);


/*_____________________________________________________________________*/


    AnnularCell cell = AnnularCell();

    std::string filename = "order_param_final_3472_06.txt";
    int numberRods = 3472;
    cell.fillAnnularCellFromFile(filename, numberRods);

//    cell.fillMissingRods();
//    while(cell.m_missingRods.size()>0)
//    {
//        for(int steps = 0; steps < 1000; steps++)
//        {
//            stepMontecarlo(cell);
//        }
//
//        cell.fillMissingRods();
//    }

    /************* Change INITIAL ^^ , FINAL vv , and SEED!!!! >> ********/
    std::string baseOrderParameter = "order_param_3472_07_";
    std::string fnOrderParameter;
    for(int walks = 0; walks < 200; walks++)
    {
        std::cout << "Walk number " << walks << std::endl;
        for(int steps = 0; steps < 10000; steps++)
        {
                stepMontecarlo(cell);
        }
        //fnOrderParameter = baseOrderParameter + std::to_string(walks) + ".txt";

//                fnOrderParameter << baseOrderParameter
//                         << std::setfill( '0' )
//                         << std::setw( 4 ) << ".txt";
        //getOrderParameters(cell, fnOrderParameter);
    }

    fnOrderParameter = "order_param_final_3472_07.txt";
    //getOrderParameters(cell,fnOrderParameter);

    return 0;
}
