#include <fstream>
#include <iostream>
#include "AnnularCell.h"
#include "montecarloRoutines.h"
#include "analysisRoutines.h"

extern const double PACKING_FRACTION;

int main()
{

/** Example: Check overlaping rod-rod **/

//    Rod rodRef = Rod(6.40671,-9.12292,-0.78476);
//    Rod rod = Rod(7.51301,-9.83874,1.48417);
//    if(rod.isTouchingRod(rodRef))
//        std::cout << "Yes" << std::endl;
//    else
//        std::cout << "No" << std::endl;


/** Example: Annular Cell **/

//    AnnularCell cell = AnnularCell();
//
//    cell.fillAnnularCell();
//
//    do
//    {
//        for(int steps = 0; steps < 1000; steps++)
//        {
//            stepMontecarlo(cell);
//        }
//        cell.fillMissingRods();
//    }while(cell.m_missingRods.size()>0);
//
//    getOrderParameters(cell, "order_parameters.txt");

    return 0;
}
