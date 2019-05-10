#include <cmath>
#include <iostream>
#include "Rod.h"

extern const double WIDTH;
extern const double DIAGONAL;
extern const double HALF_WIDTH;
extern const double HALF_LENGTH;
extern const double HALF_DIAGONAL;
const double Rod::m_HALF_WIDTH    = HALF_WIDTH;
const double Rod::m_HALF_LENGTH   = HALF_LENGTH;
const double Rod::m_HALF_DIAGONAL = HALF_DIAGONAL;

extern const double ALPHA;
extern const double HALF_PI;
extern const double PI;

bool Rod::isTouchingRod(const Rod& refRod)
{
    double xrel = m_xPos-refRod.m_xPos;
    double yrel = m_yPos-refRod.m_yPos;
    double distance = std::sqrt(xrel*xrel + yrel*yrel);
    //std::cout << "Distance = " << distance << std::endl;
    if(distance < WIDTH)
    {
        return true;
    } else if(distance > DIAGONAL)
    {
        return false;
    }

    // phi in -pi/2, pi/2
    double phi = m_angle - refRod.m_angle;
    if(phi > HALF_PI){ phi -= PI;}
    else if(phi < -HALF_PI){ phi += PI;}
    //std::cout << "Phi = " << phi << std::endl;

    // theta in -pi, pi
    double theta = std::atan2(yrel,xrel) - refRod.m_angle;
    if(theta > PI){ theta -= 2.0*PI;}
    else if(theta < -PI){ theta += 2.0*PI;}
    //std::cout << "Theta = " << theta << ", ";

    if(phi < 0.0d)
    {
        phi = -phi;
        if(theta < 0.0d)
        {
            theta = -theta;
        }
        else
        {
            theta = PI-theta;
        }
    }
    else
    {
        if(theta < 0.0d)
        {
            theta = theta+PI;
        }
    }

    double THETA0  = 0.5*phi;
    double THETAm1 = THETA0-ALPHA;
    double THETA1  = THETA0+ALPHA;
    double THETA2  = THETA0+HALF_PI;
    double THETA3  = THETAm1+PI;

    double minDist = DIAGONAL*std::cos(THETA0);

    if(theta < THETAm1)
    {
        //std::cout << "Case 0: " << 0.0 << " < " << theta << " < " << THETAm1 << std::endl;
        minDist *= std::sin(THETA1)/std::sin(phi-theta);
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 0!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }else if(theta < THETA0)
    {
        //std::cout << "Case 1: " << 0.0 << " < " << theta << " < " << THETA1 << std::endl;
        minDist *= std::cos(THETAm1)/std::cos(theta);
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 1!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }else if(theta < THETA1)
    {
        //std::cout << "Case 2: " << THETA1 << " < " << theta << " < " << THETA2 << std::endl;
        minDist *= std::cos(THETAm1)/std::cos(theta-phi);
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 2!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }else if(theta < THETA2)
    {
        //std::cout << "Case 3: " << THETA2 << " < " << theta << " < " << THETA3 << std::endl;
        minDist *= std::sin(THETA1)/std::sin(theta);
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 3!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }else if(theta < THETA3)
    {
        //std::cout << "Case 4: " << THETA3 << " < " << theta << " < " << THETA4 << std::endl;
        minDist *= std::sin(THETA1)/std::sin(theta-phi);
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 4!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }else
    {
        //std::cout << "Case 5: " << THETA4 << " < " << theta << " < " << PI << std::endl;
        minDist *= std::cos(THETAm1)/(-std::cos(theta));
        //if(minDist < 0.0d){ std::cout << "ERROR CASE 5!";}
        //std::cout << "Min distance = " << minDist
        //          << ", distance = " << distance << std::endl;
        if(distance < minDist){return true;}else{return false;}
    }

//    if(distance < minDist)
//    {
//        return true;
//    }else
//    {
//        return false;
//    }

}


