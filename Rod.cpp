/**** 
  'Rod.cpp':
  --------------------
  Implementation of Rod struct.
  For declaration details, see 'Rod.h'.
  
  Needs: 'parameters.cpp'
  
  --------------------
  A 'Rod' is represented by a rectangle with
    width  = 2*'HALF_WIDTH' and
    length = 2*'HALF_LENGTH'
  at the 2D cartesian position of its center
    (m_xPos, m_yPos)
  with orientation
    m_angle
  defined as the angle between the 'LENGTH'- and the OX-axis.
  
  Note: the file 'parameters.cpp' is needed to specify the parameters
        'HALF_WIDTH', 'HALF_LENGTH', 'HALF_DIAGONAL',
        'WIDTH', 'LENGTH', 'DIAGONAL', 
        'ALPHA', 'HALF_PI' and 'PI'.
        
  --------------------
  Last modified: 2019-05-12 
  By: M. E. Maza Cuello
****/

#include <cmath>
#include <iostream>
#include "Rod.h"

/** PARAMETERS OBTAINED FROM 'parameters.cpp' **/
// Lengths
extern const double WIDTH;
extern const double DIAGONAL;
extern const double HALF_WIDTH;
extern const double HALF_LENGTH;
extern const double HALF_DIAGONAL;
// Angles
extern const double PI;
extern const double HALF_PI;
extern const double ALPHA;

/** 'Rod' IMPLEMENTATION **/
/* Static variables */
// Length variables
const double Rod::m_HALF_WIDTH    = HALF_WIDTH;
const double Rod::m_HALF_LENGTH   = HALF_LENGTH;
const double Rod::m_HALF_DIAGONAL = HALF_DIAGONAL;

/* Methods */

bool Rod::isTouchingRod(const Rod& otherRod)
{
  /** 
    'isTouchingRod': Return 'true' if rod is overlapping 'otherRod', 'false' if not.
                     Assumes 'm_angle' of both rods to be in the interval [-HALF_PI, HALF_PI].
  **/
  
  /* First criterion: rods too far/close to eachother to overlap */
  // Relative cartesian position between this and 'refRod'
  double xrel = m_xPos-otherRod.m_xPos;
  double yrel = m_yPos-otherRod.m_yPos;
  
  // Cartesian distance between this and 'refRod' centers
  double distance = std::sqrt(xrel*xrel + yrel*yrel);

  // Check first criterion
  if(distance < WIDTH)
  {
      return true;
  } else if(distance > DIAGONAL)
  {
      return false;
  }

  /* Second criterion: 'distance' less than minimum distance 'minDist', 
                      obtained via an analytical expression which is function
                      of relatives angles 'phi' and 'theta' (defined below)   */
  
  // Relative orientation between rods orientations: 'phi'
  double phi = m_angle - otherRod.m_angle;
  
  // 'phi' in interval [-pi/2, pi/2]
  if(phi > HALF_PI){ phi -= PI;}
  else if(phi < -HALF_PI){ phi += PI;}
  
  // Relative angle between rod centers (with respect to otherRod orientation): 'theta'
  double theta = std::atan2(yrel,xrel) - otherRod.m_angle;
  
  // 'theta' in interval [-pi, pi]
  if(theta > PI){ theta -= 2.0*PI;}
  else if(theta < -PI){ theta += 2.0*PI;}

  // Rotate / reflect system to work in the first and second cartesian quadrants 
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

  // Auxiliary angles: analytical frontiers for different cheking regions
  double THETA0  = 0.5*phi;
  double THETAm1 = THETA0-ALPHA;
  double THETA1  = THETA0+ALPHA;
  double THETA2  = THETA0+HALF_PI;
  double THETA3  = THETAm1+PI;
  
  // Common factor of 'minDist'
  double minDist = DIAGONAL*std::cos(THETA0);

  // Compute specific 'minDist' 
  if(theta < THETAm1)
  {
      // First region: 'theta' in [-pi, phi/2-ALPHA]
      minDist *= std::sin(THETA1)/std::sin(phi-theta);
  }else if(theta < THETA0)
  {
      // Second region: 'theta' in [phi/2-ALPHA, phi/2]
      minDist *= std::cos(THETAm1)/std::cos(theta);
  }else if(theta < THETA1)
  {
      // Third region: 'theta' in [phi/2, phi/2+ALPHA]
      minDist *= std::cos(THETAm1)/std::cos(theta-phi);
  }else if(theta < THETA2)
  {
      // Forth region: 'theta' in [phi/2+ALPHA, phi/2+HALF_PI]
      minDist *= std::sin(THETA1)/std::sin(theta);
  }else if(theta < THETA3)
  {
      // Fifth region: 'theta' in [phi/2+HALF_PI, phi/2-ALPHA+PI]
      minDist *= std::sin(THETA1)/std::sin(theta-phi);
  }else
  {
      // Sixth region: 'theta' in [phi/2-ALPHA+PI, PI]
      minDist *= std::cos(THETAm1)/(-std::cos(theta));
  }

  // Check second criterion
  return (distance < minDist);
}


