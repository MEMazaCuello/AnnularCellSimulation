/**
  * "Rod.cpp":
  * --------------------
  * Implementation of the Rod structure.
  * For declaration details, see "Rod.h".
  *
  * Needs: "parameters.cpp", "Rod.h".
  *
  * --------------------
  * Last modified: 2020-04-29
  * By: M. E. Maza-Cuello
  */

#include <cmath>

#include "Rod.h"

// Parameters defined in "parameters.cpp" _____________________________
extern const double WIDTH;
extern const double DIAGONAL;
extern const double HALF_WIDTH;
extern const double HALF_LENGTH;
extern const double HALF_DIAGONAL;
extern const double ALPHA;
extern const double HALF_PI;
extern const double PI;

// Static variables ___________________________________________________
const double Rod::m_HALF_WIDTH    = HALF_WIDTH;
const double Rod::m_HALF_LENGTH   = HALF_LENGTH;
const double Rod::m_HALF_DIAGONAL = HALF_DIAGONAL;

// Methods ____________________________________________________________
void Rod::confineAngle()
{
  while (m_angle > HALF_PI)
  {
    m_angle -= PI;
  }

  while (m_angle < -HALF_PI)
  {
    m_angle += PI;
  }
}

void Rod::translate(const double& dx, const double& dy)
{
  m_xPos += dx;
  m_yPos += dy;
}

void Rod::rotate(const double& angle)
{
  m_angle += angle;

  confineAngle();
}

double Rod::sqDistanceTo(const double& x, const double& y)
{
  return (x-m_xPos)*(x-m_xPos)+(y-m_yPos)*(y-m_yPos);
}

double Rod::sqDistanceTo(const Rod& other)
{
  return (other.m_xPos-m_xPos)*(other.m_xPos-m_xPos)+(other.m_yPos-m_yPos)*(other.m_yPos-m_yPos);
}

bool Rod::isWithinRadius(const Rod& other, const double& radius)
{
  return ( radius*radius > sqDistanceTo(other) );
}

double Rod::thresholdDistance(const Rod& other)
{
  double xrel = m_xPos-other.m_xPos;
  double yrel = m_yPos-other.m_yPos;

  // Relative angle in interval [-pi/2,pi/2]
  double phi = m_angle - other.m_angle;
  if (phi > HALF_PI)
  {
    phi -= PI;
  }
  else if (phi < -HALF_PI)
  {
    phi += PI;
  }

  // Relative center-to-center orientation in [-pi,pi]
  double theta = std::atan2(yrel,xrel) - other.m_angle;
  if (theta > PI)
  {
    theta -= 2.0*PI;
  }
  else if (theta < -PI)
  {
    theta += 2.0*PI;
  }

  // Symmetry rotation/reflexion to first or second quadrants
  if (phi < 0.0d)
  {
    phi = -phi;
    theta = (theta < 0.0d) ? -theta : PI - theta;
  }
  else if (theta < 0.0d)
  {
    theta = theta+PI;
  }

  // Analytical computation of threshold distance
  double THETA0  = 0.5*phi;
  double THETAm1 = THETA0-ALPHA;
  double THETA1  = THETA0+ALPHA;
  double THETA2  = THETA0+HALF_PI;
  double THETA3  = THETAm1+PI;

  double thresholdDist = DIAGONAL*std::cos(THETA0);

  if (theta < THETAm1)
  {
    thresholdDist *= std::sin(THETA1)/std::sin(phi-theta);
  }
  else if (theta < THETA0)
  {
    thresholdDist *= std::cos(THETAm1)/std::cos(theta);
  }
  else if (theta < THETA1)
  {
    thresholdDist *= std::cos(THETAm1)/std::cos(theta-phi);
  }
  else if (theta < THETA2)
  {
    thresholdDist *= std::sin(THETA1)/std::sin(theta);
  }
  else if (theta < THETA3)
  {
    thresholdDist *= std::sin(THETA1)/std::sin(theta-phi);
  }
  else
  {
    thresholdDist *= std::cos(THETAm1)/(-std::cos(theta));
  }

  return thresholdDist;
}

bool Rod::isTouchingRod(const Rod& refRod)
{
  double sqDist = sqDistanceTo(refRod);

  if (sqDist < WIDTH*WIDTH)
  {
      return true;
  }
  else if (sqDist > DIAGONAL*DIAGONAL)
  {
      return false;
  }

  // Analytically computed minimum distance between the two rods
  double minDist = thresholdDistance(refRod);

  return ( sqDist < minDist*minDist );
}
