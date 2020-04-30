/**
  * "Rod.h":
  * --------------------
  * Declaration of the Rod structure.
  * For implementation details, see "Rod.cpp".
  *
  * --------------------
  * A "Rod" is represented by a rectangle with
  *  width  = 2*"m_HALF_WIDTH" and
  *  length = 2*"m_HALF_LENGTH"
  * at the 2D cartesian position of its center
  *  (m_xPos, m_yPos)
  * and with long-axis orientation
  *   m_angle
  * defined as the angle between the "LENGTH"- and the OX-axis,
  * in radians, set in the interval [-pi/2,pi/2].
  *
  * --------------------
  * Last modified: 2020-04-29
  * By: M. E. Maza-Cuello
  */

#ifndef ROD_H_INCLUDED
#define ROD_H_INCLUDED

struct Rod
{
  static const double m_HALF_WIDTH;
  static const double m_HALF_LENGTH;
  static const double m_HALF_DIAGONAL;

  double m_xPos;
  double m_yPos;
  double m_angle;

  Rod()
    : m_xPos(0.0d), m_yPos(0.0d), m_angle(0.0d) {}

  Rod(const double& xPos, const double& yPos, const double& angle)
    : m_xPos(xPos), m_yPos(yPos), m_angle(angle)
  {
    confineAngle();
  }

  /**
    * Sets m_angle in the interval [ -pi/2 , pi/2 ]
    */
  void confineAngle();

  void translate(const double& dx, const double& dy);
  void rotate(const double& angle);

  /**
    * Returns squared euclidean distance to (x,y) or (other.m_xPos,other.m_yPos)
    */
  double sqDistanceTo(const double& x, const double& y);
  double sqDistanceTo(const Rod& other);

  bool isWithinRadius(const Rod& other, const double& radius);

  /**
    * Returns analytical threshold distance between this and other rod.
    * The threshold distance is the minimum distance between the centers
    * of two rectangles of the same size, given their relative position
    * and angle, such that the rods do not overlap.
    */
  double thresholdDistance(const Rod& other);

  bool isTouchingRod(const Rod& other);
};

#endif // ROD_H_INCLUDED
