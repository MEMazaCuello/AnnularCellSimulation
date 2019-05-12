/****
  'Rod.h':
  --------------------
  Declaration of 'Rod' struct.
  For implementation details, see 'Rod.cpp'.

  --------------------
  Last modified: 2019-05-12 
  By: M. E. Maza Cuello
****/

#ifndef ROD_H_INCLUDED
#define ROD_H_INCLUDED

struct Rod
{
  /** PUBLIC PARAMETERS **/
  // Intrinsic lengths
  static const double m_HALF_WIDTH;
  static const double m_HALF_LENGTH;
  static const double m_HALF_DIAGONAL;
  // Coordinates and orientation
  double m_xPos;    // X coordinate
  double m_yPos;    // Y coordinate
  double m_angle;   // Angle between LENGTH axis and  OX axis

  /** CONSTRUCTORS **/
  // Default
  Rod()
      : m_xPos(0.0d), m_yPos(0.0d), m_angle(0.0d) {}

  // Custom
  Rod(const double& xPos, const double& yPos, const double& angle)
      : m_xPos(xPos), m_yPos(yPos), m_angle(angle)
  {}

  /** METHODS **/
  // Decide wether this rod is overlapping 'otherRod'
  bool isTouchingRod(const Rod& otherRod);
};

#endif // ROD_H_INCLUDED
