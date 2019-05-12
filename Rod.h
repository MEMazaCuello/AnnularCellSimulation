/****
  'Rod.h':
  --------------------
  Declaration of 'Rod' struct.
  For implementation details, see 'Rod.cpp'.
  
  --------------------
  A 'Rod' is represented by a rectangle with
    width  = 2*'m_HALF_WIDTH' and
    length = 2*'m_HALF_LENGTH'
  at the 2D cartesian position of its center
    (m_xPos, m_yPos)
  with orientation
    m_angle
  defined as the angle between the 'LENGTH'- and the OX-axis.
        
  --------------------
  Last modified: 2019-05-12 
  By: M. E. Maza Cuello
****/

#ifndef ROD_H_INCLUDED
#define ROD_H_INCLUDED

/** 'Rod' DECLARATION **/
struct Rod
{
  /* Static variables */
  // Length variables
  static const double m_HALF_WIDTH;
  static const double m_HALF_LENGTH;
  static const double m_HALF_DIAGONAL;
  
  /* Public variables */
  // Coordinates and orientation
  double m_xPos;    // X coordinate
  double m_yPos;    // Y coordinate
  double m_angle;   // Angle between LENGTH axis and  OX axis

  /* Constructors */
  // Default
  Rod()
      : m_xPos(0.0d), m_yPos(0.0d), m_angle(0.0d) {}

  // Custom
  Rod(const double& xPos, const double& yPos, const double& angle)
      : m_xPos(xPos), m_yPos(yPos), m_angle(angle)
  {}

  /* Methods */
  
  bool isTouchingRod(const Rod& otherRod);
  /** 
    'isTouchingRod': Return 'true' if rod is overlapping 'otherRod', 'false' if not.
                     Assumes 'm_angle' of both rods to be in the interval [-HALF_PI, HALF_PI].
  **/
};

#endif // ROD_H_INCLUDED
