#ifndef ROD_H_INCLUDED
#define ROD_H_INCLUDED

struct Rod
{
    double m_xPos;
    double m_yPos;
    double m_angle;
    static const double m_HALF_WIDTH;
    static const double m_HALF_LENGTH;
    static const double m_HALF_DIAGONAL;

    Rod()
        : m_xPos(0.0d), m_yPos(0.0d), m_angle(0.0d) {}

    Rod(const double& xPos, const double& yPos, const double& angle)
        : m_xPos(xPos), m_yPos(yPos), m_angle(angle)
    {}

    bool isTouchingRod(const Rod& otherRod);
};

#endif // ROD_H_INCLUDED
