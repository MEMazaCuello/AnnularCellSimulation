#include "GlobalParameters.hpp" // includes <cmath> and <numbers>
#include "rod.hpp"

using namespace std::numbers;

auto Rod::moveBy(const double& dx, const double& dy, const double& da) -> void
{
    x += dx;
    y += dy;
    a = std::remainder(a + da, pi); // [-pi/2, pi/2]
}

struct Vec2
{
    double x;
    double y;
};

inline static auto rotateClockwise(const Vec2& v, const Vec2& n) -> Vec2
{
    return {n.x * v.x + n.y * v.y, n.x * v.y - n.y * v.x};
}

inline static auto isInsideRect(const Vec2& v, const double& X_MAX, const double& Y_MAX) -> bool
{
    return (std::abs(v.x) < X_MAX) && (std::abs(v.y) < Y_MAX);
}

[[nodiscard]] auto Rod::overlaps(const Rod& other) const -> bool
{
    const Vec2 P{other.x - x, other.y - y};

    const double sqDist = P.x * P.x + P.y * P.y;
    if (sqDist > GP::ROD::D_SQ)
    {
        return false;
    }
    else if (sqDist < GP::ROD::W_SQ)
    {
        return true;
    }
    else
    {
        const double da = std::abs(std::remainder(other.a - a, pi)); // 0 to pi/2
        const Vec2 aux{ (1.0 + std::cos(da)), std::sin(da) };
        const double X_MAX = GP::ROD::HALF_L * aux.x + GP::ROD::HALF_W * aux.y;
        const double Y_MAX = GP::ROD::HALF_W * aux.x + GP::ROD::HALF_L * aux.y;
    
        const Vec2 n1{ std::cos(a), std::sin(a) };             // Own normal
        const Vec2 n2{ std::cos(other.a), std::sin(other.a) }; // Other normal

        return isInsideRect(rotateClockwise(P, n1), X_MAX, Y_MAX)
            && isInsideRect(rotateClockwise(P, n2), X_MAX, Y_MAX);
    }
}

[[nodiscard]] auto Rod::getCornersRadiiSq() const -> std::array<double, 4>
{
    const double base = x * x + y * y + GP::ROD::HALF_D * GP::ROD::HALF_D;
    const double aplus = a + GP::ROD::ALPHA;
    const double aminus = a - GP::ROD::ALPHA;

    return std::array<double, 4>
    { 
        base + GP::ROD::D * (x * std::cos(aplus)  + y * std::sin(aplus)),
        base - GP::ROD::D * (x * std::cos(aminus) + y * std::sin(aminus)),
        base - GP::ROD::D * (x * std::cos(aplus)  + y * std::sin(aplus)),
        base + GP::ROD::D * (x * std::cos(aminus) + y * std::sin(aminus))
    };
}
