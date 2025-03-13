#pragma once
#include <array>

struct Rod
{
	double x;
	double y;
	double a; // angle in interval [-HALF_PI, HALF_PI]
	int index;

	auto moveBy(const double& dx, const double& dy, const double& da) -> void;
	
	[[nodiscard]] auto overlaps(const Rod& other) const -> bool;
	[[nodiscard]] auto getCornersRadiiSq() const -> std::array<double, 4>;
};
