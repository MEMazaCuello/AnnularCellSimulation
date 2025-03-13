#pragma once

#include "GlobalParameters.hpp" // includes <cmath>, <numbers> and <filesystem>
#include "rod.hpp"
#include "grid.hpp"

class AnnularCell {
public:
	[[nodiscard]] auto getRod(const int idx) const -> const Rod&;
	[[nodiscard]] auto getRods() const ->  const std::array<Rod, GP::NUM_RODS>&;

	[[maybe_unused]] auto MCStep() -> double;
	[[maybe_unused]] auto thermalize() -> double;
	[[maybe_unused]] auto MCSimulation() -> double;

	/* Fills the Cell with coordinates saved in file.
		- Each line in filename is exactly of the form: x,y,a
		- Filled only up to GP::NUM_RODS number of rods.
		* FIXME: Does not check for coordinates validity.
	 */
	[[maybe_unused]] auto fillFromFile(const std::filesystem::path& filename) -> bool;
	[[maybe_unused]] auto fill() -> bool;
	[[maybe_unused]] auto save(const std::filesystem::path& filename, const int n) const -> bool;

private:
	[[nodiscard]] inline auto rodIsWithinWalls(const Rod& rod) const -> bool;

	[[nodiscard]] auto isOverlapingNeighbor(const Rod& rod) const -> bool;
	
	[[nodiscard]] inline auto positionIsValid(const Rod& rod) const -> bool;

	auto tryToMoveRod(Rod& rod, [[maybe_unused]] int& count) -> void;
	auto tryToBringRodTowardsCenter(Rod& rod, const double& dr) -> void;

private:
	std::array<Rod, GP::NUM_RODS> m_bundle{};
	Grid m_grid{};
};