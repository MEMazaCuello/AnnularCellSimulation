#include "grid.hpp"

#include <iostream>
#include <algorithm>

[[nodiscard]] auto Grid::getBoxIndexAt(const double& x, const double& y) const -> int
{
	return GP::GRID::CENTRAL_INDEX + static_cast<int>(std::round(x * GP::GRID::BOX_INV_W) - std::round(y * GP::GRID::BOX_INV_W) * GP::GRID::BOXES_PER_SIDE);
}

auto Grid::addIndexAt(const int idx, const double& x, const double& y) -> void
{
	m_boxes[getBoxIndexAt(x, y)].emplace_front(idx);
}

auto Grid::moveIndex(const int idx, const double& from_x, const double& from_y, const double& to_x, const double& to_y) -> void
{
	if (from_x != to_x || from_y != to_y)
	{
		// Removes first coincidence of idx
		auto& from = m_boxes[getBoxIndexAt(from_x, from_y)];
		auto oit = from.before_begin();
		auto it = std::next(oit);
		while (it != from.end()) {
			if ((*it) == idx)
			{
				auto& to = m_boxes[getBoxIndexAt(to_x, to_y)];
				to.splice_after(to.cbefore_begin(), from, oit);
				break;
			}
			oit = it++;
		}
	}
}
