#pragma once

#include "GlobalParameters.hpp" // includes <cmath> and <numbers>
#include <array>
#include <forward_list>

consteval std::array<std::array<int, 9>, GP::GRID::NUM_BOXES> setNeighborBoxes() 
{
    // Evaluates at compile time the neighboring index of a box
    // Note that the boundaries are *not* treated differently
    std::array<std::array<int, 9>, GP::GRID::NUM_BOXES> neighborBoxes{};
    for (int i = 0; i < GP::GRID::NUM_BOXES; ++i)
    {
        // The order is in *inverse* probability (central box last)
        // because the function getNeighborhoodOf() returns a push_front of the forward_lists
        neighborBoxes[i] = { i, i - 1, i + 1,
                             i - GP::GRID::BOXES_PER_SIDE, i + GP::GRID::BOXES_PER_SIDE,
                             i - 1 + GP::GRID::BOXES_PER_SIDE, i + 1 + GP::GRID::BOXES_PER_SIDE,
                             i - 1 - GP::GRID::BOXES_PER_SIDE, i + 1 - GP::GRID::BOXES_PER_SIDE };

    }
    return neighborBoxes;
};

class Grid
{
public:

    [[nodiscard]] auto getBoxIndexAt(const double& x, const double& y) const -> int;
    
    auto addIndexAt(const int idx, const double& x, const double& y) -> void;
    auto moveIndex(const int idx, const double& from_x, const double& from_y, const double& to_x, const double& to_y) -> void;

//private:

    std::array<std::forward_list<int>, GP::GRID::NUM_BOXES> m_boxes{};
    static constexpr std::array<std::array<int, 9>, GP::GRID::NUM_BOXES> m_neighborBoxesIndexes{ setNeighborBoxes() };

};

