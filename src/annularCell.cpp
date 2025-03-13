#include "annularCell.hpp" // includes <cmath> and <numbers>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <ranges>
#include <algorithm>

using std::numbers::pi;

std::mt19937 gen(std::random_device{}());
std::uniform_real_distribution<double> rand_dw(-GP::MC::dW, GP::MC::dW);
std::uniform_real_distribution<double> rand_dl(-GP::MC::dL, GP::MC::dL);
std::uniform_real_distribution<double> rand_da(-GP::MC::dA, GP::MC::dA);

struct Vec2
{
    double x;
    double y;
};

inline static auto coordinatesSFD(const Vec2& v, const Vec2& n) -> Vec2
{
    return {std::max(0.0, std::abs(n.x*v.x + n.y*v.y) - GP::ROD::HALF_L), 
            std::max(0.0, std::abs(n.x*v.y - n.y*v.x) - GP::ROD::HALF_W)};
}

[[nodiscard]] inline auto AnnularCell::rodIsWithinWalls(const Rod& rod) const -> bool
{
    const double sqDist = rod.x * rod.x + rod.y * rod.y;

    if (sqDist > GP::CHECKS::MAX_IN_DIST_SQ && sqDist < GP::CHECKS::MIN_OUT_DIST_SQ) [[likely]]
    {   // Clearly inside
        return true;
    }
    else if (sqDist < GP::CHECKS::MIN_IN_DIST_SQ 
          || sqDist > GP::CHECKS::MAX_OUT_DIST_SQ
          || std::ranges::any_of(rod.getCornersRadiiSq(),
                                [](const double& r_sq) { return (r_sq < GP::CELL::R_IN_SQ)
                                                             || (r_sq > GP::CELL::R_OUT_SQ); })) [[unlikely]]
    {   // Center too close to boundary || Any corner outside
        return false;
    }
    else [[unlikely]]
    {   // Chek R_IN from SDF for rounded rectangle
        const Vec2 P{ rod.x, rod.y };
        const Vec2 n{ std::cos(rod.a), std::sin(rod.a)};
        const Vec2 d = coordinatesSFD(P, n);
        return ((d.x*d.x + d.y*d.y) > GP::CELL::R_IN_SQ);
    }
}

[[nodiscard]] auto AnnularCell::isOverlapingNeighbor(const Rod& rod) const -> bool
{
    const auto isOverlaping = [&](int n) { return (n!=rod.index) && rod.overlaps(m_bundle[n]); };
    for (const auto& neighborBoxIndex : m_grid.m_neighborBoxesIndexes[m_grid.getBoxIndexAt(rod.x, rod.y)])
    {
        if (std::ranges::any_of(m_grid.m_boxes[neighborBoxIndex], isOverlaping))
        {
            return true;
        }
    }
    return false;
}

[[nodiscard]] inline auto AnnularCell::positionIsValid(const Rod& rod) const -> bool
{
    return rodIsWithinWalls(rod) && !isOverlapingNeighbor(rod);
}

auto AnnularCell::tryToMoveRod(Rod& rod, [[maybe_unused]] int& count) -> void
{
    const double dx = rand_dl(gen);
    const double dy = rand_dw(gen);

    Rod newRod(rod);
    newRod.moveBy(dx * std::cos(newRod.a) - dy * std::sin(newRod.a),
                  dx * std::sin(newRod.a) + dy * std::cos(newRod.a),
                  rand_da(gen));

    if (positionIsValid(newRod))
    {
        m_grid.moveIndex(rod.index, rod.x, rod.y, newRod.x, newRod.y);
        rod = newRod;
        ++count;
    }
}

auto AnnularCell::tryToBringRodTowardsCenter(Rod& rod, const double& dr) -> void
{
    const double theta = std::atan2(rod.y, rod.x);
    
    Rod newRod = rod;
    newRod.moveBy(-dr * std::cos(theta), -dr * std::sin(theta), 0.0);

    if (positionIsValid(newRod))
    {
        m_grid.moveIndex(rod.index, rod.x, rod.y, newRod.x, newRod.y);
        rod = newRod;
    }
}

[[nodiscard]] auto AnnularCell::getRod(const int idx) const -> const Rod&
{
    return m_bundle[idx];
}

[[nodiscard]] auto AnnularCell::getRods() const -> const std::array<Rod, GP::NUM_RODS>&
{
    return m_bundle;
}

[[maybe_unused]] auto AnnularCell::MCStep() -> double
{
    int successes{ 0 };
    std::ranges::for_each(m_bundle, [&](Rod& rod){ tryToMoveRod(rod, successes); });
    return (100.0 * successes) / GP::NUM_RODS;
}

[[maybe_unused]] auto AnnularCell::thermalize() -> double
{
    double mean_acceptance{ 0.0 };
    for (int s = 0; s < GP::MC::THERMAL_STEPS; s++)
    {
        mean_acceptance += MCStep();
    }
    return mean_acceptance / GP::MC::THERMAL_STEPS;
}

[[maybe_unused]] auto AnnularCell::MCSimulation() -> double
{
    double mean_acceptance{ 0.0 };
    for (int s = 0; s < GP::MC::MC_STEPS; s++)
    {
        mean_acceptance += MCStep();
    }
    return mean_acceptance / GP::MC::MC_STEPS;
}

[[maybe_unused]] auto AnnularCell::fillFromFile(const std::filesystem::path& filename) -> bool
{
    std::ifstream infile(filename);
    if (infile.is_open())
    {
        char c; // Only for commas
        int i = 0;
        for (std::string line; std::getline(infile, line) && i < GP::NUM_RODS; ++i)
        {
            std::istringstream ss(line);
            ss >> m_bundle[i].x >> c >> m_bundle[i].y >> c >> m_bundle[i].a;
            m_bundle[i].index = i;
            m_grid.addIndexAt(i, m_bundle[i].x, m_bundle[i].y);
        }
        infile.close();
        return true;
    }
    else
    {
        std::cout << "FILE " << filename <<" COULD NOT BE OPENED!\n";
        return false;
    }
}

[[maybe_unused]] auto AnnularCell::fill() -> bool
{
    int current_index = 0;
    Rod rod{};
    { // Fill in rings, angle is "tangent" to outer radius
        double r_max = GP::CHECKS::R_OUT_MAX;
        double offset = 0.0;
        while (r_max > GP::CHECKS::R_OUT_MIN && current_index < GP::NUM_RODS)
        {
            const double beta = 2.0 * std::atan(GP::ROD::HALF_L / (r_max - GP::ROD::HALF_W));
            const double spaces = std::floor((2.0 * pi) / beta);
            const double buffer = (2.0 * pi - spaces * beta) / spaces;
            double theta = 0.0;
            do 
            {
                rod.x = r_max * std::cos(theta + offset);
                rod.y = r_max * std::sin(theta + offset);
                rod.a = remainder(theta + offset - 0.5 * pi, pi);
                rod.index = current_index;

                if (positionIsValid(rod))
                {
                    m_bundle[current_index] = rod;
                    m_grid.addIndexAt(current_index, rod.x, rod.y);
                    ++current_index;
                }

                theta += beta + buffer;
            } while (theta < 2 * pi && current_index < GP::NUM_RODS);
            offset += 0.5 * beta;
            r_max -= GP::ROD::HALF_W;
            r_max = 0.9993 * std::sqrt(r_max * r_max + GP::ROD::HALF_D * GP::ROD::HALF_D + r_max * GP::ROD::D * std::cos(0.5 * beta + std::asin(2.0 * r_max * std::sin(0.5 * beta) / GP::ROD::D)));
        }
    }
    
    if (current_index == GP::NUM_RODS)
    {
        std::cout << "Ring filling ended succesfully.\n";

        // Try to attract rods towards the center, so that the ones in the outer rim
        // have the possibility to move a bit
        for (int n = 0; n < 20; ++n)
        {
            std::ranges::for_each(m_bundle | std::views::reverse, [&](Rod& rod) { tryToBringRodTowardsCenter(rod, (GP::ROD::W / (n + 1))); });
        }

        if (save("default_initial_Configuration.csv", GP::NUM_RODS))
        {
            std::cout << "Configuration saved as 'default_initial_Configuration.csv'.\n";
            return true;
        }
        else
        {
            std::cout << "WARNING: COULD NOT SAVE CONFIGURATION!\n";
            return false;
        }
    }
    else
    {
        std::cout << "Ring filling ended. " << GP::NUM_RODS - current_index << " rods missing.\n";
        if (save("default_initial_Configuration_INCOMPLETE.csv", current_index))
        {
            std::cout << "Current configuration saved as 'default_initial_Configuration_INCOMPLETE.csv'.\n";
        }
        else
        {
            std::cout << "WARNING: COULD NOT SAVE CONFIGURATION!\n";
        }
        
        std::cout << "Current cover fraction: " << current_index * GP::ROD::W * GP::ROD::L / (pi * (GP::CELL::R_OUT_SQ - GP::CELL::R_IN_SQ)) << '\n';
        std::cout << "Target cover fraction: " << GP::NUM_RODS * GP::ROD::W * GP::ROD::L / (pi * (GP::CELL::R_OUT_SQ - GP::CELL::R_IN_SQ)) << '\n';
        std::cout << "Starting addition in random locations. This process is not guaranteed to finish.\n";
        
        std::uniform_real_distribution<double> distR(-GP::CELL::R_OUT, GP::CELL::R_OUT);
        current_index -= 1;
        for (; current_index < GP::NUM_RODS; ++current_index)
        {
            m_bundle[current_index].index = current_index;
            do 
            {
                const double x = distR(gen);
                const double y = distR(gen);
                m_bundle[current_index].moveBy(x, y, std::atan2(y, x));
            } while (!positionIsValid(m_bundle[current_index]));

            m_grid.addIndexAt(current_index, m_bundle[current_index].x, m_bundle[current_index].y);
            std::cout << "Rod #" << current_index << " has been inserted." << std::endl;
            ++current_index;
        }
        std::cout << "Random addition finished.\n";
        if (save("default_initial_Configuration.csv", GP::NUM_RODS))
        {
            std::cout << "Configuration saved as 'default_initial_Configuration.csv'. \n";
            return true;
        }
        else
        {
            std::cout << "WARNING: COULD NOT SAVE CONFIGURATION!\n";
            return false;
        }
    }
}

[[maybe_unused]] auto AnnularCell::save(const std::filesystem::path& filename, const int n) const -> bool
{
    std::ofstream of(filename);

    if (of.is_open())
    {
        auto print = [&](const Rod& rod) {of << rod.x << "," << rod.y << "," << rod.a << '\n';};
        
        of << std::scientific << std::setprecision(15);
        std::ranges::for_each(m_bundle | std::views::take(n), print);
        of.close();

        return true;
    }
    else
    {
        std::cout << "FILE " << filename << " COULD NOT BE OPENED!\n";
        return false;
    }
}