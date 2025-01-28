#include "analysis.hpp"
#include <ranges>
#include <algorithm>
#include <utility>
#include <fstream>
#include <iostream>

auto Analysis::analize(const std::filesystem::path& file_in, const std::filesystem::path& file_out) -> void
{
    cell.fillFromFile(file_in);

    std::ofstream of(file_out);
    if (of.is_open())
    {
        of << std::scientific << std::setprecision(15);

        const auto dirs = computeLocalDirectors();
        const auto params = computeOrderParameters();

        for (int i = 0; i < GP::NUM_RODS; ++i)
        {
            const Rod& rod = cell.getRod(i);

            of << rod.x << "," << rod.y << "," << rod.a << "," 
               << dirs[i] << "," << params[i].q2 << "," << params[i].q4 << "," << params[i].qS << '\n';
        }
        of.close();
    }
}

[[nodiscard]] auto Analysis::getRegions() const -> std::vector<std::forward_list<int>>
{
    std::vector<std::forward_list<int>> regions{ GP::NUM_RODS };

    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        const Rod& ref = cell.getRod(i);
        for (int j = 0; j <= i; ++j)
        {
            const Rod& rod = cell.getRod(j);
            if ((rod.x - ref.x) * (rod.x - ref.x) + (rod.y - ref.y) * (rod.y - ref.y) < GP::ANALYSIS::R_SQ)
            {
                regions[i].push_front(j);
                regions[j].push_front(i);
            }
        }
    }

    return regions;
}

[[nodiscard]] auto Analysis::computeLocalDirectors() const -> std::vector<double>
{
    std::vector<double> directors{};
    directors.reserve(GP::NUM_RODS);

    const std::vector<std::forward_list<int>> regions = getRegions();
    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        double cos2a = 0.0;
        double sin2a = 0.0;
        for (const int j : regions[i])
        {
            const double angle = 2.0 * cell.getRod(j).a;
            // Note that indexes are mixed
            cos2a += std::cos(angle);
            sin2a += std::sin(angle);
        }
        directors[i] = std::remainder( std::atan2(std::sqrt(cos2a * cos2a + sin2a * sin2a) - cos2a, sin2a), std::numbers::pi);
    }

    return directors;
}

[[nodiscard]] auto Analysis::computeQ2() const -> std::vector<double>
{
    const auto dir = computeLocalDirectors();
    const auto regions = getRegions();

    std::vector<double> q2{ GP::NUM_RODS };
    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        int size = 0;
        for (const int j : regions[i])
        {
            q2[i] += std::cos(2.0 * (cell.getRod(j).a - dir[i]));
            ++size;
        }
        q2[i] /= size;
    }

    return q2;
}

[[nodiscard]] auto Analysis::computeQ4() const -> std::vector<double>
{
    const auto dir = computeLocalDirectors();
    const auto regions = getRegions();

    std::vector<double> q4{ GP::NUM_RODS };
    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        int size = 0;
        for (const int j : regions[i])
        {
            q4[i] += std::cos(4.0 * (cell.getRod(j).a - dir[i]));
        }
        q4[i] /= size;
    }

    return q4;
}

[[nodiscard]] auto Analysis::computeQS() const -> std::vector<double>
{
    const auto dir = computeLocalDirectors();
    const auto regions = getRegions();

    constexpr double K = 2.0 * std::numbers::pi * GP::ANALYSIS::Q;
    std::vector<double> qS{ GP::NUM_RODS };
    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        int size{ 0 };
        double cs{ 0.0 };
        double sn{ 0.0 };
        for (const int j : regions[i])
        {
            const double angle = K * (std::cos(dir[i]) * (cell.getRod(j).x - cell.getRod(i).x) 
                                    + std::sin(dir[i]) * (cell.getRod(j).y - cell.getRod(i).y));
            cs += std::cos(angle);
            sn += std::sin(angle);
            ++size;
        }
        qS[i] = std::sqrt(cs*cs + sn*sn) / size;
    }

    return qS;
}

[[nodiscard]] auto Analysis::computeOrderParameters() const -> std::vector<Params>
{
    const auto dir = computeLocalDirectors();
    const auto regions = getRegions();
    constexpr double K = 2.0 * std::numbers::pi * GP::ANALYSIS::Q;

    std::vector<Params> params{ GP::NUM_RODS };
    for (int i = 0; i < GP::NUM_RODS; ++i)
    {
        int size{ 0 };
        double cs{ 0.0 };
        double sn{ 0.0 };        
        for (const int j : regions[i])
        {
            params[i].q2 += std::cos(2.0 * (cell.getRod(j).a - dir[i]));
            params[i].q4 += std::cos(4.0 * (cell.getRod(j).a - dir[i]));

            const double angle = K * (std::cos(dir[i]) * (cell.getRod(j).x - cell.getRod(i).x)
                + std::sin(dir[i]) * (cell.getRod(j).y - cell.getRod(i).y));
            cs += std::cos(angle);
            sn += std::sin(angle);

            ++size;
        }
        params[i].q2 /= size;
        params[i].q4 /= size;
        params[i].qS = std::sqrt(cs * cs + sn * sn) / size;
    }

    return params;
}
