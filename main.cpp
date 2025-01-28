#include <iostream>
#include <format>
#include <chrono>
#include "annularCell.hpp"
#include "analysis.hpp"
#include <fstream>
#include <random>

int main()
{
    /* Create a new MC simulation __________________________________________ */
    using std::chrono::steady_clock;

    AnnularCell cell{};
    if (cell.fill())
    {
        steady_clock::time_point tic{ steady_clock::now() };
        double mean_acceptance = cell.thermalize();
        steady_clock::time_point toc{ steady_clock::now() };
        
        std::cout << std::format("Thermalization duration: {} s\n", 0.001 * std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count());
        std::cout << std::format("Mean acceptance: {}%\n", mean_acceptance);

        cell.save(GP::IO::THERMALIZED, GP::NUM_RODS);

        for (int iter = 0; iter < GP::MC::MC_ITERATIONS; ++iter)
        {
            tic = steady_clock::now();
            double mean_acceptance = cell.MCSimulation();
            toc = steady_clock::now();

            std::cout << std::format(" --- ITERATION {} OF {} --- \n",1 + iter, GP::MC::MC_ITERATIONS);
            std::cout << std::format("Duration: {} s\n", std::chrono::duration_cast<std::chrono::seconds>(toc - tic).count());
            std::cout << std::format("Mean acceptance: {}%\n\n", mean_acceptance);

            std::filesystem::path filename = GP::IO::MC_BASE;
            (filename += std::to_string(iter)) += GP::IO::MC_EXT;
            cell.save(filename, GP::NUM_RODS);
        }
    }


    /* Continue MC simulation on saved configuration _______________________ */
    // using std::chrono::steady_clock;
    //
    // AnnularCell cell{};
    // if (cell.fillFromFile(GP::IO::INITIAL))
    // {
    //     steady_clock::time_point tic{ steady_clock::now() };
    //     double mean_acceptance = cell.thermalize();
    //     steady_clock::time_point toc{ steady_clock::now() };
    //        
    //     std::cout << std::format("Thermalization duration: {} s\n", 0.001 * std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count());
    //     std::cout << std::format("Mean acceptance: {}%\n", mean_acceptance);
    //
    //     cell.save(GP::IO::THERMALIZED, GP::NUM_RODS);
    //
    //     for (int iter = 0; iter < GP::MC::MC_ITERATIONS; ++iter)
    //     {
    //         tic = steady_clock::now();
    //         double mean_acceptance = cell.MCSimulation();
    //         toc = steady_clock::now();
    //
    //         std::cout << std::format(" --- ITERATION {} OF {} --- \n",1 + iter, GP::MC::MC_ITERATIONS);
    //         std::cout << std::format("Duration: {} s\n", std::chrono::duration_cast<std::chrono::seconds>(toc - tic).count());
    //         std::cout << std::format("Mean acceptance: {}%\n\n", mean_acceptance);
    //
    //         std::filesystem::path filename = GP::IO::MC_BASE;
    //         (filename += std::to_string(iter)) += GP::IO::MC_EXT;
    //         cell.save(filename, GP::NUM_RODS);
    //     }
    // }


    /* Run analysis on saved configuration _________________________________ */
    //const std::filesystem::path FILE_IN{ "initial_configuration.csv" };
    //const std::filesystem::path FILE_OUT{ "analyzed_configuration.csv" };
    //
    //Analysis analysis{};
    //analysis.analize(FILE_IN, FILE_OUT);
}
