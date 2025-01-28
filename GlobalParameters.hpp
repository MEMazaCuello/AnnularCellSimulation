#pragma once

#include <cmath>
#include <numbers>
#include <filesystem>

/** 
 * TUNABLE PARAMETERS 
 */
namespace GP 
{
	inline constexpr unsigned int NUM_RODS{ 2983 };

	namespace ROD
	{
		inline constexpr double W{ 1.0 };
		inline constexpr double L{ 4.0 };

		// Size requirements
		static_assert(W > 0.0);
		static_assert(W <= L);
	}

	namespace CELL
	{
		inline constexpr double R_OUT{ 70.0 };
		inline constexpr double R_IN{ 20.0 };

		// Size requirements
		static_assert(R_IN > 0.0);
		static_assert(R_IN < R_OUT);
		// Space requirements
		static_assert(R_OUT * R_OUT > (R_IN + ROD::W) * (R_IN + ROD::W) + 0.25 * ROD::L * ROD::L);
		static_assert(NUM_RODS * ROD::W * ROD::L < std::numbers::pi * (R_OUT * R_OUT - R_IN * R_IN));
	}

	namespace GRID
	{
		inline constexpr int BOXES_PER_SIDE{ 35 };

		// Parity requirement
		static_assert(BOXES_PER_SIDE % 2 == 1);
	}

	namespace MC
	{
		inline constexpr double dW{ 0.01 * ROD::W };
		inline constexpr double dL{ 0.01 * ROD::L };
		inline constexpr double dA{ 0.01 * std::numbers::pi };
		inline constexpr int THERMAL_STEPS{ 1'000'000 };
		inline constexpr int MC_STEPS{ 10'000 };
		inline constexpr int MC_ITERATIONS{ 24 }; // Number of repetitions of MC_STEPS
	}

	namespace IO
	{
		inline const std::filesystem::path INITIAL{ "intial_configuration.csv" };
		inline const std::filesystem::path THERMALIZED{ "thermalized_configuration.csv" };
		inline const std::filesystem::path MC_BASE{ "configuration_" };
		inline const std::filesystem::path MC_EXT{ ".csv" };
	}

	namespace ANALYSIS
	{
		/*
		  Radius of the circular region around a rod for computing local order parameters
		*/
		inline constexpr double R_SQ{ 16.0 * ROD::L * ROD::L };

		/*
		  Expected inverse distance between smectic layers
		*/
		inline constexpr double Q{ 1.0 / (1.01 * ROD::L) };
	}
}

/**
 * ---------------------------------------------------------------------------
 * ___________________________ DERIVED CONSTANTS _____________________________ 
 * _____________________________ DO NOT CHANGE  ______________________________
 * ---------------------------------------------------------------------------
 */
namespace GP
{
	namespace ROD 
	{
		inline constexpr double W_SQ{ W * W };
		inline constexpr double L_SQ{ L * L };
		inline constexpr double HALF_W{ 0.5 * W };
		inline constexpr double HALF_L{ 0.5 * L };
		inline const double ALPHA = std::atan2(W, L); // Constexpr in C++26
		inline const double D = std::sqrt( W_SQ + L_SQ ); // Constexpr in C++26
		inline const double D_SQ = D * D; // Constexpr in C++26
		inline const double HALF_D = 0.5 * D; // Constexpr in C++26
	}

	namespace CELL 
	{
		inline constexpr double R_OUT_SQ{ R_OUT * R_OUT };
		inline constexpr double R_IN_SQ{ R_IN * R_IN };
	}

	namespace GRID
	{
		inline constexpr int NUM_BOXES{ BOXES_PER_SIDE * BOXES_PER_SIDE };
		inline constexpr int CENTRAL_BOX{ (BOXES_PER_SIDE - 1) / 2 };
		inline constexpr int CENTRAL_INDEX{ (NUM_BOXES - 1) / 2};
		inline constexpr double BOX_W{ 2.0 * CELL::R_OUT / (BOXES_PER_SIDE - 2) }; // There is an empty, 1-box-wide frame
		inline constexpr double BOX_INV_W{ 1.0 / BOX_W };
	}

	namespace CHECKS
	{
		inline constexpr double R_IN_PLUS_HALF_L{ CELL::R_IN + ROD::HALF_L };
		inline constexpr double R_IN_PLUS_HALF_W{ CELL::R_IN + ROD::HALF_W };
		inline const     double HALF_D_OVER_R_IN{ GP::ROD::HALF_D / GP::CELL::R_IN }; // Constexpr in C++26
		
		inline const	 double D_OVER_D_OUT{ GP::ROD::HALF_D / GP::CELL::R_OUT }; // Constexpr in C++26
		inline const	 double D_OVER_D_OUT_SQ{ D_OVER_D_OUT * D_OVER_D_OUT }; // Constexpr in C++26
		inline const	 double ONE_MINUS_D_OVER_D_OUT_SQ{ 1.0 - D_OVER_D_OUT_SQ }; // Constexpr in C++26
		
		inline const	 double R_OUT_MAX{ std::sqrt(CELL::R_OUT_SQ - ROD::HALF_L * ROD::HALF_L) - GP::ROD::HALF_W }; // Constexpr in C++26
		inline const	 double R_OUT_MIN{ std::sqrt((CELL::R_IN + ROD::W) * (CELL::R_IN + ROD::W) + ROD::HALF_L * ROD::HALF_L) };
		
		inline const     double MAX_OUT_DIST_SQ{ R_OUT_MAX * R_OUT_MAX };
		inline const     double MIN_OUT_DIST_SQ{ (CELL::R_OUT - ROD::HALF_D) * (CELL::R_OUT - ROD::HALF_D) }; // Constexpr in C++26
		inline constexpr double MIN_IN_DIST_SQ{ (CELL::R_IN + ROD::HALF_W) * (CELL::R_IN + ROD::HALF_W) };
		inline const     double MAX_IN_DIST_SQ{ (CELL::R_IN + ROD::HALF_D) * (CELL::R_IN + ROD::HALF_D) }; // Constexpr in C++26
		
		inline const	 double MIN_AUX_ANGLE{ std::atan2(GP::ROD::HALF_W, GP::CHECKS::R_IN_PLUS_HALF_L) }; // Constexpr in C++26
		inline const	 double MAX_AUX_ANGLE{ std::atan2(GP::CHECKS::R_IN_PLUS_HALF_W, GP::ROD::HALF_L) }; // Constexpr in C++26
	}
}