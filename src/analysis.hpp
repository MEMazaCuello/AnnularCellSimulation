#pragma once

#include "annularCell.hpp"

struct Params {
	double q2;
	double q4;
	double qS;
};

struct Analysis {
	AnnularCell cell;

	auto analize(const std::filesystem::path& file_in, const std::filesystem::path& file_out) -> void;

	[[nodiscard]] auto getRegions() const->std::vector<std::forward_list<int>>;

	[[nodiscard]] auto computeLocalDirectors() const -> std::vector<double>;

	[[nodiscard]] auto computeQ2() const -> std::vector<double>;
	[[nodiscard]] auto computeQ4() const -> std::vector<double>;
	[[nodiscard]] auto computeQS() const -> std::vector<double>;

	[[nodiscard]] auto computeOrderParameters() const -> std::vector<Params>;
};
