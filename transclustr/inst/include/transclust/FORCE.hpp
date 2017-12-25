#ifndef FORCE_HPP
#define FORCE_HPP
#include <vector>
#include <iostream>
#include "transclust/Common.hpp"
#include "transclust/ConnectedComponent.hpp"
#include "transclust/ClusteringResult.hpp"

namespace FORCE{
	float dist(std::vector<std::vector<float>>& pos,unsigned i, unsigned j);

	void layout(
			const ConnectedComponent& cc,
			std::vector<std::vector<float>>& pos,
			float p,
			float f_att,
			float f_rep,
			const unsigned R,
			float start_t,
			const unsigned dim,
			unsigned seed = 42);

	void partition(
			const ConnectedComponent& cc,
			std::vector<std::vector<float>>& pos,
			ClusteringResult& cs,
			float d_init,
			float d_maximal,
			float s_init,
			float f_s);

	std::vector<std::vector<unsigned>> geometricLinking(
			std::vector<std::vector<float>>& pos,
			const float maxDist,
			const std::vector<std::vector<unsigned>>& objects);
}
#endif
