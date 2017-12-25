#ifndef FINDCONNECTEDCOMPONENTS_HPP
#define FINDCONNECTEDCOMPONENTS_HPP
#include <vector>
#include <deque>
#include <queue>
#include "transclust/Common.hpp"
#include "transclust/ConnectedComponent.hpp"

namespace FCC{

	/****************************************************************************
	 * FIND CCs IN CC W. THRESHOLD
	 ***************************************************************************/
	void findConnectedComponents(
			TCC::TransClustParams& tcp,
			const ConnectedComponent &cc,
			std::deque<ConnectedComponent> &ccs,
			const float threshold);

	/* Breadth-first search */
	void BFS_cc(
			std::vector<std::vector<unsigned>>& membership,
			const ConnectedComponent &cc,
			const float threshold);
}
#endif
