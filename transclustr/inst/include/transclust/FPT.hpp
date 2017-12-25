#ifndef FPT_HPP
#define FPT_HPP
#include <iomanip>
#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <limits>
#include <chrono>
#include "transclust/Common.hpp"
#include "transclust/ConnectedComponent.hpp"
#include "transclust/ClusteringResult.hpp"


class FPT{

 	public:
		struct Node {
			std::vector<std::vector<unsigned>> nodeParents;
			std::vector<std::vector<float>> edgeCost;
			unsigned size;
			float cost;
		};

		FPT(
			ConnectedComponent& cc,
			float time_limit,
			float stepSize,
			float mockInf
		);

		void cluster(ClusteringResult &cr);

	private:

		ConnectedComponent &cc;
		std::chrono::time_point<std::chrono::system_clock> start; 
		float time_limit;
		float stepSize;
		float maxK;
		float inf;
		int LEVEL = 0;

		bool solution_found;
		float solution_cost;
		std::vector<std::vector<unsigned>> solution_nodeParents;
		std::vector<std::vector<float>> solution_edgeCost;


		inline float getDeltaTime()
		{
			std::chrono::duration<float> diff = (std::chrono::system_clock::now()-start);
			return diff.count();
		}
		void reduce(Node& fptn);
		float costSetForbidden(
				Node& fptn, 
				unsigned node_i,
				unsigned node_j);

		float costSetPermanent(
				Node& fptn, 
				unsigned node_i,
				unsigned node_j);
		void find_solution(Node& fptn0);
		void mergeNodes(
			Node& fptn, 
			unsigned i,
			unsigned j, 
			float costForMerging);
		
		void buildSolution(ClusteringResult &cr);
		void clone_node(Node& fptn0,Node& fptn1);
};
#endif
