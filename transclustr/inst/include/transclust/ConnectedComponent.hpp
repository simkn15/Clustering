#ifndef CONNECTEDCOMPONENT_HPP
#define CONNECTEDCOMPONENT_HPP
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "transclust/Common.hpp"
#include "transclust/TriangularMatrix.hpp"

static unsigned _cc_id(0);

class ConnectedComponent
{
	public:
		//ConnectedComponent(const std::string &filename,bool use_custom_fallback,float sim_fallback,std::string ft);
		ConnectedComponent(const std::string &filename,TCC::TransClustParams& tcp);
		//ConnectedComponent(std::vector<float>& sim_matrix_1d,unsigned num_o,bool use_custom_fallback,float sim_fallback);
		ConnectedComponent(std::vector<float>& sim_matrix_1d,unsigned num_o,TCC::TransClustParams& tcp);//,bool use_custom_fallback,float sim_fallback);
		ConnectedComponent(const ConnectedComponent& cc,const std::vector<unsigned>& objects, float th,TCC::TransClustParams& tcp);
		/*
			ConnectedComponent(const std::vector<std::vector<float>>& pos,float th);
			*/
		inline const TriangularMatrix& getMatrix() const { return m; }
		inline const unsigned size()const { return m.getNumObjects(); }
		inline const float getMinSimilarity() const{ return m.getMinValue(); }
		inline const float getMaxSimilarity() const{ return m.getMaxValue(); }
		inline const std::vector<std::string>& getObjectNames() const
		{
			return m.getObjectNames();
		}
		inline const std::string getObjectName(unsigned i)const { return m.getObjectName(i);}

		// get cost for some threshold
		inline const float getCost(unsigned i, unsigned j, float t) const { return TCC::round(m(i, j) - t);}

		// get the value of the edge in the cc 
		inline const float at(unsigned i, unsigned j, bool normalized = true) const
		{
			float sim = m(i, j) - threshold;
			if (normalized)
			{

				if( sim > 0){
					return TCC::round(sim / normalization_context_positive);
				}else{
					return TCC::round(sim / normalization_context_negative);
				}
			}
			else
			{
				return TCC::round(sim);
			}
		}
		inline const float getThreshold() const
		{
			return threshold;
		};

		inline const unsigned getObjectId(unsigned i) const
		{
			return m.getObjectId(i);
		};
		inline const unsigned getId() const
		{
			return id;
		};
		inline const unsigned getMatrixSize() const
		{
			return m.getMatrixSize();
		};

		void dump();

	private:
		inline void init_normalization_context(TCC::TransClustParams& tcp)
		{
			float max_val = std::max(
					TCC::round(std::abs(m.getMaxValue()-threshold)),
					TCC::round(std::abs(m.getMinValue()-threshold))
					);
			normalization_context_positive = max_val;
			normalization_context_negative = max_val;
		
		}
		unsigned id;
		TriangularMatrix m;
		float threshold;
		float normalization_context_positive;
		float normalization_context_negative;
		float cost;

		inline static unsigned getNewId()
		{
			_cc_id++;
			return _cc_id;
		}
};
#endif
