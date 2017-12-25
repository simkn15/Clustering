#ifndef TRANSCLUST_COMMON_HPP
#define TRANSCLUST_COMMON_HPP
#include <string>
#include <cmath>
#include <vector>

// TransClust Commons
namespace TCC
{
	/* Config struct with fluint interface */
	struct TransClustParams {
	  	 std::string   file_type              =   "SIMPLE"; // LEGACY
	  	 float        sim_fallback           =   0.0;
	  	 float        threshold                 =   0.0;
	  	 float        p                      =   1.0;
	  	 float        f_att                  =   100.0;
	  	 float        f_rep                  =   100.0;
	  	 unsigned      R                      =   100;
	  	 unsigned      dim                    =   3;
	  	 float        start_t                =   100;
	  	 float        d_init                 =   0.01;
	  	 float        d_maximal              =   5.0;
	  	 float        s_init                 =   0.01;
	  	 float        f_s                    =   0.01;
	  	 float        fpt_time_limit         =   5;
	  	 float        fpt_max_cost           =   5000;
	  	 float        fpt_step_size          =   10;
	  	 bool          disable_force          =   false;
	  	 bool          disable_fpt            =   false;
	  	 unsigned      seed                   =   42;
		TransClustParams& set_file_type(std::string val){file_type = val;return *this;}
		TransClustParams& set_sim_fallback(float val){sim_fallback = val;return *this;}
		TransClustParams& set_threshold(float val){threshold = val;return *this;}
		TransClustParams& set_p(float val){p = val;return *this;}
		TransClustParams& set_f_att(float val){f_att = val;return *this;}
		TransClustParams& set_f_rep(float val){f_rep = val;return *this;}
		TransClustParams& set_R(float val){R = val;return *this;}
		TransClustParams& set_dim(float val){dim = val;return *this;}
		TransClustParams& set_start_t(float val){start_t = val;return *this;}
		TransClustParams& set_d_init(float val){d_init = val;return *this;}
		TransClustParams& set_d_maximal(float val){d_maximal = val;return *this;}
		TransClustParams& set_s_init(float val){s_init = val;return *this;}
		TransClustParams& set_f_s(float val){f_s = val;return *this;}
		TransClustParams& set_fpt_time_limit(float val){fpt_time_limit = val;return *this;}
		TransClustParams& set_fpt_max_cost(float val){fpt_max_cost = val;return *this;}
		TransClustParams& set_fpt_step_size(float val){fpt_step_size = val;return *this;}
		TransClustParams& set_disable_force(bool val){disable_force = val;return *this;}
		TransClustParams& set_disable_fpt(bool val){disable_fpt = val;return *this;}
		TransClustParams& set_seed(bool val){seed = val;return *this;}
	};

	inline float round(float d){
		return std::rint(d*100000)/100000;

	}

	inline float dist(std::vector<std::vector<float>>& pos,unsigned i, unsigned j)
	{
		float res = 0;
		for(unsigned d = 0; d < pos[0].size(); d++){
			float side = pos[i][d] - pos[j][d];
			res += side*side;
		}
		return std::sqrt(res);
	}
}
#endif
