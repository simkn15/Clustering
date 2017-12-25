#ifdef _OPENMP
#	include <omp.h>
#endif
#include "transclust/TransClust.hpp"
#include "transclust/ConnectedComponent.hpp"
#include "transclust/FindConnectedComponents.hpp"
#include "transclust/FORCE.hpp"
#include "transclust/FPT.hpp"
#include "transclust/ClusteringResult.hpp"
#include "transclust/Result.hpp"

TransClust::TransClust(
   const std::string& filename,
   TCC::TransClustParams& _tcp
)
{
   tcp = _tcp;
   // Read input similarity file
   ConnectedComponent sim_matrix(filename,tcp);
   id2object = sim_matrix.getObjectNames();
   FCC::findConnectedComponents(tcp,sim_matrix,ccs,tcp.threshold);
}

TransClust::TransClust(
   std::vector<float>& sim_matrix_1d,
   unsigned num_o,
   TCC::TransClustParams& _tcp
)
{
   tcp = _tcp;
   ConnectedComponent sim_matrix(sim_matrix_1d,num_o,tcp);
   id2object = sim_matrix.getObjectNames();
   FCC::findConnectedComponents(tcp,sim_matrix,ccs,tcp.threshold);
}


clustering TransClust::cluster()
{

   Result result(id2object);
   #pragma omp parallel
   {
      #pragma omp for schedule(dynamic)
      for(unsigned i = 0; i < ccs.size(); i++)
      {
         ConnectedComponent& cc = ccs.at(i);

         ClusteringResult cr;

         // set initial cost to negativ, indicating 'no solution found (yet)'
         cr.cost = -1;

         // if cc is at least a conflict tripple
         if(cc.size() > 2){

            /*******************************************************************
             * Cluster using FORCE
             ******************************************************************/
            if(!tcp.disable_force){
               // init position array
               std::vector<std::vector<float>> pos;
               pos.resize(cc.size(), std::vector<float>(tcp.dim,0));

               // layout
               FORCE::layout(cc, pos, tcp.p, tcp.f_att, tcp.f_rep, tcp.R, tcp.start_t, tcp.dim);
               // partition
               FORCE::partition(cc, pos, cr, tcp.d_init, tcp.d_maximal, tcp.s_init, tcp.f_s);
            }

            /*******************************************************************
             * Cluster using FPT
             ******************************************************************/
            if(cr.cost <= tcp.fpt_max_cost && !tcp.disable_fpt){

               float tmp_force_cost = cr.cost;
               FPT fpt(cc,tcp.fpt_time_limit,tcp.fpt_step_size,cr.cost+1);
               fpt.cluster(cr);

               if(cr.cost < 0){
                  cr.cost = tmp_force_cost;
               }
            }
         }else{
            // cc consist of 1 or 2 nodes and is a cluster
            cr.cost = 0;
            cr.membership = std::vector<unsigned>(cc.size(),0);
         }
         #pragma omp critical
         {
            result.add(cc,cr);
         }
      }
   }
   return result.get();
}

//void TransClust::init(TCC::TransClustParams& tcp){
//   // set class vars from tcp
//   // general vars
//   use_custom_fallback = tcp.use_custom_fallback;
//   sim_fallback = tcp.sim_fallback;
//   use_custom_range = tcp.use_default_interval;
//   threshold_min = tcp.th_min;
//   threshold_max = tcp.th_max;
//   threshold_step = tcp.th_step;
//
//   // Layout values
//   p = tcp.p;
//   f_att = tcp.f_att;
//   f_rep = tcp.f_rep;
//   R = tcp.R;
//   start_t = tcp.start_t;
//   dim = tcp.dim;
//
//   // partitioning values
//   d_init = tcp.d_init;
//   d_maximal = tcp.d_maximal;
//   s_init = tcp.s_init;
//   f_s = tcp.f_s;
//
//   // FPT values
//   fpt_time_limit = tcp.fpt_time_limit;
//   fpt_max_cost = tcp.fpt_max_cost;
//   fpt_step_size = tcp.fpt_step_size;
//
//   disable_force = tcp.disable_force;
//   disable_fpt = tcp.disable_fpt;
//
//   seed = tcp.seed;
//
//}
