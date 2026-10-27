#ifndef SINK_BOUND_CHECK_HPP_
#define SINK_BOUND_CHECK_HPP_
#ifdef AMREX_USE_GPU
#error "This isolated Sink binding experiment is CPU only."
#endif
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "hydro/hydro_system.hpp"
#include "fundamental_constants.H"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>

namespace quokka::SinkBoundCheck {
inline bool configured=false,require_bound=false,audit=false;
inline amrex::Long step_candidates=0,step_bound=0,step_invalid=0;
inline FILE *candidate_file=nullptr;
inline void configure() {
 if(configured)return;
 configured=true;amrex::ParmParse pp("particles");
 pp.query("sink_creation_require_bound",require_bound);
 pp.query("sink_creation_bound_audit",audit);
}
struct Energy {
 double mass=0.,thermal=0.,thermal_conservative=0.,kinetic=0.,gravity=0.;
 int cells=0;bool finite=true;
 bool bound()const{return finite && thermal+kinetic+gravity<0.;}
};
struct Cell {double mass,vx,vy,vz,x,y,z;};
template<class problem_t>
inline Energy evaluate(amrex::Array4<const amrex::Real> const &s,int i,int j,int k,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const &dx) {
 static_assert(!Physics_Traits<problem_t>::is_mhd_enabled,"Gas-only binding test excludes MHD.");
 Energy e;
 // The 3-cell control sphere is fully inside the existing four ghost cells.
 // Cubic cells are required for the cube self-energy coefficient below.
 AMREX_ALWAYS_ASSERT(std::abs(dx[1]/dx[0]-1.)<1.e-12 && std::abs(dx[2]/dx[0]-1.)<1.e-12);
 std::array<Cell,123> c{};double pvx=0.,pvy=0.,pvz=0.;const double vol=dx[0]*dx[1]*dx[2];
 for(int a=-3;a<=3;++a)for(int b=-3;b<=3;++b)for(int d=-3;d<=3;++d) {
  if(a*a+b*b+d*d>9)continue;
  const int ii=i+a,jj=j+b,kk=k+d;
  double rho=s(ii,jj,kk,HydroSystem<problem_t>::density_index);
  double aux=s(ii,jj,kk,HydroSystem<problem_t>::internalEnergy_index);
  double et=s(ii,jj,kk,HydroSystem<problem_t>::energy_index);
  double px=s(ii,jj,kk,HydroSystem<problem_t>::x1Momentum_index);
  double py=s(ii,jj,kk,HydroSystem<problem_t>::x2Momentum_index);
  double pz=s(ii,jj,kk,HydroSystem<problem_t>::x3Momentum_index);
  if(!std::isfinite(rho)||rho<=0.||!std::isfinite(aux)||aux<=0.||!std::isfinite(et)||
     !std::isfinite(px)||!std::isfinite(py)||!std::isfinite(pz)){e.finite=false;return e;}
  double m=rho*vol,vx=px/rho,vy=py/rho,vz=pz/rho;
  c[e.cells++]={m,vx,vy,vz,a*dx[0],b*dx[1],d*dx[2]};
  e.mass+=m;pvx+=m*vx;pvy+=m*vy;pvz+=m*vz;
  // Use the stored physical thermal energy, not unresolved bulk kinetic energy.
  e.thermal+=aux*vol;
  e.thermal_conservative+=(et-.5*(px*px+py*py+pz*pz)/rho)*vol;
 }
 AMREX_ALWAYS_ASSERT(e.cells==123);
 const double vx=pvx/e.mass,vy=pvy/e.mass,vz=pvz/e.mass;
 for(int a=0;a<e.cells;++a) {
  const auto &q=c[a];double x=q.vx-vx,y=q.vy-vy,z=q.vz-vz;
  e.kinetic+=.5*q.mass*(x*x+y*y+z*z);
  // Uniform cubic-cell self energy; pair interactions use cell centres.
  e.gravity-=.941156*C::Gconst*q.mass*q.mass/dx[0];
  for(int b=0;b<a;++b) {
   double xx=q.x-c[b].x,yy=q.y-c[b].y,zz=q.z-c[b].z;
   e.gravity-=C::Gconst*q.mass*c[b].mass/std::sqrt(xx*xx+yy*yy+zz*zz);
  }
 }
 e.finite=std::isfinite(e.mass)&&std::isfinite(e.thermal)&&std::isfinite(e.thermal_conservative)&&
          std::isfinite(e.kinetic)&&std::isfinite(e.gravity);
 return e;
}
inline void record(double time,int i,int j,int k,double dx,const Energy &e) {
 if(!candidate_file) {
  char name[128];std::snprintf(name,sizeof(name),"binding_candidates_rank%04d.tsv",amrex::ParallelDescriptor::MyProc());
  candidate_file=std::fopen(name,"a");AMREX_ALWAYS_ASSERT(candidate_file!=nullptr);
  std::fprintf(candidate_file,"# time i j k dx cells mass Eth_aux Eth_cons Ekin W finite bound\n");
 }
 ++step_candidates;step_bound+=e.bound();step_invalid+=!e.finite;
 std::fprintf(candidate_file,"%.17g %d %d %d %.17g %d %.17g %.17g %.17g %.17g %.17g %d %d\n",
    time,i,j,k,dx,e.cells,e.mass,e.thermal,e.thermal_conservative,e.kinetic,e.gravity,int(e.finite),int(e.bound()));
}
inline void flushStep(double time) {
 configure();
 if(!require_bound&&!audit)return;
 amrex::Long counts[3]={step_candidates,step_bound,step_invalid};
 amrex::ParallelDescriptor::ReduceLongSum(counts,3);
 if(counts[0])amrex::Print()<<"SINK_BINDING time="<<time<<" eligible="<<counts[0]<<" bound="<<counts[1]
    <<" unbound="<<(counts[0]-counts[1]-counts[2])<<" invalid="<<counts[2]<<" enabled="<<require_bound<<"\n";
 if(candidate_file)std::fflush(candidate_file);
 step_candidates=step_bound=step_invalid=0;
}
} // namespace quokka::SinkBoundCheck
#endif
