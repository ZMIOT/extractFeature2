#include "GpuSurfaceNormals.h"

#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value,const char *message){if(!value)throw std::runtime_error(message);}
void environment(const char *name,const char *value){
#ifdef _WIN32
  _putenv_s(name,value);
#else
  if(*value)setenv(name,value,1);else unsetenv(name);
#endif
}
}

int main() try {
  environment("WELD_SEAM_GPU","1");
  weld::CloudT source;
  std::vector<std::vector<std::size_t>> groups(4);
  for(int group=0;group<4;++group)for(int i=0;i<180;++i) {
    const float x=.02f*i;
    const float y=.18f*group+.15f*x+.003f*std::sin(.13f*i+group);
    groups[group].push_back(source.size());source.push_back(weld::PointT(0,x,y));
  }
  weld::detail::GpuSurfaceNormals gpu;
  if(!weld::detail::computeSurfaceNormalsGpu(source,groups,1,2,.75f,gpu)) {
    std::cout<<"OpenCL GPU unavailable; surface-normal equivalence skipped.\n";
    return 0;
  }
  require(gpu.normals.size()==groups.size(),"GPU group count changed");
  for(std::size_t g=0;g<groups.size();++g) {
    weld::CloudT::Ptr profile(new weld::CloudT);
    for(auto index:groups[g])profile->push_back(weld::PointT(source[index].y,source[index].z,0));
    pcl::KdTreeFLANN<weld::PointT> tree;tree.setInputCloud(profile);
    std::vector<int> ids;std::vector<float> distances;
    for(std::size_t i=0;i<profile->size();++i) {
      tree.radiusSearch((*profile)[i],.75f,ids,distances,64);
      const bool supported=ids.size()>=5;
      require((gpu.supported[g][i]!=0)==supported,"GPU neighbourhood support differs from PCL");
      if(!supported)continue;
      Eigen::Vector2f center=Eigen::Vector2f::Zero();
      for(int id:ids)center+=Eigen::Vector2f((*profile)[id].x,(*profile)[id].y);
      center/=float(ids.size());
      Eigen::Matrix2f covariance=Eigen::Matrix2f::Zero();
      for(int id:ids){const Eigen::Vector2f d=Eigen::Vector2f((*profile)[id].x,(*profile)[id].y)-center;covariance.noalias()+=d*d.transpose();}
      covariance/=float(ids.size());
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> fit(covariance);
      const bool valid=fit.info()==Eigen::Success&&fit.eigenvalues()[1]>=1.e-5f&&
          fit.eigenvalues()[0]<=.15f*fit.eigenvalues().sum();
      if(!valid){require(gpu.normals[g][i].squaredNorm()<1.e-12f,"GPU accepted invalid covariance");continue;}
      require(std::abs(gpu.normals[g][i].dot(fit.eigenvectors().col(0)))>.9999f,
              "GPU normal differs from PCL covariance normal");
    }
  }
  std::cout<<"OpenCL surface normals match PCL on "<<gpu.device<<".\n";
  return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
