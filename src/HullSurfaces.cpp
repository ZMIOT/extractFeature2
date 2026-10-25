// Ported from pointcloudExtractFeature/algorithms/extract_feature2/core/
// HullSurfaces.cpp. The only host adaptation is the cancellation exception.
#include "HullSurfaces.h"
#include <pcl/common/point_tests.h>
#include <pcl/common/common.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <set>
#include <Eigen/Cholesky>

namespace weld {
namespace {
using Basis = Eigen::Matrix<double,6,1>;
Basis basis(double x,double y) {Basis b;b<<1,x,y,x*x,x*y,y*y;return b;}
float cross(const Eigen::Vector2f& a,const Eigen::Vector2f& b) {return a.x()*b.y()-a.y()*b.x();}
std::vector<Eigen::Vector2f> convexBoundary(std::vector<Eigen::Vector2f> points) {
  std::sort(points.begin(),points.end(),[](const auto& a,const auto& b){return a.x()==b.x()?a.y()<b.y():a.x()<b.x();});
  points.erase(std::unique(points.begin(),points.end(),[](const auto& a,const auto& b){return (a-b).squaredNorm()<1e-12f;}),points.end());
  if(points.size()<3)return {};
  std::vector<Eigen::Vector2f> hull;
  for(const auto& p:points){while(hull.size()>1&&cross(hull.back()-hull[hull.size()-2],p-hull.back())<=0)hull.pop_back();hull.push_back(p);}
  const auto lower=hull.size();
  for(int i=int(points.size())-2;i>=0;--i){const auto& p=points[i];while(hull.size()>lower&&cross(hull.back()-hull[hull.size()-2],p-hull.back())<=0)hull.pop_back();hull.push_back(p);}
  hull.pop_back();return hull;
}
struct SurfacePatch {
  Eigen::Vector3f center,normal,tangent,bitangent;
  Basis coefficients=Basis::Zero();
  std::vector<Eigen::Vector2f> boundary,bodyBoundary;
  int family=1;
  bool contains(const Eigen::Vector2f& uv,float margin)const {
    if(boundary.size()<3)return false;
    for(std::size_t i=0;i<boundary.size();++i){const Eigen::Vector2f edge=boundary[(i+1)%boundary.size()]-boundary[i];if(cross(edge,uv-boundary[i]) < -margin*edge.norm())return false;}
    for(std::size_t i=0;i<bodyBoundary.size();++i){const Eigen::Vector2f edge=bodyBoundary[(i+1)%bodyBoundary.size()]-bodyBoundary[i];if(cross(edge,uv-bodyBoundary[i]) < -margin*edge.norm())return false;}
    return true;
  }
};
// Fit bounded overlapping charts from FIXED body seeds. Recovered points never
// become seeds, so a scan trail cannot extend the model through repeated growth.
std::vector<SurfacePatch> buildPatches(const CloudT::Ptr& anchors,
    const std::vector<int>& anchorIds,const std::vector<Eigen::Vector3f>& normals,
    const std::vector<int>& selected,pcl::KdTreeFLANN<PointT>& tree,float radius,float distance,
    const std::function<void()>& check,const CloudT::Ptr& body={},pcl::KdTreeFLANN<PointT>* bodyTree=nullptr,int widthAxis=1,
    const std::vector<int>* bodySourceIds=nullptr) {
  std::vector<SurfacePatch> patches;std::set<std::array<int,4>> visited;
  std::vector<int> ids;std::vector<float> distances;
  for(std::size_t i=0;i<anchors->size();++i) {
    if((i&255u)==0)check();
    if(selected[anchorIds[i]]<0)continue;
    const auto center=toVector((*anchors)[i]);Eigen::Vector3f normal=normals[anchorIds[i]];
    const float spacing=radius/3;
    const std::array<int,4> key{int(std::floor(center.x()/spacing)),int(std::floor(center.y()/spacing)),int(std::floor(center.z()/spacing)),selected[anchorIds[i]]};
    if(!visited.insert(key).second)continue;
    std::vector<int> bodyIds;std::vector<float> bodyDistances;
    if(bodyTree){bodyTree->radiusSearch((*anchors)[i],radius,bodyIds,bodyDistances);if(bodyIds.size()<24)continue;}
    tree.radiusSearch((*anchors)[i],radius,ids,distances);
    float initialOffset=0;
    if(bodyTree) {
      // The missing strip's point normals are exactly what failed upstream.
      // Estimate a floor-oriented plane from robust raw triples instead of
      // requiring these unreliable normals to agree before fitting anything.
      if(ids.size()<24)continue;
      std::uint32_t state=std::uint32_t(i)+0x9e3779b9u;
      auto pick=[&](){state=1664525u*state+1013904223u;return ids[state%ids.size()];};
      int best=0;const std::size_t stride=std::max<std::size_t>(1,ids.size()/256);
      for(int trial=0;trial<96;++trial) {
        const Eigen::Vector3f a=toVector((*anchors)[pick()]),b=toVector((*anchors)[pick()]),c=toVector((*anchors)[pick()]);
        Eigen::Vector3f candidate=(b-a).cross(c-a);const float length=candidate.norm();if(length<1e-5f)continue;candidate/=length;
        if(std::abs(candidate[widthAxis])>.7f)continue;
        int support=0;for(std::size_t k=0;k<ids.size();k+=stride)support+=std::abs(candidate.dot(toVector((*anchors)[ids[k]])-a))<=distance;
        if(support>best){best=support;normal=candidate;initialOffset=normal.dot(a-center);}
      }
      if(best<12||best*6<int((ids.size()+stride-1)/stride))continue;
    }
    SurfacePatch patch;patch.center=center;patch.normal=normal;patch.tangent=normal.unitOrthogonal();patch.bitangent=normal.cross(patch.tangent);patch.family=selected[anchorIds[i]];
    patch.coefficients[0]=initialOffset;
    std::vector<Eigen::Vector3f> local;
    for(int neighbor:ids) {
      if(!bodyTree&&std::abs(normal.dot(normals[anchorIds[neighbor]]))<.85f)continue;
      const Eigen::Vector3f delta=toVector((*anchors)[neighbor])-center;
      if(std::abs(normal.dot(delta)-initialOffset)>std::max(4*distance,radius*.3f))continue;
      local.emplace_back(delta.dot(patch.tangent)/radius,delta.dot(patch.bitangent)/radius,delta.dot(normal));
    }
    if(local.size()<24)continue;
    // Huber reweighting prevents a few stray seeds from pulling a floor chart.
    for(int iteration=0;iteration<4;++iteration) {
      Eigen::Matrix<double,6,6> matrix=Eigen::Matrix<double,6,6>::Zero();Basis rhs=Basis::Zero();
      for(const auto& p:local) {
        const Basis b=basis(p.x(),p.y());const double residual=std::abs(p.z()-patch.coefficients.dot(b));
        const double weight=iteration==0&&!bodyTree?1.:std::min(1.,double(distance)/std::max(residual,1e-9));
        matrix.noalias()+=weight*b*b.transpose();rhs.noalias()+=weight*b*p.z();
      }
      matrix.diagonal().array()+=1e-7*local.size();patch.coefficients=matrix.ldlt().solve(rhs);
    }
    if(!patch.coefficients.allFinite()||std::abs(patch.coefficients[0])>(bodyTree?radius/3:2*distance)||patch.coefficients.tail<3>().norm()>2*radius)continue;
    std::vector<Eigen::Vector2f> support;
    for(const auto& p:local)if(std::abs(p.z()-patch.coefficients.dot(basis(p.x(),p.y())))<=distance)support.emplace_back(p.x(),p.y());
    if(support.size()<24||double(support.size())<local.size()*(bodyTree?.25:.75))continue;
    patch.boundary=convexBoundary(std::move(support));
    if(patch.boundary.size()<3)continue;
    double area=0;for(std::size_t j=0;j<patch.boundary.size();++j)area+=cross(patch.boundary[j],patch.boundary[(j+1)%patch.boundary.size()]);
    if(area<.12)continue;
    if(bodyTree) {
      // A floor chart fitted from raw samples is allowed only when established
      // body points agree with it on surrounding sides. A nearby parallel
      // interference sheet cannot pass by being inside a 3-D search sphere.
      std::vector<Eigen::Vector2f> contacts;int sameSurfaceContacts=0;
      for(int j:bodyIds) {
        const Eigen::Vector3f delta=toVector((*body)[j])-center;
        const Eigen::Vector2f uv(delta.dot(patch.tangent)/radius,delta.dot(patch.bitangent)/radius);
        const double residual=std::abs(delta.dot(normal)-patch.coefficients.dot(basis(uv.x(),uv.y())));
        const bool sameSurface=residual<=3*distance+radius/18;
        sameSurfaceContacts+=sameSurface;
        // At a floor/wall junction the wall MUST NOT be required to lie on the
        // floor plane. A nearby differently oriented body face can bound the
        // other bank, but cannot replace all same-surface contact evidence.
        const bool joiningFace=bodySourceIds && residual<=radius*.5f &&
          std::abs(normal.dot(normals[(*bodySourceIds)[j]]))<.7f;
        if(sameSurface||joiningFace)contacts.push_back(uv);
      }
      if(sameSurfaceContacts<12||contacts.size()<24)continue;
      patch.bodyBoundary=convexBoundary(std::move(contacts));
      if(patch.bodyBoundary.size()<3||!patch.contains(Eigen::Vector2f::Zero(),0))continue;
    }
    patches.push_back(std::move(patch));
  }
  return patches;
}
}
HullSurfacesResult extractHullSurfaces(const CloudT::ConstPtr &input,
                                     const DetectorConfig &config,
                                     const std::function<bool()> &cancelled) {
  HullSurfacesResult result;
  auto check = [&] { if(cancelled && cancelled()) throw std::runtime_error("Hull extraction cancelled"); };
  const float grid = config.hullSurfaceVoxelSize;
  const int axis=int(config.longitudinalAxis),u=axis==0?1:0;
  if(!input || !std::isfinite(grid) || grid<=0 || axis<0 || axis>2 ||
     !std::isfinite(config.hullSurfaceCurvature) || config.hullSurfaceCurvature<=0 || config.hullSurfaceCurvature>=.34f ||
     !std::isfinite(config.hullSurfaceDistance) || config.hullSurfaceDistance<=0 ||
     !std::isfinite(config.hullPreservationMargin) || config.hullPreservationMargin<0) {
    result.error="Invalid hull surface input or parameters";return result;
  }
  CloudT::Ptr finite(new CloudT);
  for(std::size_t i=0;i<input->size();++i) {
    if((i&4095u)==0)check();
    if(pcl::isFinite((*input)[i]))finite->push_back((*input)[i]);
  }
  if(finite->size()<100) {result.error="Insufficient surface support";return result;}
  PointT minimum,maximum;pcl::getMinMax3D(*finite,minimum,maximum);
  const Eigen::Vector3d span=(toVector(maximum)-toVector(minimum)).cast<double>();
  const double cells=(std::floor(span.x()/grid)+2)*(std::floor(span.y()/grid)+2)*(std::floor(span.z()/grid)+2);
  if(!std::isfinite(cells)||cells>double(std::numeric_limits<std::int32_t>::max())) {
    result.error="Analysis grid is too fine for the input extent";return result;
  }
  CloudT::Ptr samples(new CloudT);
  pcl::VoxelGrid<PointT> voxel;
  voxel.setInputCloud(finite);voxel.setLeafSize(grid,grid,grid);voxel.filter(*samples);
  result.analysisPoints=samples->size();
  pcl::KdTreeFLANN<PointT> tree;tree.setInputCloud(samples);
  const std::size_t n=samples->size();
  std::vector<Eigen::Vector3f> normals(n,Eigen::Vector3f::Zero());
  std::vector<float> curvature(n,1), spread(n,0);
  std::vector<int> ids;std::vector<float> distances;
  const float radius=3.5f*grid;
  for(std::size_t i=0;i<n;++i) {
    if((i&255u)==0)check();
    tree.radiusSearch((*samples)[i],radius,ids,distances,64);
    if(ids.size()<8)continue;
    Eigen::Vector3f center=Eigen::Vector3f::Zero();
    for(int j:ids)center+=toVector((*samples)[j]);center/=float(ids.size());
    Eigen::Matrix3f cov=Eigen::Matrix3f::Zero();
    for(int j:ids){const Eigen::Vector3f d=toVector((*samples)[j])-center;cov+=d*d.transpose();}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(cov);
    if(solver.info()!=Eigen::Success)continue;
    const auto e=solver.eigenvalues();normals[i]=solver.eigenvectors().col(0);
    curvature[i]=e[0]/std::max(1e-10f,e.sum());spread[i]=e[1]/std::max(1e-10f,e[2]);
  }
  std::vector<int> category(n,0);
  for(std::size_t i=0;i<n;++i)category[i]=std::abs(normals[i][u])>.7f?1:0;
  std::vector<int> labels(n,-1);std::vector<std::vector<int>> regions;
  for(std::size_t seed=0;seed<n;++seed) {
    check();
    if(labels[seed]>=0 || curvature[seed]>config.hullSurfaceCurvature || spread[seed]<.08f)continue;
    const int label=int(regions.size());std::vector<int> region{int(seed)};labels[seed]=label;
    for(std::size_t q=0;q<region.size();++q) {
      if((q&255u)==0)check();
      const int current=region[q];const auto &p=(*samples)[current];
      tree.radiusSearch(p,3.5f*grid,ids,distances,64);
      for(int neighbor:ids) {
        if(labels[neighbor]>=0 || curvature[neighbor]>config.hullSurfaceCurvature || spread[neighbor]<.08f)continue;
        // Connectivity belongs to the body, not one normal direction.
        labels[neighbor]=label;region.push_back(neighbor);
      }
    }
    regions.push_back(std::move(region));
  }
  std::vector<int> order(regions.size());std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[&](int a,int b){return regions[a].size()>regions[b].size();});
  for(int id:order) {
    result.regionSizes.push_back(regions[id].size());
    Eigen::Vector3f center=Eigen::Vector3f::Zero(),normal=Eigen::Vector3f::Zero();
    for(int i:regions[id]){center+=toVector((*samples)[i]);normal+=normals[i].cwiseAbs();}
    result.regionCenters.push_back(center/float(regions[id].size()));
    result.regionNormals.push_back(normal/float(regions[id].size()));
  }
  std::vector<int> selected(n,-1);
  // Keep exactly ONE largest spatially continuous surface body. In particular,
  // never fill missing sides by adding the second/third largest components.
  // Rank voxel support rather than raw scan density (repeated scan lines must
  // not let a small, oversampled interference plane win).
  if(order.empty() || regions[order.front()].size()<100){result.error="No supported surface body";return result;}
  const auto& body=regions[order.front()];
  std::vector<float> transverse;transverse.reserve(body.size());
  for(int i:body)transverse.push_back(toVector((*samples)[i])[u]);
  std::nth_element(transverse.begin(),transverse.begin()+transverse.size()/2,transverse.end());
  const float middle=transverse[transverse.size()/2];
  // These labels only colour local normal families; they never select blocks.
  for(int i:body)selected[i]=category[i]==0?1:(toVector((*samples)[i])[u]<middle?0:2);
  // Grow only a short distance from established surface patches. This restores
  // high-curvature intersections and measured boundaries without chaining out
  // along an attached scan trail.
  CloudT::Ptr anchors(new CloudT);std::vector<int> anchorIds;
  for(std::size_t i=0;i<n;++i)if(selected[i]>=0){anchors->push_back((*samples)[i]);anchorIds.push_back(int(i));}
  if(anchors->empty()){result.error="No supported hull surfaces";return result;}
  pcl::KdTreeFLANN<PointT> anchorTree;anchorTree.setInputCloud(anchors);
  // The chart must reach seeds on both banks of a normal-support gap;
  // its convex support boundary, not the search sphere, controls acceptance.
  const float patchRadius=18*grid;
  const float preservationDistance=std::max(config.hullSurfaceDistance,config.hullPreservationMargin);
  auto patches=buildPatches(anchors,anchorIds,normals,selected,anchorTree,patchRadius,config.hullSurfaceDistance,check);
  // Long floor strips may have no accepted seeds at all. Fit candidate charts
  // from measured samples there, but require independent surrounding contact
  // with the already selected main body before allowing any recovery.
  std::vector<int> rawIds(n),floorFamilies(n,-1);std::iota(rawIds.begin(),rawIds.end(),0);
  for(std::size_t i=0;i<n;++i)if(selected[i]<0&&normals[i].squaredNorm()>.5f)floorFamilies[i]=1;
  auto floorPatches=buildPatches(samples,rawIds,normals,floorFamilies,tree,patchRadius,config.hullSurfaceDistance,check,anchors,&anchorTree,u,&anchorIds);
  patches.insert(patches.end(),std::make_move_iterator(floorPatches.begin()),std::make_move_iterator(floorPatches.end()));
  result.referencePatches=patches.size();
  CloudT::Ptr patchCenters(new CloudT);
  for(const auto& patch:patches)patchCenters->push_back(PointT(patch.center.x(),patch.center.y(),patch.center.z()));
  pcl::KdTreeFLANN<PointT> patchTree;if(!patchCenters->empty())patchTree.setInputCloud(patchCenters);
  std::vector<int> patchIds;std::vector<float> patchDistances;
  for(std::size_t i=0;i<input->size();++i) {
    if((i&1023u)==0)check();const auto &p=(*input)[i];if(!pcl::isFinite(p))continue;
    anchorTree.nearestKSearch(p,1,ids,distances);
    if(ids.empty())continue;
    const int j=anchorIds[ids[0]];
    int family=selected[j];
    bool retained=distances[0]<=4.f*grid*grid && std::abs(normals[j].dot(toVector(p)-toVector((*samples)[j])))<=config.hullSurfaceDistance;
    const bool originalSupport=retained;
    bool strictSupport=retained;
    bool nearSurface=false;
    if(!retained&&!patches.empty()&&distances[0]<=patchRadius*patchRadius) {
      // Evaluate all bounded charts in range: adding a nearby wall/raw chart
      // must not crowd a valid floor chart out of a fixed nearest-K shortlist.
      patchTree.radiusSearch(p,patchRadius,patchIds,patchDistances);
      for(std::size_t k=0;k<patchIds.size();++k) {
        if(patchDistances[k]>patchRadius*patchRadius)break;
        const auto& patch=patches[patchIds[k]];const Eigen::Vector3f delta=toVector(p)-patch.center;
        const Eigen::Vector2f uv(delta.dot(patch.tangent)/patchRadius,delta.dot(patch.bitangent)/patchRadius);
        // Only interpolate within measured support, with half a voxel of edge
        // tolerance; never use an unbounded plane to pull in external clutter.
        if(!patch.contains(uv,.5f*grid/patchRadius))continue;
        const double residual=std::abs(delta.dot(patch.normal)-patch.coefficients.dot(basis(uv.x(),uv.y())));
        nearSurface|=residual<3*config.hullSurfaceDistance;
        // A cleanup mask is NOT a fitted surface. Retain original points in
        // this bounded band, including rough weld toes and local deformation.
        // Keep searching for strict support to classify uncertain points
        // independently of patch iteration order. No recovered point grows ROI.
        if(residual<=preservationDistance){retained=true;family=patch.family;}
        if(residual<=config.hullSurfaceDistance){strictSupport=true;break;}
      }
    }
    if(!retained){result.rejectedNearSurfacePoints+=nearSurface;continue;}
    if(!originalSupport)result.recoveredIndices.push_back(i);
    if(!strictSupport)result.protectedUncertainIndices.push_back(i);
    result.cloud->push_back(p);result.sourceIndices.push_back(i);result.surfaceIds.push_back(std::uint32_t(family));
  }
  if(result.cloud->size()<100)result.error="Surface support too small";
  return result;
}
}
