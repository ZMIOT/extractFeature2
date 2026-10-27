#include "SurfaceClouds.h"
#include "SurfaceFitting.h"
#include "GpuSurfaceNormals.h"
#include <pcl/common/point_tests.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <map>
#include <array>
#include <chrono>
#include <limits>

namespace weld {
namespace {
float median(std::vector<float> values) {
  auto middle = values.begin() + values.size()/2;
  std::nth_element(values.begin(), middle, values.end());
  return *middle;
}

struct SidePlane {
  Eigen::Vector2f normal=Eigen::Vector2f::Zero();
  float offset=0, begin=0, end=0;
  bool valid=false;
  bool toePatch=false;
};

// Fit the whole supported wall in a longitudinal section. This transverse
// line defines a local side plane containing the section's longitudinal axis.
// Independent sections follow a bent workpiece without a single unbounded
// global plane. RANSAC rejects clutter; orthogonal refitting uses all inliers.
SidePlane fitSidePlane(const std::vector<Eigen::Vector2f> &source,
                       const Eigen::Vector2f &toe, float tolerance, float radius,
                       float minimumSlope, float maximumToeDistance) {
  SidePlane result;
  if(source.size()<6) return result;
  std::vector<Eigen::Vector2f> samples;
  const auto stride=std::max<std::size_t>(1,(source.size()+511)/512);
  for(std::size_t i=0;i<source.size();i+=stride)samples.push_back(source[i]);
  std::size_t best=0;
  float bestToe=std::numeric_limits<float>::max();
  std::uint32_t state=0x6a09e667u;
  auto pick=[&](){state=1664525u*state+1013904223u;return samples[(state>>8)%samples.size()];};
  for(int trial=0;trial<128;++trial) {
    const auto first=pick();
    const Eigen::Vector2f delta=pick()-first;
    if(delta.norm()<radius)continue;
    Eigen::Vector2f normal(-delta.y(),delta.x());normal.normalize();
    if(normal.x()<0)normal=-normal;
    const float offset=normal.dot(first), distance=std::abs(normal.dot(toe)-offset);
    if(normal.x()<minimumSlope || distance>maximumToeDistance)continue;
    std::size_t count=0;
    for(const auto &p:samples)count+=std::abs(normal.dot(p)-offset)<=tolerance;
    if(count>best || (count==best && distance<bestToe)) {
      best=count;bestToe=distance;result.normal=normal;result.offset=offset;
    }
  }
  if(best<6 || best<samples.size()/5)return SidePlane();
  for(int iteration=0;iteration<3;++iteration) {
    std::vector<Eigen::Vector2f> inliers;
    Eigen::Vector2f center=Eigen::Vector2f::Zero();
    for(const auto &p:source)if(std::abs(result.normal.dot(p)-result.offset)<=tolerance) {
      inliers.push_back(p);center+=p;
    }
    if(inliers.size()<6)return SidePlane();
    center/=float(inliers.size());
    Eigen::Matrix2f covariance=Eigen::Matrix2f::Zero();
    for(const auto &p:inliers) {const Eigen::Vector2f d=p-center;covariance.noalias()+=d*d.transpose();}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> solver(covariance);
    if(solver.info()!=Eigen::Success)return SidePlane();
    Eigen::Vector2f normal=solver.eigenvectors().col(0);
    if(normal.x()<0)normal=-normal;
    const float offset=normal.dot(center);
    if(normal.x()<minimumSlope || std::abs(normal.dot(toe)-offset)>maximumToeDistance)return SidePlane();
    result.normal=normal;result.offset=offset;
    const Eigen::Vector2f tangent(-normal.y(),normal.x());
    result.begin=std::numeric_limits<float>::max();result.end=-result.begin;
    for(const auto &p:inliers) {const float t=tangent.dot(p);result.begin=std::min(result.begin,t);result.end=std::max(result.end,t);}
  }
  result.valid=result.end-result.begin>2*radius;
  return result;
}

// An emerging groove can have an accepted lower turn before its wall is long
// or steep enough for the broad-wall fit. Fit only the measured transition next
// to that turn. The bounded patch must rise away from the measured floor chord;
// a flat continuation of the floor, a remote shoulder or a detached sheet is not
// evidence for a side. No neighboring frame supplies output points or offsets.
SidePlane fitToePatch(const CloudT &profile, const std::vector<bool> &supported,
                      const Eigen::Vector2f &toe, const Eigen::Vector2f &outward,
                      const Eigen::Vector2f &up, float radius, float tolerance) {
  const float reach=2*radius;
  const float fitTolerance=std::min(tolerance,.25f*radius);
  std::vector<Eigen::Vector2f> samples;
  for(std::size_t j=0;j<profile.size();++j) {
    if(!supported[j])continue;
    const Eigen::Vector2f p(profile[j].x,profile[j].y),delta=p-toe;
    const float along=outward.dot(delta),height=up.dot(delta);
    if(along>=0 && along<=reach && height>=-tolerance && height<=reach)
      samples.push_back(p);
  }
  auto plane=fitSidePlane(samples,toe,fitTolerance,.2f*radius,0,tolerance);
  if(!plane.valid)return SidePlane();
  Eigen::Vector2f tangent(-plane.normal.y(),plane.normal.x());
  if(tangent.dot(outward)<0)tangent=-tangent;
  // Require a resolved rise relative to the local floor (at least 3 degrees),
  // not a global height/slope test on a tilted scan.
  if(tangent.dot(up)<.052336f || tangent.dot(outward)<=0)return SidePlane();
  plane.toePatch=true;
  return plane;
}
}

void extractSurfaceClouds(const CloudT::ConstPtr &source,
                          const std::vector<SurfaceSection> &sections,
                          const DetectorConfig &config, SeamResult &result,
                          const std::vector<ScanFrameRange> *frames,
                          bool generateSurfaceModels) {
  result.leftSurface.reset(new CloudT);
  result.rightSurface.reset(new CloudT);
  result.weldSurface.reset(new CloudT);
  result.leftSurfaceModel=BSplineSurfaceModel{};
  result.weldSurfaceModel=BSplineSurfaceModel{};
  result.rightSurfaceModel=BSplineSurfaceModel{};
  if (!source || sections.empty()) return;
  const auto classificationStarted=std::chrono::steady_clock::now();
  const int a = int(config.longitudinalAxis), u = a == 0 ? 1 : 0,
            v = a == 2 ? 1 : 2;
  const float radius = std::max(.5f, config.profileBinSize);
  const float tolerance = std::max(config.hullSurfaceDistance, config.profileWallTolerance);
  // Retry a short patch with a baseline wider than measurement uncertainty
  // when fine profile binning leaves too little single-frame fit support.
  const float shortRadius = std::max(radius,2*tolerance);
  const float binSize = std::max(.1f, .5f*config.profileBinSize);
  std::vector<std::vector<std::size_t>> groups(sections.size());
  if (frames && !frames->empty()) {
    // Use identities, not physical longitudinal coordinates: tilted/reversed
    // scans may overlap in XYZ while belonging to different sections.
    std::map<int, std::size_t> byFrame;
    for (std::size_t i=0; i<sections.size(); ++i) byFrame.emplace(sections[i].frameIndex, i);
    for (const auto &frame : *frames) {
      const auto found = byFrame.find(frame.frameIndex);
      if (found == byFrame.end()) continue;
      for (std::size_t j=frame.begin; j<frame.begin+frame.count; ++j)
        if (pcl::isFinite((*source)[j])) groups[found->second].push_back(j);
    }
  } else {
    for (std::size_t j=0; j<source->size(); ++j) {
      const auto &p = (*source)[j];
      if (!pcl::isFinite(p)) continue;
      const float station = p.data[a];
      auto next = std::lower_bound(sections.begin(), sections.end(), station,
          [](const SurfaceSection &s, float x) { return s.station < x; });
      std::size_t i = next == sections.end() ? sections.size()-1 : next-sections.begin();
      if (i && std::abs(sections[i-1].station-station) < std::abs(sections[i].station-station)) --i;
      // Half a station on either side, even at a gap/end. Each raw point is
      // assigned once; overlapping detector slices cannot duplicate output.
      if (std::abs(sections[i].station-station) <= .5f*config.sliceStep + 1.e-4f)
        groups[i].push_back(j);
    }
  }
  std::vector<unsigned char> labels(source->size(), 0);
  detail::GpuSurfaceNormals gpuNormalBatch;
  const bool gpuNormals=detail::computeSurfaceNormalsGpu(
      *source,groups,u,v,radius,gpuNormalBatch);
  if(gpuNormals)
    result.stageMilliseconds["gpu_surface_normals"]=gpuNormalBatch.milliseconds;
  for (std::size_t s=0; s<sections.size(); ++s) {
    const auto &section = sections[s];
    const auto &indices = groups[s];
    if (indices.size()<6 || section.left.data[u]>=section.right.data[u]) continue;
    const Eigen::Vector2f left(section.left.data[u], section.left.data[v]);
    const Eigen::Vector2f right(section.right.data[u], section.right.data[v]);
    const float middle = .5f*(left.x()+right.x());
    CloudT::Ptr profile(new CloudT);
    for (auto i : indices) {
      const auto &p = (*source)[i];
      profile->push_back(PointT(p.data[u], p.data[v], 0));
    }
    std::vector<Eigen::Vector2f> normals(indices.size(),Eigen::Vector2f::Zero());
    std::vector<bool> supported(indices.size(),false);
    std::array<std::vector<Eigen::Vector2f>,2> wallSeeds;
    if(gpuNormals) {
      for(std::size_t j=0;j<indices.size();++j) {
        supported[j]=gpuNormalBatch.supported[s][j]!=0;
        normals[j]=gpuNormalBatch.normals[s][j];
      }
    } else {
      pcl::KdTreeFLANN<PointT> tree;
      tree.setInputCloud(profile);
      std::vector<int> neighbors;
      std::vector<float> distances;
      for (std::size_t j=0; j<indices.size(); ++j) {
        const auto &p = (*profile)[j];
        tree.radiusSearch(p, radius, neighbors, distances, 64);
        if (neighbors.size()<5) continue;
        supported[j]=true;
        Eigen::Vector2f center = Eigen::Vector2f::Zero();
        for (int k : neighbors) center += Eigen::Vector2f((*profile)[k].x, (*profile)[k].y);
        center /= float(neighbors.size());
        Eigen::Matrix2f covariance = Eigen::Matrix2f::Zero();
        for (int k : neighbors) {
          const Eigen::Vector2f delta = Eigen::Vector2f((*profile)[k].x, (*profile)[k].y)-center;
          covariance.noalias() += delta*delta.transpose();
        }
        covariance /= float(neighbors.size());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> fit(covariance);
        if (fit.info()!=Eigen::Success || fit.eigenvalues()[1]<1.e-5f ||
            fit.eigenvalues()[0]>.15f*fit.eigenvalues().sum()) continue;
        normals[j]=fit.eigenvectors().col(0);
      }
    }
    for (std::size_t j=0; j<indices.size(); ++j) {
      if(!supported[j])continue;
      const auto &p = (*profile)[j];
      const Eigen::Vector2f normal=normals[j];
      const Eigen::Vector2f point(p.x, p.y);
      const bool isLeft = p.x < middle;
      const auto &toe = isLeft ? left : right;
      Eigen::Vector2f wall = isLeft ? section.leftNormal : section.rightNormal;
      const bool knownWall = wall.squaredNorm()>.5f;
      if (knownWall) wall.normalize();
      // Broad wall support, including points far above the toe; local normal
      // gates exclude horizontal outer shoulders before fitting the big face.
      const float minimumSlope=std::max(std::sin(config.minSurfaceAngleDegrees*.01745329252f),
                                        knownWall ? .75f*std::abs(wall.x()) : .7f);
      const bool wallDirection = (!knownWall || std::abs(normal.dot(wall))>.5f) &&
                                 std::abs(normal.x())>minimumSlope;
      Eigen::Vector2f tangent(-wall.y(), wall.x());
      const int opening = config.interiorOpeningDirection;
      if ((opening && tangent.y()*opening<0) || (!opening && tangent.x()*(isLeft?-1:1)<0)) tangent = -tangent;
      const bool outward=isLeft ? p.x<=toe.x()+tolerance : p.x>=toe.x()-tolerance;
      const bool aboveToe=opening ? (p.y-toe.y())*opening>tolerance :
                                    knownWall && tangent.dot(point-toe)>tolerance;
      if (wallDirection && outward && aboveToe) {
        wallSeeds[isLeft?0:1].push_back(point);
      }
    }
    const float minimumSlope=std::sin(config.minSurfaceAngleDegrees*.01745329252f);
    const float maximumToeDistance=std::max(3*radius,config.profileCurvatureRadius);
    std::array<SidePlane,2> planes{{
      fitSidePlane(wallSeeds[0],left,tolerance,radius,minimumSlope,maximumToeDistance),
      fitSidePlane(wallSeeds[1],right,tolerance,radius,minimumSlope,maximumToeDistance)}};
    const Eigen::Vector2f floorDirection=(right-left).normalized();
    Eigen::Vector2f up(-floorDirection.y(),floorDirection.x());
    if(config.interiorOpeningDirection<0)up=-up;
    // Without a specified opening, use the detector's wall direction when it
    // is available. Do not guess the missing side of an ambiguous flat scan.
    bool knownOpening=config.interiorOpeningDirection!=0;
    if(!knownOpening) {
      for(int side=0;side<2 && !knownOpening;++side) {
        const auto &normal=side==0?section.leftNormal:section.rightNormal;
        if(normal.squaredNorm()<.5f)continue;
        Eigen::Vector2f tangent(-normal.y(),normal.x());
        if(tangent.dot(floorDirection)*(side==0?-1:1)<0)tangent=-tangent;
        if(std::abs(tangent.dot(up))<.052336f)continue;
        if(tangent.dot(up)<0)up=-up;
        knownOpening=true;
      }
    }
    if(knownOpening)for(int side=0;side<2;++side)if(!planes[side].valid) {
      planes[side]=fitToePatch(*profile,supported,side==0?left:right,
          side==0?Eigen::Vector2f(-floorDirection):floorDirection,up,radius,tolerance);
      if(!planes[side].valid && shortRadius>radius)
        planes[side]=fitToePatch(*profile,supported,side==0?left:right,
            side==0?Eigen::Vector2f(-floorDirection):floorDirection,up,shortRadius,tolerance);
    }
    const bool shortSection=planes[0].toePatch || planes[1].toePatch;
    auto inside=[&](const Eigen::Vector2f &p) {
      if(!planes[0].valid || !planes[1].valid)return false;
      if(shortSection) {
        const float along=floorDirection.dot(p-left);
        return along>0 && along<(right-left).norm() &&
               std::abs(up.dot(p-left))<=2*shortRadius;
      }
      return planes[0].normal.dot(p)-planes[0].offset>tolerance &&
             planes[1].normal.dot(p)-planes[1].offset<-tolerance;
    };
    // The measured body between supported sides is the weld face. Shallow
    // patches use a bounded toe-to-toe corridor: extending nearly horizontal
    // side planes would cut away the floor. A per-bin robust contour band
    // removes isolated floating samples without generating a synthetic plane.
    std::map<int,std::vector<float>> heights;
    for(const auto &p:*profile)if(inside(Eigen::Vector2f(p.x,p.y)))
      heights[int(std::floor((p.x-left.x())/binSize))].push_back(p.y);
    std::map<int,float> floor;
    for(const auto &bin:heights)floor.emplace(bin.first,median(bin.second));
    for(std::size_t j=0;j<indices.size();++j) {
      // Weld material can be rough: a stable wall normal is required for
      // assigning a side face, not for retaining supported points between it.
      if(!supported[j])continue;
      const auto &p=(*profile)[j];
      const Eigen::Vector2f point(p.x,p.y);
      const int side=p.x<middle?0:1;
      const auto &plane=planes[side];
      const Eigen::Vector2f tangent(-plane.normal.y(),plane.normal.x());
      const float t=tangent.dot(point);
      const bool wallNormal=plane.valid && std::abs(normals[j].dot(plane.normal))>.9f &&
                            (plane.toePatch || std::abs(normals[j].x())>minimumSlope);
      const auto &toe=side==0?left:right;
      const bool outsidePatch=!plane.toePatch ||
          (floorDirection.dot(point-toe)*(side==0?-1:1)>=0 &&
           (point-toe).norm()<=3*shortRadius);
      if(plane.valid && std::abs(plane.normal.dot(point)-plane.offset)<=2*tolerance &&
         wallNormal && outsidePatch &&
         t>=plane.begin-tolerance && t<=plane.end+tolerance) {
        labels[indices[j]]=side==0?1:2;
      } else if(inside(point)) {
        const bool outsideToe=side==0?p.x<toe.x()+tolerance:p.x>toe.x()-tolerance;
        const bool upperWall=config.interiorOpeningDirection &&
            (p.y-toe.y())*config.interiorOpeningDirection>2*radius;
        // A slightly bowed upper wall may miss the local plane's residual
        // band. It must not become weld material merely by falling inside it.
        if(outsideToe && upperWall)continue;
        const auto bin=floor.find(int(std::floor((p.x-left.x())/binSize)));
        if(bin!=floor.end() && std::abs(p.y-bin->second)<=2*tolerance)
          labels[indices[j]]=3;
      }
    }
  }
  // Preserve original order and exact XYZ, including legitimate repeated
  // samples. Distinct output clouds never share a source index.
  for (std::size_t i=0; i<labels.size(); ++i) {
    if (labels[i]==1) result.leftSurface->push_back((*source)[i]);
    else if (labels[i]==2) result.rightSurface->push_back((*source)[i]);
    else if (labels[i]==3) result.weldSurface->push_back((*source)[i]);
  }
  result.stageMilliseconds["surface_cloud_classification"]=
      std::chrono::duration<double,std::milli>(
          std::chrono::steady_clock::now()-classificationStarted).count();
  if (generateSurfaceModels) {
    const auto modelStarted=std::chrono::steady_clock::now();
    result.leftSurfaceModel=fitSurfaceModel(*source,groups,sections,labels,1,config);
    result.weldSurfaceModel=fitSurfaceModel(*source,groups,sections,labels,3,config);
    result.rightSurfaceModel=fitSurfaceModel(*source,groups,sections,labels,2,config);
    result.stageMilliseconds["surface_model_fitting"]=
        std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-modelStarted).count();
  }
}
}
