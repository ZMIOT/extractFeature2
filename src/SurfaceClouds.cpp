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

// Follow the measured floor inward from a supported wall foot when the other
// trajectory is absent. Each bin must connect to the preceding measured band;
// neither a nominal groove width nor a neighboring frame supplies its extent.
std::map<int,float> traceOpenFloor(const CloudT &profile,
    const std::vector<bool> &supported, const std::vector<Eigen::Vector2f> &normals,
    const Eigen::Vector2f &toe, int inward, float binSize, float radius, float tolerance) {
  std::map<int,std::vector<std::size_t>> bins;
  for(std::size_t j=0;j<profile.size();++j) {
    const float distance=inward*(profile[j].x-toe.x());
    if(supported[j] && distance>0)
      bins[int(std::floor(distance/binSize))].push_back(j);
  }
  std::map<int,float> floor;
  float lastDistance=0,lastHeight=toe.y();
  for(const auto &bin:bins) {
    const float distance=(bin.first+.5f)*binSize;
    if(distance-lastDistance>2*radius)break;
    float slope=0;
    if(floor.size()>1) {
      auto previous=floor.lower_bound(bin.first-int(std::ceil(3*radius/binSize)));
      if(previous!=floor.end()) {
        const float span=lastDistance-(previous->first+.5f)*binSize;
        if(span>=radius)slope=(lastHeight-previous->second)/span;
      }
    }
    const float predicted=lastHeight+slope*(distance-lastDistance);
    std::vector<float> heights;
    for(auto j:bin.second) {
      // A floor may be bowed or rough, but an upright return is its boundary.
      if(std::abs(normals[j].y())<.707107f)continue;
      if(std::abs(profile[j].y-predicted)<=2*tolerance)heights.push_back(profile[j].y);
    }
    if(heights.empty())continue;
    const float height=median(heights);
    // Check a resolved baseline rather than amplifying noise in adjacent bins.
    auto anchor=floor.lower_bound(bin.first-int(std::ceil(3*radius/binSize)));
    const float anchorDistance=anchor==floor.end()?0:(anchor->first+.5f)*binSize;
    const float anchorHeight=anchor==floor.end()?toe.y():anchor->second;
    if(std::abs(height-anchorHeight)>.7f*(distance-anchorDistance)+tolerance)break;
    floor.emplace(bin.first,height);
    lastDistance=distance;lastHeight=height;
  }
  if(floor.size()<3 || lastDistance<2*radius)floor.clear();
  return floor;
}

bool measuredWallToe(const CloudT &profile, const std::vector<bool> &supported,
    const std::vector<Eigen::Vector2f> &normals, const SidePlane &plane,
    int side, int opening, float radius, float tolerance, float toeReach, Eigen::Vector2f &toe) {
  if(!plane.valid || !opening)return false;
  const Eigen::Vector2f tangent(-plane.normal.y(),plane.normal.x());
  const Eigen::Vector2f lower=plane.normal*plane.offset+tangent*(opening>0?plane.begin:plane.end);
  float best=std::numeric_limits<float>::max();
  for(std::size_t j=0;j<profile.size();++j) {
    if(!supported[j] || std::abs(normals[j].y())<.866025f)continue;
    const Eigen::Vector2f point(profile[j].x,profile[j].y);
    const float distance=(plane.normal.dot(point)-plane.offset)*(side==0?1:-1);
    if(distance<-tolerance || distance>3*radius || (point-lower).norm()>toeReach)continue;
    if((point.y()-lower.y())*opening>tolerance)continue;
    const float score=std::abs(distance);
    if(score<best) {best=score;toe=point;}
  }
  return std::isfinite(best) && best<std::numeric_limits<float>::max();
}

SidePlane recoverMeasuredWall(std::vector<Eigen::Vector2f> seeds, const CloudT &profile,
    const std::vector<bool> &supported, const std::vector<Eigen::Vector2f> &normals,
    int side, int opening, float radius, float tolerance, float minimumSlope, float toeReach,
    Eigen::Vector2f &toe) {
  for(int attempt=0;attempt<3 && seeds.size()>=6;++attempt) {
    const auto plane=fitSidePlane(seeds,toe,tolerance,radius,minimumSlope,std::numeric_limits<float>::max());
    if(!plane.valid)return {};
    if(measuredWallToe(profile,supported,normals,plane,side,opening,radius,tolerance,toeReach,toe))return plane;
    // A detached wall may have more samples than the real face. Without a
    // measured floor-side contact it cannot bound the weld; try the next face.
    seeds.erase(std::remove_if(seeds.begin(),seeds.end(),[&](const Eigen::Vector2f &p) {
      return std::abs(plane.normal.dot(p)-plane.offset)<=2*tolerance;
    }),seeds.end());
  }
  return {};
}
}

SurfaceFootRecovery selectMeasuredFootContinuation(const CloudT &track,
    const std::vector<float> &stations, const std::vector<SurfaceEnd> &feet,
    std::size_t firstRejected, const DetectorConfig &config) {
  SurfaceFootRecovery selected;
  const auto count=track.size();
  if(firstRejected>=count || stations.size()!=count || feet.size()!=count)return selected;
  const auto begin=std::size_t(std::lower_bound(stations.begin(),stations.end(),
      stations[firstRejected]-config.endpointTrimWindow)-stations.begin());
  for(std::size_t i=0;i<begin;++i)selected.kept.push_back(i);
  const float radius=std::max(.5f,config.profileBinSize);
  std::size_t previous=count,distant=0;
  float distantBegin=0;bool displaced=false;
  for(auto i=begin;i<count;++i) {
    const auto &toe=feet[i];
    bool supported=toe.valid;
    if(supported && previous<count)
      supported=(toVector(toe.point)-toVector(feet[previous].point)).norm()<=
          3*(stations[i]-stations[previous])+2*radius;
    if(!supported) {
      if(i<firstRejected)selected.kept.push_back(i);
      if(previous==count || stations[i]-stations[previous]>.25f*config.endpointTrimWindow)break;
      distant=0;continue;
    }
    const float distance=(toVector(toe.point)-toVector(track[i])).norm();
    if(previous==count && distance>std::max(3*radius,3*config.endpointTrimDeviation))break;
    previous=i;selected.kept.push_back(i);selected.recovered.push_back(i);
    if(distance>std::max(2*radius,2*config.endpointTrimDeviation)) {
      if(!distant)distantBegin=stations[i];
      ++distant;
      if(stations[i]-distantBegin>=.25f*config.endpointTrimWindow)displaced=true;
    } else distant=0;
  }
  if(!displaced || previous==count || previous<=firstRejected ||
      stations[previous]-stations[firstRejected]<config.endpointTrimWindow)return {};
  return selected;
}

void extractSurfaceClouds(const CloudT::ConstPtr &source,
                          const std::vector<SurfaceSection> &sections,
                          const DetectorConfig &config, SeamResult &result,
                          const std::vector<ScanFrameRange> *frames,
                          bool generateSurfaceModels,
                          const std::array<SurfaceEnd,2> *ends,
                          std::vector<std::array<SurfaceEnd,2>> *measuredToes) {
  if(measuredToes)measuredToes->assign(sections.size(),{});
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
  std::vector<std::array<float,2>> floorBounds(sections.size(),
      {{-std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity()}});
  detail::GpuSurfaceNormals gpuNormalBatch;
  const bool gpuNormals=detail::computeSurfaceNormalsGpu(
      *source,groups,u,v,radius,gpuNormalBatch);
  if(gpuNormals)
    result.stageMilliseconds["gpu_surface_normals"]=gpuNormalBatch.milliseconds;
  for (std::size_t s=0; s<sections.size(); ++s) {
    const auto &section = sections[s];
    const auto &indices = groups[s];
    if (indices.size()<6 || section.left.data[u]>=section.right.data[u]) continue;
    Eigen::Vector2f left(section.left.data[u], section.left.data[v]);
    Eigen::Vector2f right(section.right.data[u], section.right.data[v]);
    float middle = .5f*(left.x()+right.x());
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
      const bool uncertain=isLeft?section.leftBoundaryUncertain:section.rightBoundaryUncertain;
      const auto &toe = isLeft ? left : right;
      Eigen::Vector2f wall = isLeft ? section.leftNormal : section.rightNormal;
      const bool knownWall = !uncertain && wall.squaredNorm()>.5f;
      if (knownWall) wall.normalize();
      // Broad wall support, including points far above the toe; local normal
      // gates exclude horizontal outer shoulders before fitting the big face.
      const float minimumSlope=std::max(std::sin(config.minSurfaceAngleDegrees*.01745329252f),
                                        knownWall ? .75f*std::abs(wall.x()) : uncertain ? .5f : .7f);
      const bool wallDirection = (!knownWall || std::abs(normal.dot(wall))>.5f) &&
                                 std::abs(normal.x())>minimumSlope;
      Eigen::Vector2f tangent(-wall.y(), wall.x());
      const int opening = config.interiorOpeningDirection;
      if ((opening && tangent.y()*opening<0) || (!opening && tangent.x()*(isLeft?-1:1)<0)) tangent = -tangent;
      const bool outward=isLeft ? p.x<=toe.x()+tolerance : p.x>=toe.x()-tolerance;
      const bool aboveToe=opening ? (p.y-toe.y())*opening>tolerance :
                                    knownWall && tangent.dot(point-toe)>tolerance;
      const Eigen::Vector2f oriented=normal.x()<0?-normal:normal;
      const bool risesOutward=oriented.y()*opening*(isLeft?1:-1)>=-.1f;
      if (wallDirection && outward && (uncertain?risesOutward:aboveToe)) {
        wallSeeds[isLeft?0:1].push_back(point);
      }
    }
    const float minimumSlope=std::sin(config.minSurfaceAngleDegrees*.01745329252f);
    const float maximumToeDistance=std::max(3*radius,config.profileCurvatureRadius);
    // A trimmed trajectory may have moved onto a weak interior ripple. In that
    // case use the current frame's broad rising wall, then locate its measured
    // floor-side foot. Never reuse the rejected point as a short patch anchor.
    std::array<SidePlane,2> planes;
    for(int side=0;side<2;++side) {
      auto &toe=side==0?left:right;
      const bool uncertain=side==0?section.leftBoundaryUncertain:section.rightBoundaryUncertain;
      planes[side]=uncertain
          ? recoverMeasuredWall(wallSeeds[side],*profile,supported,normals,side,config.interiorOpeningDirection,radius,tolerance,minimumSlope,std::max(4*radius,maximumToeDistance),toe)
          : fitSidePlane(wallSeeds[side],toe,tolerance,radius,minimumSlope,maximumToeDistance);
    }
    middle=.5f*(left.x()+right.x());
    if(measuredToes)for(int side=0;side<2;++side) {
      const bool uncertain=side?section.rightBoundaryUncertain:section.leftBoundaryUncertain;
      const bool valid=side?section.rightValid:section.leftValid;
      if(!uncertain || !valid || !planes[side].valid)continue;
      const auto &toe=side?right:left;
      for(std::size_t j=0;j<indices.size();++j) {
        if((Eigen::Vector2f((*profile)[j].x,(*profile)[j].y)-toe).squaredNorm()>1.e-8f)continue;
        auto &measured=(*measuredToes)[s][side];
        measured.valid=true;measured.station=section.station;
        measured.point=(*source)[indices[j]];
        break;
      }
    }
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
      if(side==0?section.leftBoundaryUncertain:section.rightBoundaryUncertain)continue;
      planes[side]=fitToePatch(*profile,supported,side==0?left:right,
          side==0?Eigen::Vector2f(-floorDirection):floorDirection,up,radius,tolerance);
      if(!planes[side].valid && shortRadius>radius)
        planes[side]=fitToePatch(*profile,supported,side==0?left:right,
            side==0?Eigen::Vector2f(-floorDirection):floorDirection,up,shortRadius,tolerance);
    }
    if(!section.leftValid)planes[0].valid=false;
    if(!section.rightValid)planes[1].valid=false;
    // A two-wall floor ends at the measured feet, even when the infinite
    // fitted wall planes would also admit an exterior sheet or shoulder.
    if(planes[0].valid)floorBounds[s][0]=left.x();
    if(planes[1].valid)floorBounds[s][1]=right.x();
    const bool recoveredSection=section.leftBoundaryUncertain || section.rightBoundaryUncertain;
    const bool singleSide=recoveredSection ? planes[0].valid!=planes[1].valid : section.leftValid!=section.rightValid;
    const int presentSide=recoveredSection? (planes[0].valid?0:1) : (section.leftValid?0:1);
    const auto &floorToe=presentSide==0?left:right;
    const int inward=presentSide==0?1:-1;
    const auto openFloor=(singleSide || recoveredSection) && planes[presentSide].valid
        ? traceOpenFloor(*profile,supported,normals,floorToe,inward,binSize,radius,tolerance)
        : std::map<int,float>{};
    const bool shortSection=planes[0].toePatch || planes[1].toePatch;
    auto inside=[&](const Eigen::Vector2f &p) {
      if(singleSide || recoveredSection) {
        if(!singleSide && (p.x()<=left.x() || p.x()>=right.x()))return false;
        const float distance=inward*(p.x()-floorToe.x());
        if(distance<=0)return false;
        const auto bin=openFloor.find(int(std::floor(distance/binSize)));
        return bin!=openFloor.end() && std::abs(p.y()-bin->second)<=2*tolerance;
      }
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
        if((singleSide || recoveredSection) && normals[j].squaredNorm()>.5f && std::abs(normals[j].y())<.707107f)continue;
        const bool outsideToe=side==0?p.x<toe.x()+tolerance:p.x>toe.x()-tolerance;
        const bool upperWall=config.interiorOpeningDirection &&
            (p.y-toe.y())*config.interiorOpeningDirection>2*radius;
        // A slightly bowed upper wall may miss the local plane's residual
        // band. It must not become weld material merely by falling inside it.
        if(!singleSide && outsideToe && upperWall)continue;
        const auto bin=floor.find(int(std::floor((p.x-left.x())/binSize)));
        if(bin!=floor.end() && std::abs(p.y-bin->second)<=2*tolerance)
          labels[indices[j]]=3;
      }
    }
  }
  // Recover raw measurements rejected only by unstable local normals. Use
  // fixed, same-section face seeds and interpolate within their contour span;
  // never grow recursively into shoulders or borrow another face's points.
  const auto seeds=labels;
  for(std::size_t s=0;s<groups.size();++s)for(unsigned char face=1;face<=3;++face) {
    if((face==1 && !sections[s].leftValid) || (face==2 && !sections[s].rightValid))continue;
    std::vector<Eigen::Vector2f> contour;
    for(auto index:groups[s])if(seeds[index]==face)
      contour.emplace_back((*source)[index].data[u],(*source)[index].data[v]);
    if(contour.size()<6)continue;
    Eigen::Vector2f center=Eigen::Vector2f::Zero();
    for(const auto &p:contour)center+=p;
    center/=float(contour.size());
    Eigen::Matrix2f covariance=Eigen::Matrix2f::Zero();
    for(const auto &p:contour) {const Eigen::Vector2f d=p-center;covariance.noalias()+=d*d.transpose();}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> fit(covariance);
    if(fit.info()!=Eigen::Success || fit.eigenvalues()[1]<1.e-6f)continue;
    const Eigen::Vector2f direction=fit.eigenvectors().col(1);
    std::sort(contour.begin(),contour.end(),[&](const Eigen::Vector2f &p,const Eigen::Vector2f &q) {
      return direction.dot(p)<direction.dot(q);
    });
    const float reach=std::max(2*radius,4*tolerance);
    for(auto index:groups[s])if(!labels[index]) {
      const Eigen::Vector2f p((*source)[index].data[u],(*source)[index].data[v]);
      const float t=direction.dot(p);
      if(t<=direction.dot(contour.front())+.25f*radius ||
         t>=direction.dot(contour.back())-.25f*radius)continue;
      auto next=std::lower_bound(contour.begin(),contour.end(),t,
          [&](const Eigen::Vector2f &q,float x){return direction.dot(q)<x;});
      if(next==contour.begin() || next==contour.end())continue;
      const auto &before=*(next-1);
      const Eigen::Vector2f delta=*next-before;
      const float span=direction.dot(delta);
      if(span<1.e-6f || delta.norm()>reach)continue;
      const Eigen::Vector2f expected=before+delta*((t-direction.dot(before))/span);
      if((p-expected).norm()<=std::min(tolerance,.25f*radius))labels[index]=face;
    }
  }
  // Use detector stations, not raw X or frame-vector offsets: original scans
  // may be tilted or reversed. Each wall keeps its own end; a two-wall floor
  // stops at their common end instead of projecting past the shorter wall.
  if(config.limitSurfaceEndsToTracks && ends) {
    const auto &left=(*ends)[0], &right=(*ends)[1];
    for(std::size_t s=0;s<groups.size();++s)for(auto index:groups[s]) {
      const auto label=labels[index];
      bool limited=false;
      float stop=0;
      if(label==1 && left.valid) {stop=left.station;limited=true;}
      if(label==2 && right.valid) {stop=right.station;limited=true;}
      if(label==3 && (left.valid || right.valid)) {
        limited=true;
        stop=left.valid?left.station:right.station;
        if(left.valid && right.valid)stop=std::min(left.station,right.station);
      }
      if(limited && sections[s].station>stop+1.e-4f)labels[index]=0;
    }
  }
  // On a two-wall workpiece, the floor only exists over the common measured
  // wall span. Apply this after recovery and endpoint clipping: otherwise a
  // single-side floor trace can leave a tongue beyond either wall at either
  // end. An actual single-wall scan keeps its connected measured floor.
  std::array<std::size_t,2> first{{groups.size(),groups.size()}},last{{0,0}};
  for(std::size_t s=0;s<groups.size();++s)for(auto index:groups[s]) {
    const auto label=labels[index];
    if(label==1 || label==2) {
      first[label-1]=std::min(first[label-1],s);
      last[label-1]=std::max(last[label-1],s);
    }
  }
  const bool twoWalls=first[0]<groups.size() && first[1]<groups.size();
  const auto floorBegin=std::max(first[0],first[1]),floorEnd=std::min(last[0],last[1]);
  if(twoWalls)for(int side=0;side<2;++side) {
    std::size_t previous=groups.size();
    for(std::size_t s=0;s<groups.size();++s)if(std::isfinite(floorBounds[s][side])) {
      // Only bridge a short loss of wall support inside two measured feet.
      // Never extrapolate a missing side from the other side's open floor.
      if(previous<groups.size() && s>previous+1) {
        const float span=sections[s].station-sections[previous].station;
        if(span>0 && span<=std::max(4.5f*config.sliceStep,config.profileCurvatureRadius))
          for(auto i=previous+1;i<s;++i) {
            const float t=(sections[i].station-sections[previous].station)/span;
            floorBounds[i][side]=(1-t)*floorBounds[previous][side]+t*floorBounds[s][side];
          }
      }
      previous=s;
    }
  }
  for(std::size_t s=0;s<groups.size();++s)for(auto index:groups[s])if(labels[index]==3) {
    const float across=(*source)[index].data[u];
    if((twoWalls && (s<floorBegin || s>floorEnd ||
         !std::isfinite(floorBounds[s][0]) || !std::isfinite(floorBounds[s][1]))) ||
        across<=floorBounds[s][0] || across>=floorBounds[s][1])labels[index]=0;
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
