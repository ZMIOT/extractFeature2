#include "SeamDirection.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <pcl/common/point_tests.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace weld {
namespace {
double median(std::vector<double> values) {
  const auto mid = values.begin() + values.size()/2;
  std::nth_element(values.begin(), mid, values.end());
  return *mid;
}
// Equal transverse-bin weights avoid mistaking dense scan margins for the floor.
// A quadratic term measures concavity after removing offset and linear tilt.
int openingVote(const std::vector<Eigen::Vector2d> &points, double minDepth) {
  if (points.size() < 20) return 0;
  double low=points.front().x(), high=low;
  for (const auto &p:points) {low=std::min(low,p.x());high=std::max(high,p.x());}
  if (high-low < 1.e-6) return 0;
  std::array<std::vector<double>,41> bins;
  for (const auto &p:points) {
    const int i=std::min(40,int(41*(p.x()-low)/(high-low)));
    bins[i].push_back(p.y());
  }
  std::vector<Eigen::Vector3d> rows;
  std::vector<double> heights;
  for (int i=0;i<41;++i) if (!bins[i].empty()) {
    const double u=2*(i+.5)/41-1;
    rows.emplace_back(1,u,u*u);heights.push_back(median(bins[i]));
  }
  if (rows.size()<12) return 0;
  const double origin=median(heights);
  for (auto &v:heights) v-=origin;
  std::vector<double> weights(rows.size(),1);
  Eigen::Vector3d fit=Eigen::Vector3d::Zero();
  double residual=0;
  for (int pass=0;pass<3;++pass) {
    Eigen::Matrix3d normal=Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs=Eigen::Vector3d::Zero();
    for (std::size_t i=0;i<rows.size();++i) {
      normal+=weights[i]*rows[i]*rows[i].transpose();
      rhs+=weights[i]*rows[i]*heights[i];
    }
    fit=normal.ldlt().solve(rhs);
    if (!fit.allFinite()) return 0;
    std::vector<double> errors;
    for (std::size_t i=0;i<rows.size();++i) errors.push_back(std::abs(heights[i]-rows[i].dot(fit)));
    residual=1.4826*median(errors);
    for (std::size_t i=0;i<rows.size();++i)
      weights[i]=std::min(1.,std::max(.01,2*residual)/std::max(1.e-9,errors[i]));
  }
  if (std::abs(fit.z())<std::max(minDepth,2*residual)) return 0;
  // A cropped single shoulder can have curvature of the opposite sign.
  // Require a supported interior extremum, with room for both flanks.
  if (std::abs(fit.y())>1.4*std::abs(fit.z())) return 0;
  return fit.z()>0 ? 1 : -1;
}
}

SeamDirectionEstimate estimateSeamDirection(const CloudT &cloud,
    const DetectorConfig &config, const std::vector<ScanFrameRange> *frames) {
  SeamDirectionEstimate result;result.axis=config.longitudinalAxis;
  Eigen::Vector3d center=Eigen::Vector3d::Zero();
  Eigen::Matrix3d scatter=Eigen::Matrix3d::Zero();
  std::size_t count=0;
  // Welford covariance in double precision remains stable far from the origin.
  for (const auto &p:cloud) if (pcl::isFinite(p)) {
    const Eigen::Vector3d v=toVector(p).cast<double>();
    const Eigen::Vector3d delta=v-center;
    center+=delta/double(++count);
    scatter+=delta*(v-center).transpose();
  }
  if (count>=3) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(scatter/double(count));
    if (solver.info()==Eigen::Success && solver.eigenvalues()[2]>1.e-12) {
      result.principalDirection=solver.eigenvectors().col(2).cast<float>();
      result.principalGap=float((solver.eigenvalues()[2]-solver.eigenvalues()[1])/solver.eigenvalues()[2]);
      bool frameProgressionReliable=false;
      // A short scan may be wider than it is long, and its scan plane may be
      // oblique to travel. PCA of retained frame centers measures progression,
      // rather than confusing cross-section extent or plane normal with travel.
      if (frames && frames->size()>=4) {
        std::vector<Eigen::Vector3d> centers,velocities;
        std::vector<int> frameIds;
        const auto n=std::min<std::size_t>(17,frames->size());
        for (std::size_t j=0;j<n;++j) {
          const auto &f=(*frames)[j*(frames->size()-1)/(n-1)];
          if (f.count<20 || f.begin>cloud.size() || f.count>cloud.size()-f.begin) continue;
          Eigen::Vector3d frameCenter=Eigen::Vector3d::Zero();std::size_t samples=0;
          for (std::size_t i=f.begin;i<f.begin+f.count;++i) if (pcl::isFinite(cloud[i])) {
            frameCenter+=weld::toVector(cloud[i]).cast<double>();++samples;
          }
          if (samples<20) continue;
          centers.push_back(frameCenter/double(samples));frameIds.push_back(f.frameIndex);
        }
        if (centers.size()>=3) {
          Eigen::Vector3d mean=Eigen::Vector3d::Zero();
          for (const auto &p:centers) mean+=p;
          mean/=double(centers.size());
          Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();
          for (const auto &p:centers) {const Eigen::Vector3d delta=p-mean;covariance+=delta*delta.transpose();}
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> progression(covariance);
          if (progression.info()==Eigen::Success && progression.eigenvalues()[2]>1.e-12 &&
              progression.eigenvalues()[1]<.2*progression.eigenvalues()[2]) {
            result.principalDirection=progression.eigenvectors().col(2).cast<float>();
            result.principalGap=float(1-progression.eigenvalues()[1]/progression.eigenvalues()[2]);
            frameProgressionReliable=true;
          }
        }
        for (std::size_t i=1;i<centers.size();++i) {
          const double gap=double(frameIds[i])-double(frameIds[i-1]);
          if (gap>0) velocities.push_back((centers[i]-centers[i-1])/gap);
        }
        if (velocities.size()>=3) {
          Eigen::Vector3d typical;
          for (int a=0;a<3;++a) {
            std::vector<double> components;
            for (const auto &v:velocities) components.push_back(v[a]);
            typical[a]=median(components);
          }
          std::vector<double> deviations;
          for (const auto &v:velocities) deviations.push_back((v-typical).norm());
          const double gate=std::max(1.e-6,std::max(3*median(deviations),.1*typical.norm()));
          Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();int retained=0;
          // A sudden wall-height/coverage change moves the center without
          // moving the seam. Reject such displacement spikes before the PCA.
          // Use directions about the origin, so constant travel is not removed.
          for (std::size_t i=0;i<velocities.size();++i) if (deviations[i]<=gate) {
            covariance+=velocities[i]*velocities[i].transpose();++retained;
          }
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> progression(covariance);
          if (retained>=3 && progression.info()==Eigen::Success && progression.eigenvalues()[2]>1.e-12 &&
              progression.eigenvalues()[1]<.2*progression.eigenvalues()[2] &&
              (!frameProgressionReliable || (retained<int(velocities.size()) &&
               progression.eigenvalues()[1]<.001*progression.eigenvalues()[2]))) {
            result.principalDirection=progression.eigenvectors().col(2).cast<float>();
            result.principalGap=float(1-progression.eigenvalues()[1]/progression.eigenvalues()[2]);
            frameProgressionReliable=true;
          }
        }
      }
      Eigen::Index dominant=0;
      result.principalDirection.cwiseAbs().maxCoeff(&dominant);
      if (result.principalDirection[dominant]<0) result.principalDirection*=-1;
      result.axisReliable=result.principalGap>=.1f &&
          (!frames || frames->empty() || frameProgressionReliable);
      if (config.autoLongitudinalAxis && result.axisReliable) result.axis=Axis(dominant);
    }
  }
  if (!config.autoOpeningDirection) return result;
  const int a=int(result.axis),u=a==0?1:0,v=a==2?1:2;
  std::array<std::vector<Eigen::Vector2d>,17> sections;
  if (frames && !frames->empty()) {
    const auto n=std::min<std::size_t>(sections.size(),frames->size());
    for (std::size_t j=0;j<n;++j) {
      const auto &f=(*frames)[n==1?0:j*(frames->size()-1)/(n-1)];
      if (f.begin>cloud.size() || f.count>cloud.size()-f.begin) continue;
      for (std::size_t i=f.begin;i<f.begin+f.count;++i) if (pcl::isFinite(cloud[i]))
        sections[j].emplace_back(cloud[i].data[u],cloud[i].data[v]);
    }
  } else {
    double low=std::numeric_limits<double>::max(),high=-low;
    for (const auto &p:cloud) if (pcl::isFinite(p)) {low=std::min(low,double(p.data[a]));high=std::max(high,double(p.data[a]));}
    if (high-low>1.e-6) for (const auto &p:cloud) if (pcl::isFinite(p)) {
      const int i=std::min(16,int(17*(p.data[a]-low)/(high-low)));
      sections[i].emplace_back(p.data[u],p.data[v]);
    }
  }
  for (const auto &section:sections) {
    const int vote=openingVote(section,std::max(.1,double(config.minimumFeatureDepth)*.25));
    result.positiveVotes+=vote>0;result.negativeVotes+=vote<0;
  }
  const int votes=result.positiveVotes+result.negativeVotes;
  result.openingReliable=votes>=3 && 3*std::max(result.positiveVotes,result.negativeVotes)>=2*votes;
  if (result.openingReliable) result.opening=result.positiveVotes>result.negativeVotes?1:-1;
  return result;
}
}
