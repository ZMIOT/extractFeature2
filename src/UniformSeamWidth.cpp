#include "UniformSeamWidth.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

namespace weld {
namespace {
struct Arc {
  const std::vector<PointT> *points=nullptr;
  std::vector<double> lengths;
  int axis=1;
  double center=0;
  bool init(const UniformWidthSection &section, int transverse) {
    points=&section.contour; axis=transverse;
    if (points->size()<2 || !std::isfinite(section.bottomCenterU)) return false;
    lengths.push_back(0);
    for (std::size_t i=0;i<points->size();++i) {
      if (!toVector((*points)[i]).allFinite()) return false;
      if (i) {
        if ((*points)[i].data[axis]<=(*points)[i-1].data[axis]) return false;
        lengths.push_back(lengths.back()+(toVector((*points)[i])-toVector((*points)[i-1])).cast<double>().norm());
      }
    }
    center=atU(section.bottomCenterU);
    return center>0 && center<lengths.back();
  }
  double atU(double u) const {
    if (u<points->front().data[axis] || u>points->back().data[axis]) return -1;
    auto it=std::lower_bound(points->begin(),points->end(),u,[&](const PointT &p,double v){return p.data[axis]<v;});
    if (it==points->begin()) return 0;
    if (it==points->end()) return lengths.back();
    const auto i=std::size_t(it-points->begin());
    const double t=(u-(*points)[i-1].data[axis])/((*points)[i].data[axis]-(*points)[i-1].data[axis]);
    return lengths[i-1]+t*(lengths[i]-lengths[i-1]);
  }
  Eigen::Vector3d at(double s) const {
    auto it=std::lower_bound(lengths.begin(),lengths.end(),s);
    if (it==lengths.begin()) return toVector(points->front()).cast<double>();
    if (it==lengths.end()) return toVector(points->back()).cast<double>();
    const auto i=std::size_t(it-lengths.begin());
    const double t=(s-lengths[i-1])/(lengths[i]-lengths[i-1]);
    return (1-t)*toVector((*points)[i-1]).cast<double>()+t*toVector((*points)[i]).cast<double>();
  }
  bool place(double width, UniformWidthMetric metric, PointT &left, PointT &right) const {
    const double available=std::min(center,lengths.back()-center);
    double half=width/2;
    if (metric==UniformWidthMetric::Chord) {
      // Each knot interval has affine endpoints. The norm is convex there;
      // the first crossing is found without assuming global monotonicity.
      std::vector<double> knots{0,available};
      for (double s:lengths) if (std::abs(s-center)<available) knots.push_back(std::abs(s-center));
      std::sort(knots.begin(),knots.end());
      auto distance=[&](double r){return (at(center+r)-at(center-r)).norm();};
      bool found=false;
      for (std::size_t i=1;i<knots.size();++i) if (distance(knots[i])>=width) {
        double low=knots[i-1],high=knots[i];
        for (int j=0;j<40;++j) { const double mid=(low+high)/2; if(distance(mid)<width)low=mid;else high=mid; }
        half=(low+high)/2; found=true; break;
      }
      if (!found) return false;
    }
    if (half>available || half<=0) return false;
    left=toPoint(at(center-half).cast<float>()); right=toPoint(at(center+half).cast<float>());
    return left.data[axis]<at(center)[axis] && right.data[axis]>at(center)[axis];
  }
};
template<class T> void select(std::vector<T> &values,const std::vector<std::size_t> &keep) {
  if(values.empty())return;
  std::vector<T> out; out.reserve(keep.size());
  for(auto i:keep) if(i<values.size()) out.push_back(values[i]);
  values=std::move(out);
}

// Jointly filter the pair center and half-width vector in physical station
// coordinates. Normalize the fitted vector only afterwards, so a fixed chord
// width survives smoothing. Never project back onto the noisy bottom samples.
void filterUniformPairs(SeamResult &result,const DetectorConfig &config,
                        const std::vector<UniformWidthSection> &sections,
                        const std::vector<std::size_t> &keep) {
  if(!config.uniformWidthFilterEnabled || config.uniformWidthMetric!=UniformWidthMetric::Chord || keep.size()<3)return;
  const auto count=keep.size();
  auto starts=result.leftSegmentStarts;
  starts.push_back(0);
  std::vector<double> gaps;
  for(std::size_t i=1;i<count;++i) {
    const double gap=sections[keep[i]].station-sections[keep[i-1]].station;
    if(std::isfinite(gap) && gap>1.e-6)gaps.push_back(gap);
    else starts.push_back(i);
  }
  if(gaps.empty())return;
  std::sort(gaps.begin(),gaps.end());
  const double maxGap=gaps[gaps.size()/2]*std::max(2.,double(config.maxInterpolatedGapStations)+1.5);
  for(std::size_t i=1;i<count;++i)
    if(sections[keep[i]].station-sections[keep[i-1]].station>maxGap)starts.push_back(i);
  starts.push_back(count);std::sort(starts.begin(),starts.end());
  starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
  result.leftSegmentStarts=result.rightSegmentStarts=std::vector<std::size_t>(starts.begin(),starts.end()-1);
  const int transverse=config.longitudinalAxis==Axis::X?1:0;
  const double lambda=std::pow(double(config.uniformWidthFilterLength)/(2*3.141592653589793),4);
  double totalCorrection=0;
  for(std::size_t segment=1;segment<starts.size();++segment) {
    const auto begin=starts[segment-1];const int n=int(starts[segment]-begin);
    if(n<3)continue;
    Eigen::MatrixXd observed(n,6),fitted;
    Eigen::VectorXd weights(n),base(n);
    const Eigen::Vector3d origin=.5*(toVector((*result.leftSeam)[begin]).cast<double>()+toVector((*result.rightSeam)[begin]).cast<double>());
    for(int j=0;j<n;++j) {
      const auto i=begin+j;
      const Eigen::Vector3d l=toVector((*result.leftSeam)[i]).cast<double>(),r=toVector((*result.rightSeam)[i]).cast<double>();
      observed.block<1,3>(j,0)=(.5*(l+r)-origin).transpose();
      observed.block<1,3>(j,3)=(.5*(r-l)).transpose();
      const double score=i<result.leftPointConfidence.size() && i<result.rightPointConfidence.size() ?
          std::min(result.leftPointConfidence[i],result.rightPointConfidence[i]):1.;
      base[j]=std::isfinite(score)?std::max(.1,std::min(1.,score)):.1;
    }
    std::vector<Eigen::Triplet<double>> penalty;
    for(int j=1;j+1<n;++j) {
      const double h0=sections[keep[begin+j]].station-sections[keep[begin+j-1]].station;
      const double h1=sections[keep[begin+j+1]].station-sections[keep[begin+j]].station;
      const double d[3]={2/(h0*(h0+h1)),-2/(h0*h1),2/(h1*(h0+h1))};
      for(int row=0;row<3;++row)for(int col=0;col<3;++col)
        penalty.emplace_back(j+row-1,j+col-1,lambda*d[row]*d[col]);
    }
    weights=base;bool solved=false;
    for(int pass=0;pass<8;++pass) {
      auto entries=penalty;for(int j=0;j<n;++j)entries.emplace_back(j,j,weights[j]);
      Eigen::SparseMatrix<double> system(n,n);system.setFromTriplets(entries.begin(),entries.end());
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(system);
      if(solver.info()!=Eigen::Success){solved=false;break;}
      fitted=solver.solve(weights.asDiagonal()*observed);
      solved=solver.info()==Eigen::Success && fitted.allFinite();if(!solved)break;
      for(int j=0;j<n;++j)weights[j]=base[j]*std::min(1.,.6/std::max(1.e-9,(fitted.row(j)-observed.row(j)).norm()));
    }
    if(!solved)continue;
    result.uniformWidthFilterApplied=true;
    for(int j=0;j<n;++j) {
      const auto i=begin+j;const auto &section=sections[keep[i]];
      const Eigen::Vector3d oldL=toVector((*result.leftSeam)[i]).cast<double>(),oldR=toVector((*result.rightSeam)[i]).cast<double>();
      Eigen::Vector3d center=origin+fitted.block<1,3>(j,0).transpose();
      Eigen::Vector3d direction=fitted.block<1,3>(j,3).transpose();
      // Preserve each original scan plane, even when world X is not monotonic.
      Eigen::Vector3d mean=Eigen::Vector3d::Zero();
      for(const auto &p:section.contour)mean+=toVector(p).cast<double>();
      mean/=section.contour.size();Eigen::Matrix3d cov=Eigen::Matrix3d::Zero();
      for(const auto &p:section.contour){const Eigen::Vector3d d=toVector(p).cast<double>()-mean;cov+=d*d.transpose();}
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
      if(eig.info()==Eigen::Success && eig.eigenvalues()[1]>1.e-8) {
        const Eigen::Vector3d normal=eig.eigenvectors().col(0);
        center-=normal*normal.dot(center-mean);direction-=normal*normal.dot(direction);
      } else {
        // A straight contour only supports translation along that line.
        const Eigen::Vector3d tangent=(oldR-oldL).normalized();
        center=.5*(oldL+oldR)+tangent*tangent.dot(center-.5*(oldL+oldR));
        direction=tangent;
      }
      if(direction.norm()<1.e-6 || direction[transverse]<=0)continue;
      direction*=.5*result.uniformWidthTarget/direction.norm();
      const Eigen::Vector3d l=center-direction,r=center+direction;
      if(l[transverse]>=section.bottomCenterU || r[transverse]<=section.bottomCenterU)continue;
      (*result.leftSeam)[i]=toPoint(l.cast<float>());(*result.rightSeam)[i]=toPoint(r.cast<float>());
      const double dl=(l-oldL).norm(),dr=(r-oldR).norm();totalCorrection+=dl+dr;
      result.uniformWidthFilterMaxCorrection=std::max(result.uniformWidthFilterMaxCorrection,float(std::max(dl,dr)));
      for(auto *quality:{&result.leftPointQuality,&result.rightPointQuality})if(i<quality->size())
        (*quality)[i].flags|=QualityFinalCurveSmoothed;
    }
  }
  result.uniformWidthFilterMeanCorrection=float(totalCorrection/(2*count));
}
}

void applyUniformSeamWidth(SeamResult &result,const DetectorConfig &config,
                           const std::vector<UniformWidthSection> &sections) {
  if(!config.uniformWidthEnabled)return;
  const auto count=result.leftSeam->size();
  const int axis=config.longitudinalAxis==Axis::X?1:0;
  std::vector<Arc> arcs(count);
  std::vector<bool> valid(count,false);
  std::vector<double> widths;
  for(std::size_t i=0;i<count && i<sections.size();++i) {
    if(!arcs[i].init(sections[i],axis))continue;
    const double left=arcs[i].atU(sections[i].sourceLeft.data[axis]);
    const double right=arcs[i].atU(sections[i].sourceRight.data[axis]);
    if(left<0 || right<=left || left>=arcs[i].center || right<=arcs[i].center)continue;
    const double width=config.uniformWidthMetric==UniformWidthMetric::BottomArc ? right-left :
        (toVector(sections[i].sourceRight)-toVector(sections[i].sourceLeft)).cast<double>().norm();
    if(!std::isfinite(width)||width<=0)continue;
    widths.push_back(width);valid[i]=true;
  }
  std::sort(widths.begin(),widths.end());
  result.uniformWidthSamples=widths.size();
  const auto trim=widths.size()/10;
  result.uniformWidthTrimEach=trim;
  const double target=widths.empty()?0:std::accumulate(widths.begin()+trim,widths.end()-trim,0.)/(widths.size()-2*trim);
  result.uniformWidthTarget=float(target);
  std::vector<std::size_t> keep,starts;
  CloudT left,right;
  bool split=true;
  for(std::size_t i=0;i<count;++i) {
    if(std::binary_search(result.leftSegmentStarts.begin(),result.leftSegmentStarts.end(),i) ||
       std::binary_search(result.rightSegmentStarts.begin(),result.rightSegmentStarts.end(),i))split=true;
    PointT l,r;
    if(!valid[i]||!arcs[i].place(target,config.uniformWidthMetric,l,r)){split=true;continue;}
    if(split)starts.push_back(keep.size());
    split=false;keep.push_back(i);left.push_back(l);right.push_back(r);
  }
  *result.leftSeam=std::move(left);*result.rightSeam=std::move(right);
  result.uniformWidthAdjusted=keep.size();result.uniformWidthOmitted=count-keep.size();
  select(result.leftFrameIndices,keep);select(result.rightFrameIndices,keep);
  select(result.leftPointConfidence,keep);select(result.rightPointConfidence,keep);
  select(result.leftPointQuality,keep);select(result.rightPointQuality,keep);
  result.leftSegmentStarts=result.rightSegmentStarts=starts;
  filterUniformPairs(result,config,sections,keep);
  // The endpoints now belong to the bottom contour, not the old wall fit.
  for(auto *quality:{&result.leftPointQuality,&result.rightPointQuality})for(auto &q:*quality) {
    q.flags|=QualityUniformWidth|QualityWallUnavailable|QualityFloorUnavailable;
    q.wallAccepted=q.floorAccepted=q.trackedWall=false;q.supportDistance=-1;
  }
  auto mean=[](const std::vector<float>&v){return v.empty()?0.f:std::accumulate(v.begin(),v.end(),0.f)/v.size();};
  result.leftConfidence=mean(result.leftPointConfidence);result.rightConfidence=mean(result.rightPointConfidence);
  result.confidence=.5f*(result.leftConfidence+result.rightConfidence);
}
}
