#include "ProfileCurvature.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <limits>

namespace weld {

int wallSupportedTurnCentroid(const std::vector<ProfileCurvatureSample> &samples,
    const Eigen::Vector2f &center, const Eigen::Vector2f &normal, float offset,
    float tolerance, float radius, int opening, float midpoint, int side) {
  if (samples.size()<3 || !center.allFinite() || !normal.allFinite() ||
      !std::isfinite(offset) || !(tolerance>0) || !(radius>0)) return -1;
  if (!opening) opening=samples[samples.size()/2].point.y() <=
      .5f*(samples.front().point.y()+samples.back().point.y()) ? 1 : -1;
  std::vector<double> arc(samples.size(),0);
  std::vector<bool> supported(samples.size(),false);
  double mass=0, firstMoment=0;
  for (const auto &sample:samples) if (!sample.point.allFinite()) return -1;
  for (std::size_t i=1;i<samples.size();++i)
    arc[i]=arc[i-1]+(samples[i].point-samples[i-1].point).norm();
  for (std::size_t i=1;i+1<samples.size();++i) {
    const auto &p=samples[i];
    const float distance=(p.point-center).norm();
    if (!p.valid || !samples[i-1].valid || !samples[i+1].valid ||
        !std::isfinite(p.curvature) || p.curvature*opening<=0 || distance>=radius ||
        (side ? p.point.x()<=midpoint : p.point.x()>=midpoint) ||
        std::fabs(normal.dot(p.point)-offset)>tolerance) continue;
    supported[i]=true;
    const double taper=1-distance*distance/(radius*radius);
    const double weight=p.curvature*opening * .5*(arc[i+1]-arc[i-1])*taper*taper;
    mass+=weight; firstMoment+=weight*arc[i];
  }
  if (mass<1.e-6) return -1;
  const double target=firstMoment/mass;
  int best=-1; double distance=std::numeric_limits<double>::infinity();
  for (std::size_t i=1;i+1<samples.size();++i)
    if (supported[i] && std::fabs(arc[i]-target)<distance) {
      best=int(i); distance=std::fabs(arc[i]-target);
    }
  return best;
}

int supportedProfileTurnRepresentative(
    const std::vector<ProfileCurvatureSample> &samples, const ProfileTurnRegion &region,
    const Eigen::Vector2f &normal, float offset, const Eigen::Vector2f &foot,
    float wallTolerance, float footRadius) {
  if (region.begin < 0 || region.end >= int(samples.size()) ||
      region.representative < region.begin || region.representative > region.end ||
      !normal.allFinite() || !foot.allFinite() || !std::isfinite(offset) ||
      !(wallTolerance >= 0) || !(footRadius > 0)) return -1;
  const auto &nominal = samples[region.representative];
  int best = -1;
  float distance = std::numeric_limits<float>::infinity();
  for (int i=region.begin; i<=region.end; ++i) {
    const auto &sample = samples[i];
    if (!sample.valid || !sample.point.allFinite() ||
        sample.curvature*nominal.curvature <= 0 ||
        std::fabs(normal.dot(sample.point)-offset) > wallTolerance ||
        (sample.point-foot).norm() > footRadius) continue;
    const float d = (sample.point-nominal.point).squaredNorm();
    if (d < distance) { best=i; distance=d; }
  }
  return best;
}

float profileBottomMidpoint(const std::vector<Eigen::Vector2f> &points,
                            int openingDirection, float heightTolerance) {
  if (points.empty()) return std::numeric_limits<float>::quiet_NaN();
  const int opening = openingDirection ? openingDirection :
      (points[points.size()/2].y() <= .5f*(points.front().y()+points.back().y()) ? 1 : -1);
  std::vector<float> height(points.size(), std::numeric_limits<float>::infinity());
  for (int i=0;i<int(points.size());++i) {
    if (!points[i].allFinite()) continue;
    std::vector<float> local;
    for (int j=std::max(0,i-2);j<=std::min(int(points.size())-1,i+2);++j)
      if (points[j].allFinite()) local.push_back(opening*points[j].y());
    std::sort(local.begin(),local.end());
    height[i]=local[local.size()/2];
  }
  const float minimum=*std::min_element(height.begin(),height.end());
  if (!std::isfinite(minimum)) return std::numeric_limits<float>::quiet_NaN();
  float maximum=minimum;
  for (float h:height) if (std::isfinite(h)) maximum=std::max(maximum,h);
  // Do not let a noise tolerance swallow an entire shallow notch or bead.
  const float tolerance=std::min(std::max(.01f,heightTolerance),
                                  std::max(.001f,.25f*(maximum-minimum)));
  float left=std::numeric_limits<float>::infinity(), right=-left;
  for (std::size_t i=0;i<points.size();++i)
    if (std::isfinite(height[i]) && height[i]<=minimum+tolerance) {
      left=std::min(left,points[i].x()); right=std::max(right,points[i].x());
    }
  return .5f*(left+right);
}

bool profilePairStraddlesBottom(float leftU, float rightU, float midpointU) {
  return std::isfinite(leftU) && std::isfinite(rightU) && std::isfinite(midpointU) &&
         leftU < midpointU && rightU > midpointU;
}

std::vector<unsigned> supplementalProfileCandidates(
    const std::vector<ProfileCurvatureSample> &samples, float windowLength,
    float prominence, int openingDirection) {
  const int n = static_cast<int>(samples.size());
  std::vector<unsigned> kinds(n, 0);
  if (n < 5 || !std::isfinite(windowLength) || windowLength <= 0 ||
      !std::isfinite(prominence) || prominence <= 0) return kinds;
  std::vector<float> arc(n, 0), change(n, 0);
  std::vector<int> left(n, -1), right(n, -1);
  for (int i = 1; i < n; ++i)
    arc[i] = arc[i-1] + (samples[i].point-samples[i-1].point).norm();
  for (int i = 1; i+1 < n; ++i) {
    if (!samples[i].valid) continue;
    int l = i, r = i;
    while (l > 0 && samples[l-1].valid && arc[i]-arc[l] < windowLength*.5f) --l;
    while (r+1 < n && samples[r+1].valid && arc[r]-arc[i] < windowLength*.5f) ++r;
    if (l == i || r == i || arc[i]-arc[l] < windowLength*.5f ||
        arc[r]-arc[i] < windowLength*.5f) continue;
    left[i] = l; right[i] = r;
    change[i] = std::fabs(samples[r].curvature-samples[l].curvature) /
                std::max(1.e-6f, arc[r]-arc[l]);
    for (int sign : {-1, 1}) {
      if (openingDirection && sign != openingDirection) continue;
      const float v = sign*samples[i].point.y();
      if (sign*samples[l].point.y()-v < prominence ||
          sign*samples[r].point.y()-v < prominence) continue;
      bool minimum = true;
      for (int j = l; j <= r; ++j)
        if (sign*samples[j].point.y() < v-1.e-6f ||
            (j < i && std::fabs(sign*samples[j].point.y()-v) <= 1.e-6f)) {
          minimum = false; break;
        }
      if (minimum) kinds[i] |= 2;
    }
  }
  for (int i = 1; i+1 < n; ++i) {
    if (left[i] < 0 || change[i] < 1.e-5f) continue;
    bool maximum = true;
    for (int j = left[i]; j <= right[i]; ++j)
      if (change[j] > change[i]+1.e-6f ||
          (j < i && std::fabs(change[j]-change[i]) <= 1.e-6f)) {
        maximum = false; break;
      }
    if (maximum) kinds[i] |= 4;
  }
  return kinds;
}
std::vector<bool> wallProfileSpikeMask(const std::vector<Eigen::Vector2f> &points,
                                      float radius, float maximumWidth, float tolerance) {
  std::vector<bool> removed(points.size(),false);
  if (points.size()<9 || radius<=0 || maximumWidth<=0 || tolerance<=0) return removed;
  for (const auto &p:points) if (!p.allFinite()) return removed;
  for (std::size_t i=1;i<points.size();++i) if (points[i].x()<points[i-1].x()) return removed;
  for (int pass=0;pass<3;++pass) {
    std::vector<std::size_t> active;
    for (std::size_t i=0;i<points.size();++i) if (!removed[i]) active.push_back(i);
    auto next=removed; bool changed=false;
    for (std::size_t k=1;k+1<active.size();++k) {
      const auto &q=points[active[k]];
      const float before=points[active[k-1]].y(),after=points[active[k+1]].y();
      if (q.y()>=std::min(before,after)-.5f*tolerance &&
          q.y()<=std::max(before,after)+.5f*tolerance) continue;
      // The wall may curve gradually: first check a single excursion against
      // two independent short flank lines instead of one long straight line.
      if (k>=3 && k+3<active.size() &&
          points[active[k+1]].x()-points[active[k-1]].x()<=maximumWidth) {
        double slopes[2]={}; bool supported=true;
        for (int side=0;side<2;++side) {
          double sx=0,sy=0,sxx=0,sxy=0;
          const auto first=side ? k+1 : k-3;
          for (auto j=first;j<first+3;++j) {
            const double x=points[active[j]].x()-q.x(),y=points[active[j]].y()-q.y();
            sx+=x;sy+=y;sxx+=x*x;sxy+=x*y;
          }
          if (3*sxx-sx*sx<1.e-9) {supported=false;break;}
          slopes[side]=(3*sxy-sx*sy)/(3*sxx-sx*sx);
          const double offset=(sy-slopes[side]*sx)/3,norm=std::hypot(1.0,slopes[side]);
          double error=0;
          for (auto j=first;j<first+3;++j) {
            const double d=(points[active[j]].y()-q.y()-offset-slopes[side]*(points[active[j]].x()-q.x()))/norm;
            error+=d*d;
          }
          if (std::fabs(slopes[side])/norm<.45 || std::sqrt(error/3)>.5*tolerance) supported=false;
        }
        const auto chord=points[active[k+1]]-points[active[k-1]];
        const auto delta=q-points[active[k-1]];
        const float distance=std::fabs(chord.x()*delta.y()-chord.y()*delta.x())/std::max(1.e-6f,chord.norm());
        if (supported && slopes[0]*slopes[1]>0 &&
            std::fabs(std::atan(slopes[0])-std::atan(slopes[1]))<20*3.141592653589793/180 && distance>tolerance) {
          next[active[k]]=true;changed=true;continue;
        }
      }
      std::size_t b=k,e=k+1;
      while (b>0 && q.x()-points[active[b-1]].x()<=radius) --b;
      while (e<active.size() && points[active[e]].x()-q.x()<=radius) ++e;
      if (e-b<9 || q.x()-points[active[b]].x()<.5f*radius ||
          points[active[e-1]].x()-q.x()<.5f*radius) continue;
      std::vector<double> weights(e-b,1);
      double intercept=0,slope=0;
      bool valid=true;
      for (int iteration=0;iteration<5;++iteration) {
        double sw=0,sx=0,sy=0,sxx=0,sxy=0;
        for (auto j=b;j<e;++j) {
          const double x=points[active[j]].x()-q.x(),y=points[active[j]].y()-q.y(),w=weights[j-b];
          sw+=w; sx+=w*x; sy+=w*y; sxx+=w*x*x; sxy+=w*x*y;
        }
        const double den=sw*sxx-sx*sx;
        if (den<1.e-9) { valid=false; break; }
        slope=(sw*sxy-sx*sy)/den; intercept=(sy-slope*sx)/sw;
        const double norm=std::hypot(1.0,slope);
        for (auto j=b;j<e;++j) {
          const double residual=std::fabs(points[active[j]].y()-q.y()-intercept-
                                        slope*(points[active[j]].x()-q.x()))/norm;
          weights[j-b]=std::min(1.0,.25*tolerance/std::max(1.e-9,residual));
        }
      }
      const double norm=std::hypot(1.0,slope);
      if (!valid || std::fabs(slope)/norm<.45 || std::fabs(intercept)/norm<=tolerance) continue;
      std::size_t supported=0; double flankSlope[2]={}; bool flanks=true;
      for (int side=0;side<2;++side) {
        double n=0,sx=0,sy=0,sxx=0,sxy=0,lo=radius,hi=-radius;
        for (auto j=b;j<e;++j) {
          const double x=points[active[j]].x()-q.x(),y=points[active[j]].y()-q.y();
          if (std::fabs(y-intercept-slope*x)/norm>.5*tolerance) continue;
          if (side==0) ++supported;
          if (side==0 ? x>=-tolerance : x<=tolerance) continue;
          ++n; sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; lo=std::min(lo,x);hi=std::max(hi,x);
        }
        if (n<3 || hi-lo<.25*radius || n*sxx-sx*sx<1.e-9) {flanks=false;break;}
        flankSlope[side]=(n*sxy-sx*sy)/(n*sxx-sx*sx);
      }
      if (!flanks || supported<.65*(e-b) || flankSlope[0]*flankSlope[1]<=0 ||
          std::fabs(std::atan(flankSlope[0])-std::atan(flankSlope[1]))>20*3.141592653589793/180) continue;
      next[active[k]]=true; changed=true;
    }
    removed=std::move(next);
    if (!changed) break;
  }
  // A long excursion is a feature requiring geometry analysis, not a spike.
  for (std::size_t b=0;b<points.size();) {
    if (!removed[b]) {++b;continue;}
    std::size_t e=b+1;
    while (e<points.size() && removed[e]) ++e;
    if (b==0 || e==points.size() || points[e].x()-points[b-1].x()>maximumWidth)
      for (auto j=b;j<e;++j) removed[j]=false;
    b=e;
  }
  return removed;
}

std::vector<ProfileCurvatureSample> fitParametricProfileSpline(
    const std::vector<Eigen::Vector2f> &points, float knotSpacing, float sampleStep) {
  std::vector<ProfileCurvatureSample> result;
  if (points.size()<7 || !std::isfinite(knotSpacing) || knotSpacing<=0 ||
      !std::isfinite(sampleStep) || sampleStep<=0) return result;
  std::vector<double> arc(points.size(),0);
  for (std::size_t i=0;i<points.size();++i) {
    if (!points[i].allFinite()) return result;
    if (i) arc[i]=arc[i-1]+(points[i]-points[i-1]).cast<double>().norm();
  }
  const double length=arc.back();
  if (length<1.e-5) return result;
  const int spans=std::max(1,std::min(int(points.size())-3,int(std::ceil(length/knotSpacing))));
  const int count=spans+3;
  std::vector<double> knots(count+4,length);
  for (int i=0;i<4;++i) knots[i]=0;
  for (int i=1;i<spans;++i) knots[i+3]=length*i/spans;
  auto basis=[&](double parameter) {
    const double s=std::max(0.0,std::min(length-1.e-9,parameter));
    Eigen::MatrixXd values=Eigen::MatrixXd::Zero(count,3);
    const int span=int(std::upper_bound(knots.begin(),knots.end(),s)-knots.begin())-1;
    if (span<3 || span>=count) return values;
    const int first=span-3;
    // A cubic span has only four nonzero basis functions. Retain the same
    // recurrence and dense result/reduction order, but omit all zero branches.
    // Two extra zero columns support the derivative expressions at the edge.
    double n[4][6] = {};
    n[0][3]=1;
    auto ratio=[](double a,double b){return b>1.e-12 ? a/b : 0.0;};
    for (int degree=1;degree<=3;++degree)
      for (int j=3-degree;j<4;++j) {
        const int i=first+j;
        n[degree][j]=ratio(s-knots[i],knots[i+degree]-knots[i])*n[degree-1][j]+
            ratio(knots[i+degree+1]-s,knots[i+degree+1]-knots[i+1])*n[degree-1][j+1];
      }
    for (int j=0;j<4;++j) {
      const int i=first+j;
      const double a=ratio(3,knots[i+3]-knots[i]), b=ratio(3,knots[i+4]-knots[i+1]);
      const double d0=ratio(2,knots[i+2]-knots[i])*n[1][j]-ratio(2,knots[i+3]-knots[i+1])*n[1][j+1];
      const double d1=ratio(2,knots[i+3]-knots[i+1])*n[1][j+1]-ratio(2,knots[i+4]-knots[i+2])*n[1][j+2];
      values(i,0)=n[3][j]; values(i,1)=a*n[2][j]-b*n[2][j+1]; values(i,2)=a*d0-b*d1;
    }
    return values;
  };
  std::vector<std::vector<std::pair<int,double>>> rows(points.size());
  for (std::size_t i=0;i<points.size();++i) {
    const auto b=basis(arc[i]);
    for (int j=0;j<count;++j) if (b(j,0)>0) rows[i].emplace_back(j,b(j,0));
  }
  const Eigen::Vector2d origin=points.front().cast<double>();
  std::vector<double> weights(points.size(),1);
  Eigen::MatrixXd controls;
  double rms=0;
  for (int pass=0;pass<3;++pass) {
    std::vector<Eigen::Triplet<double>> entries;
    Eigen::MatrixXd rhs=Eigen::MatrixXd::Zero(count,2);
    for (std::size_t i=0;i<points.size();++i) {
      // Equal arc-length influence rather than letting the dense floor dominate.
      const double support=i==0 ? arc[1]-arc[0] : i+1==points.size() ? arc[i]-arc[i-1] : .5*(arc[i+1]-arc[i-1]);
      const double w=weights[i]*std::max(.05,std::min(4.0,support/(length/(points.size()-1))));
      const Eigen::RowVector2d observed=(points[i].cast<double>()-origin).transpose();
      for (const auto &a:rows[i]) {
        rhs.row(a.first)+=w*a.second*observed;
        for (const auto &b:rows[i]) entries.emplace_back(a.first,b.first,w*a.second*b.second);
      }
    }
    const double d[3]={1,-2,1};
    for (int i=1;i+1<count;++i) for (int a=0;a<3;++a) for (int b=0;b<3;++b)
      entries.emplace_back(i+a-1,i+b-1,.15*d[a]*d[b]);
    for (int i=0;i<count;++i) entries.emplace_back(i,i,1.e-10);
    Eigen::SparseMatrix<double> normal(count,count);
    normal.setFromTriplets(entries.begin(),entries.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(normal);
    if (solver.info()!=Eigen::Success) return {};
    controls=solver.solve(rhs);
    if (solver.info()!=Eigen::Success || !controls.allFinite()) return {};
    rms=0;
    for (std::size_t i=0;i<points.size();++i) {
      Eigen::Vector2d fitted=origin;
      for (const auto &b:rows[i]) fitted+=b.second*controls.row(b.first).transpose();
      const double error=(fitted-points[i].cast<double>()).norm();
      rms+=error*error;
      weights[i]=std::min(1.0,std::max(.03,.1*double(knotSpacing))/std::max(1.e-9,error));
    }
    rms=std::sqrt(rms/points.size());
  }
  for (double s=0;s<=length;s+=sampleStep) {
    const auto b=basis(s);
    const Eigen::Vector2d point=origin+controls.transpose()*b.col(0);
    const Eigen::Vector2d velocity=controls.transpose()*b.col(1), acceleration=controls.transpose()*b.col(2);
    const double speed=velocity.norm();
    ProfileCurvatureSample sample;
    sample.point=point.cast<float>(); sample.rms=float(rms);
    if (speed>.1) {
      sample.tangent=(velocity/speed).cast<float>();
      sample.curvature=float((velocity.x()*acceleration.y()-velocity.y()*acceleration.x())/(speed*speed*speed));
      sample.valid=std::isfinite(sample.curvature);
    }
    result.push_back(sample);
  }
  return result;
}

std::vector<ProfileTurnRegion> locateProfileTurnRegions(
    const std::vector<ProfileCurvatureSample> &samples, bool angleWeighted,
    float relativeThreshold, float maxArcLength) {
  std::vector<ProfileTurnRegion> regions;
  const int n = int(samples.size());
  if (n < 3 || !std::isfinite(relativeThreshold) || relativeThreshold <= 0 ||
      relativeThreshold > 1 || !std::isfinite(maxArcLength) || maxArcLength <= 0) return regions;
  std::vector<double> arc(n, 0);
  auto valid = [&](int i) {
    return samples[i].valid && samples[i].point.allFinite() && std::isfinite(samples[i].curvature);
  };
  for (int i = 1; i < n; ++i)
    arc[i] = arc[i-1] + (valid(i) && valid(i-1) ? double((samples[i].point-samples[i-1].point).norm()) : 0.0);
  std::vector<int> peaks;
  for (int i = 0; i < n;) {
    int end = i;
    while (end+1 < n && valid(end) && valid(end+1) &&
           samples[end+1].curvature == samples[i].curvature) ++end;
    if (i > 0 && end+1 < n && valid(i-1) && valid(i) && valid(end+1) &&
        std::fabs(samples[i].curvature) > 1e-6f &&
        std::fabs(samples[i].curvature) > std::fabs(samples[i-1].curvature) &&
        std::fabs(samples[i].curvature) > std::fabs(samples[end+1].curvature)) {
      const double center = .5*(arc[i]+arc[end]);
      int peak = i;
      while (peak < end && std::fabs(arc[peak+1]-center) < std::fabs(arc[peak]-center)) ++peak;
      peaks.push_back(peak);
    }
    i = end+1;
  }
  std::stable_sort(peaks.begin(), peaks.end(), [&](int a, int b) {
    return std::fabs(samples[a].curvature) > std::fabs(samples[b].curvature);
  });
  std::vector<bool> claimed(n, false);
  for (int peak : peaks) {
    if (claimed[peak]) continue;
    const float height = std::fabs(samples[peak].curvature);
    auto eligible = [&](int i) {
      return valid(i) && !claimed[i] && samples[i].curvature*samples[peak].curvature > 0 &&
          std::fabs(samples[i].curvature) >= relativeThreshold*height;
    };
    int begin = peak, end = peak;
    while (begin > 0 && eligible(begin-1) && arc[peak]-arc[begin-1] <= .5*maxArcLength) --begin;
    while (end+1 < n && eligible(end+1) && arc[end+1]-arc[peak] <= .5*maxArcLength) ++end;
    std::vector<double> mass(end-begin+1, 0);
    for (int i = begin+1; i <= end; ++i)
      mass[i-begin] = mass[i-begin-1] + .5*(std::fabs(samples[i-1].curvature)+std::fabs(samples[i].curvature))*(arc[i]-arc[i-1]);
    const double total = angleWeighted ? mass.back() : arc[end]-arc[begin];
    int representative = begin;
    double best = std::numeric_limits<double>::infinity();
    for (int i = begin; i <= end; ++i) {
      const double distance = std::fabs((angleWeighted ? mass[i-begin] : arc[i]-arc[begin])-.5*total);
      if (distance < best) { best = distance; representative = i; }
      claimed[i] = true;
    }
    ProfileTurnRegion region;
    region.begin = begin; region.end = end; region.peak = peak; region.representative = representative;
    region.arcLength = float(arc[end]-arc[begin]); region.turningAngle = float(mass.back());
    regions.push_back(region);
  }
  std::sort(regions.begin(), regions.end(), [](const ProfileTurnRegion &a, const ProfileTurnRegion &b) {
    return a.representative < b.representative;
  });
  return regions;
}

int profileTurnFloorBoundary(
    const std::vector<ProfileCurvatureSample> &samples,
    const ProfileTurnRegion &region, int openingDirection) {
  if (samples.empty() || region.begin < 0 || region.end < region.begin ||
      region.end >= int(samples.size()))
    return region.representative;
  int opening = openingDirection;
  if (opening == 0) {
    const float edgeV = .5f * (samples.front().point.y() +
                               samples.back().point.y());
    opening = samples[samples.size() / 2].point.y() <= edgeV ? 1 : -1;
  }
  const float beginDepth = opening * samples[region.begin].point.y();
  const float endDepth = opening * samples[region.end].point.y();
  if (std::fabs(beginDepth - endDepth) > 1e-6f)
    return beginDepth < endDepth ? region.begin : region.end;

  // A locally level transition can have equal endpoint V. Its inner endpoint
  // is the one closer to the transverse centre of the complete profile.
  const float centerU = .5f * (samples.front().point.x() +
                               samples.back().point.x());
  return std::fabs(samples[region.begin].point.x() - centerU) <=
                 std::fabs(samples[region.end].point.x() - centerU)
             ? region.begin
             : region.end;
}

std::vector<ProfileCurvatureSample> fitProfileCurvature(
    const std::vector<Eigen::Vector2f> &points, float radius) {
  std::vector<ProfileCurvatureSample> result(points.size());
  if (points.size() < 7 || !std::isfinite(radius) || radius <= 0) return result;
  std::vector<double> arc(points.size(), 0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!points[i].allFinite()) return result;
    if (i) arc[i] = arc[i-1] + (points[i]-points[i-1]).cast<double>().norm();
  }
  for (std::size_t i = 3; i + 3 < points.size(); ++i) {
    const auto first = std::lower_bound(arc.begin(), arc.end(), arc[i]-radius);
    const auto last = std::upper_bound(arc.begin(), arc.end(), arc[i]+radius);
    const int begin = int(first-arc.begin()), end = int(last-arc.begin());
    if (int(i)-begin < 3 || end-int(i)-1 < 3) continue;
    const double length = arc[end-1]-arc[begin];
    if (length < 1.e-6) continue;
    const int n = end-begin;
    Eigen::MatrixXd basis(n,4), observed(n,2);
    Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
    for (int j = 0; j < n; ++j) {
      const double t = (arc[begin+j]-arc[begin])/length, q = 1-t;
      basis.row(j) << q*q*q, 3*t*q*q, 3*t*t*q, t*t*t;
      observed.row(j) = (points[begin+j]-points[i]).cast<double>().transpose();
    }
    Eigen::Matrix<double,4,2> controls;
    bool valid = true;
    double rms = 0;
    for (int pass = 0; pass < 3; ++pass) {
      const Eigen::Matrix4d normal = basis.transpose()*weights.asDiagonal()*basis;
      Eigen::LDLT<Eigen::Matrix4d> solver(normal);
      if (solver.info() != Eigen::Success || solver.vectorD().minCoeff() < 1.e-10) {
        valid = false; break;
      }
      controls = solver.solve(basis.transpose()*weights.asDiagonal()*observed);
      if (!controls.allFinite()) { valid = false; break; }
      const Eigen::MatrixXd residual = basis*controls-observed;
      rms = std::sqrt(residual.squaredNorm()/n);
      // Bounded influence for isolated scan spikes; do not force interpolation.
      const double scale = std::max(.02*double(radius), rms);
      for (int j = 0; j < n; ++j)
        weights[j] = std::min(1.0, scale/std::max(1.e-9,residual.row(j).norm()));
    }
    if (!valid) continue;
    const double t = (arc[i]-arc[begin])/length, q = 1-t;
    Eigen::RowVector4d b, d, dd;
    b << q*q*q, 3*t*q*q, 3*t*t*q, t*t*t;
    d << -3*q*q, 3*q*q-6*t*q, 6*t*q-3*t*t, 3*t*t;
    dd << 6*q, -12+18*t, 6-18*t, 6*t;
    const Eigen::RowVector2d velocity = d*controls/length;
    const Eigen::RowVector2d acceleration = dd*controls/(length*length);
    const double speed = velocity.norm();
    if (speed < .1) continue;
    auto &sample = result[i];
    sample.point = points[i]+(b*controls).transpose().cast<float>();
    sample.tangent = (velocity/speed).transpose().cast<float>();
    sample.curvature = float((velocity[0]*acceleration[1]-velocity[1]*acceleration[0])/
                             (speed*speed*speed));
    sample.rms = float(rms);
    sample.valid = std::isfinite(sample.curvature) && sample.point.allFinite();
  }
  return result;
}
} // namespace weld

namespace weld {
std::vector<ProfileTurnSample> measureProfileTurns(
    const std::vector<Eigen::Vector2f> &points, float span, float step,
    const std::vector<bool> *eligibleIndices) {
  std::vector<ProfileTurnSample> turns;
  if (points.size()<3 || !std::isfinite(span) || span<=0 ||
      !std::isfinite(step) || step<=0) return turns;
  if (eligibleIndices && eligibleIndices->size()!=points.size()) return turns;
  std::vector<double> arc(points.size(),0);
  for (std::size_t i=0;i<points.size();++i) {
    if (!points[i].allFinite()) return turns;
    if (i) arc[i]=arc[i-1]+(points[i]-points[i-1]).cast<double>().norm();
  }
  auto at = [&](double s) -> Eigen::Vector2d {
    const auto next = std::upper_bound(arc.begin(),arc.end(),s);
    if (next==arc.end()) return points.back().cast<double>();
    const auto j=std::size_t(next-arc.begin());
    if (!j) return points.front().cast<double>();
    const double t=(s-arc[j-1])/std::max(1.e-9,arc[j]-arc[j-1]);
    return (1-t)*points[j-1].cast<double>()+t*points[j].cast<double>();
  };
  auto tangent = [&](double start, double end, Eigen::Vector2f &direction,
                     Eigen::Vector2f &lineCenter, float &rms) {
    Eigen::Matrix<double,7,2> samples;
    const Eigen::Vector2d origin=at(start);
    for (int j=0;j<7;++j) samples.row(j)=(at(start+(end-start)*j/6)-origin).transpose();
    const Eigen::RowVector2d center=samples.colwise().mean();
    lineCenter=(center.transpose()+origin).cast<float>();
    samples.rowwise()-=center;
    const Eigen::Matrix2d covariance=samples.transpose()*samples/7;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(covariance);
    if (eig.info()!=Eigen::Success || eig.eigenvalues()[1]<1.e-8) return false;
    direction=eig.eigenvectors().col(1).cast<float>();
    if (direction.cast<double>().dot(at(end)-origin)<0) direction=-direction;
    rms=float(std::sqrt(std::max(0.0,eig.eigenvalues()[0])));
    return true;
  };
  // Include original vertices as well as sub-bin probes; a sharp corner need
  // not coincide with an arbitrary resampling origin.
  std::vector<double> probes;
  for (double s=span; s+span<=arc.back(); s+=step) probes.push_back(s);
  for (double s:arc) if (s>=span && s+span<=arc.back()) probes.push_back(s);
  std::sort(probes.begin(),probes.end());
  probes.erase(std::unique(probes.begin(),probes.end(),[](double a,double b){return std::abs(a-b)<1.e-6;}),probes.end());
  for (double s:probes) {
    ProfileTurnSample turn;
    auto next=std::lower_bound(arc.begin(),arc.end(),s);
    turn.index=int(next-arc.begin());
    if (turn.index>0 && s-arc[turn.index-1]<arc[turn.index]-s) --turn.index;
    // Preserve the original probe positions and nearest-index rule exactly;
    // only calculate expensive one-sided PCA tangents at eligible peaks.
    if (eligibleIndices && !(*eligibleIndices)[turn.index]) continue;
    if (!tangent(s-span,s,turn.before,turn.beforeCenter,turn.beforeRms) ||
        !tangent(s,s+span,turn.after,turn.afterCenter,turn.afterRms)) continue;
    turn.point=at(s).cast<float>();
    turn.angle=std::atan2(turn.before.x()*turn.after.y()-turn.before.y()*turn.after.x(),
                          turn.before.dot(turn.after));
    turn.curvature=turn.angle/span;
    turn.rms=std::max(turn.beforeRms,turn.afterRms);
    turns.push_back(turn);
  }
  return turns;
}
} // namespace weld
