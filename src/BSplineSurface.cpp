#include "SurfaceFitting.h"
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <functional>

namespace weld {
namespace {
std::vector<double> basis(const std::vector<double> &knots,int degree,double t) {
  const int count=int(knots.size())-degree-1;
  std::vector<double> values(knots.size()-1,0.0);
  // Use the final nonempty span at the closed right endpoint.
  int span=int(knots.size())-2;
  if(t==knots.back()) {
    while(span>0 && knots[span]==knots.back())--span;
    values[span]=1;
  } else {
    for(int i=0;i<int(values.size());++i)
      if(t>=knots[i] && t<knots[i+1])values[i]=1;
  }
  for(int p=1;p<=degree;++p)for(int i=0;i<int(knots.size())-p-1;++i) {
    const double a=knots[i+p]-knots[i], b=knots[i+p+1]-knots[i+1];
    values[i]=(a>0?(t-knots[i])*values[i]/a:0)+
              (b>0?(knots[i+p+1]-t)*values[i+1]/b:0);
  }
  values.resize(count);return values;
}
std::vector<double> derivative(const std::vector<double> &knots,int degree,double t) {
  const auto lower=basis(knots,degree-1,t);
  std::vector<double> result(lower.size()-1);
  for(int i=0;i<int(result.size());++i) {
    const double a=knots[i+degree]-knots[i],b=knots[i+degree+1]-knots[i+1];
    result[i]=degree*((a>0?lower[i]/a:0)-(b>0?lower[i+1]/b:0));
  }
  return result;
}
bool validKnots(const std::vector<double> &k,int count,int degree) {
  if(degree!=3 || count<4 || count>4096 || k.size()!=std::size_t(count+degree+1))return false;
  for(double x:k)if(!std::isfinite(x))return false;
  if(!std::is_sorted(k.begin(),k.end()) || k.front()!=0 || k.back()!=1)return false;
  for(int i=0;i<=degree;++i)if(k[i]!=0 || k[k.size()-1-i]!=1)return false;
  for(int i=degree+1;i<count;++i)if(k[i]<=k[i-1] || k[i]>=1)return false;
  return true;
}
std::vector<double> makeKnots(int count) {
  std::vector<double> k(count+4,0);
  for(int i=4;i<count;++i)k[i]=double(i-3)/(count-3);
  std::fill(k.begin()+count,k.end(),1.);return k;
}
Eigen::Vector3d value(const BSplineSurfacePatch &patch,double u,double v) {
  const auto bu=basis(patch.knotsU,3,u),bv=basis(patch.knotsV,3,v);
  Eigen::Vector3d point=Eigen::Vector3d::Zero();
  for(int i=0;i<patch.countU;++i)if(bu[i]!=0)
    for(int j=0;j<patch.countV;++j)if(bv[j]!=0) {
      const auto &p=patch.controlPoints[i*patch.countV+j];
      point+=bu[i]*bv[j]*Eigen::Vector3d(p[0],p[1],p[2]);
    }
  return point;
}
struct Row {
  std::size_t section=0;
  std::vector<Eigen::Vector3d> points;
  std::vector<double> arc;
  Eigen::Vector3d center=Eigen::Vector3d::Zero();
  double u=0;
};
Eigen::Vector3d interpolate(const Row &row,double t) {
  const auto it=std::lower_bound(row.arc.begin(),row.arc.end(),t);
  if(it==row.arc.begin())return row.points.front();
  if(it==row.arc.end())return row.points.back();
  const auto j=std::size_t(it-row.arc.begin());
  const double weight=(t-row.arc[j-1])/(row.arc[j]-row.arc[j-1]);
  return (1-weight)*row.points[j-1]+weight*row.points[j];
}
bool fitRun(std::vector<Row> &rows,const std::vector<SurfaceSection> &sections,
            const DetectorConfig &config,BSplineSurfacePatch &out) {
  if(rows.size()<4)return false;
  // Projection chooses a consistent contour direction even for vertical walls.
  const int a=int(config.longitudinalAxis),b=a==0?1:0,c=a==2?1:2;
  Eigen::Matrix2d covariance=Eigen::Matrix2d::Zero();
  for(auto &row:rows) {
    row.center.setZero();row.u=0;
    for(const auto &p:row.points)row.center+=p;
    row.center/=double(row.points.size());
    for(const auto &p:row.points) {
      const Eigen::Vector2d d(p[b]-row.center[b],p[c]-row.center[c]);
      covariance.noalias()+=d*d.transpose()/double(row.points.size());
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(covariance);
  if(eigen.info()!=Eigen::Success || eigen.eigenvalues()[1]<1e-12)return false;
  Eigen::Vector2d direction=eigen.eigenvectors().col(1);
  if(direction[std::abs(direction[0])>=std::abs(direction[1])?0:1]<0)direction=-direction;
  double maximumWidth=0;
  for(std::size_t r=0;r<rows.size();++r) {
    auto &row=rows[r];
    auto projection=[&](const Eigen::Vector3d &p){return direction[0]*p[b]+direction[1]*p[c];};
    std::sort(row.points.begin(),row.points.end(),[&](const Eigen::Vector3d &p,const Eigen::Vector3d &q){
      const double x=projection(p),y=projection(q);
      if(x!=y)return x<y;
      for(int k=0;k<3;++k)if(p[k]!=q[k])return p[k]<q[k];
      return false;
    });
    row.points.erase(std::unique(row.points.begin(),row.points.end(),
        [](const Eigen::Vector3d &p,const Eigen::Vector3d &q){return (p-q).squaredNorm()<1e-16;}),row.points.end());
    if(row.points.size()<4)return false;
    // Chord-length parameterization in the section plane. Original XYZ is fitted.
    row.arc.assign(row.points.size(),0);
    for(std::size_t j=1;j<row.points.size();++j) {
      const auto d=row.points[j]-row.points[j-1];
      row.arc[j]=row.arc[j-1]+std::hypot(d[b],d[c]);
    }
    if(row.arc.back()<1e-6)return false;
    maximumWidth=std::max(maximumWidth,row.arc.back());
    const double length=row.arc.back();for(auto &t:row.arc)t/=length;
    if(r) {
      // Transverse variation of a cropped wall must not invent longitudinal support.
      const double step=(row.center-rows[r-1].center).norm();
      if(step<1e-6)return false;
      row.u=rows[r-1].u+step;
    }
  }
  const double length=rows.back().u;
  if(length<1e-6)return false;
  for(auto &row:rows)row.u/=length;
  // A single tilted scan (or repeated scans) cannot define a surface.
  Eigen::Vector3d mean=Eigen::Vector3d::Zero();std::size_t n=0;
  for(const auto &row:rows)for(const auto &p:row.points){mean+=p;++n;}
  mean/=double(n);
  Eigen::Matrix3d cov=Eigen::Matrix3d::Zero();
  for(const auto &row:rows)for(const auto &p:row.points){const Eigen::Vector3d d=p-mean;cov.noalias()+=d*d.transpose();}
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> rank(cov);
  if(rank.info()!=Eigen::Success || rank.eigenvalues()[1]<1e-10*rank.eigenvalues()[2])return false;
  const double spacing=std::max(1.,8.*config.profileBinSize);
  out.countU=std::min(int(rows.size()),std::max(4,std::min(32,int(std::ceil(length/spacing))+3)));
  out.countV=std::max(4,std::min(12,int(std::ceil(maximumWidth/spacing))+3));
  out.knotsU=makeKnots(out.countU);out.knotsV=makeKnots(out.countV);
  struct Sample {std::vector<std::pair<int,double>> weights;Eigen::Vector3d p;double robust=1;};
  std::vector<Sample> samples;
  const int rowCount=std::min(256,int(rows.size())),columnCount=48,controlCount=out.countU*out.countV;
  for(int r=0;r<rowCount;++r) {
    const auto &row=rows[std::size_t(r)*(rows.size()-1)/(rowCount-1)];
    const auto bu=basis(out.knotsU,3,row.u);
    for(int k=0;k<columnCount;++k) {
      const double v=double(k)/(columnCount-1);const auto bv=basis(out.knotsV,3,v);
      Sample s;s.p=interpolate(row,v)-mean;
      for(int i=0;i<out.countU;++i)if(bu[i]>0)
        for(int j=0;j<out.countV;++j)if(bv[j]>0)s.weights.emplace_back(i*out.countV+j,bu[i]*bv[j]);
      samples.push_back(std::move(s));
    }
  }
  Eigen::MatrixXd controls;
  const double robustScale=std::max(.01,double(std::max(config.profileWallTolerance,config.hullSurfaceDistance)));
  const double smoothness=1e-4*samples.size()/controlCount;
  for(int iteration=0;iteration<3;++iteration) {
    std::vector<Eigen::Triplet<double>> terms;terms.reserve(samples.size()*256+controlCount*20);
    Eigen::MatrixXd rhs=Eigen::MatrixXd::Zero(controlCount,3);
    for(const auto &s:samples)for(const auto &x:s.weights) {
      rhs.row(x.first)+=s.robust*x.second*s.p.transpose();
      for(const auto &y:s.weights)terms.emplace_back(x.first,y.first,s.robust*x.second*y.second);
    }
    auto smooth=[&](int i,int j,int k) {
      const int ids[]={i,j,k};const double w[]={1,-2,1};
      for(int x=0;x<3;++x)for(int y=0;y<3;++y)terms.emplace_back(ids[x],ids[y],smoothness*w[x]*w[y]);
    };
    for(int i=0;i<out.countU;++i)for(int j=0;j<out.countV;++j) {
      const int index=i*out.countV+j;
      if(i+2<out.countU)smooth(index,index+out.countV,index+2*out.countV);
      if(j+2<out.countV)smooth(index,index+1,index+2);
      terms.emplace_back(index,index,1e-10);
    }
    Eigen::SparseMatrix<double> matrix(controlCount,controlCount);matrix.setFromTriplets(terms.begin(),terms.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;solver.compute(matrix);
    if(solver.info()!=Eigen::Success)return false;
    controls=solver.solve(rhs);
    if(solver.info()!=Eigen::Success || !controls.allFinite())return false;
    for(auto &s:samples) {
      Eigen::Vector3d fitted=Eigen::Vector3d::Zero();
      for(const auto &x:s.weights)fitted+=x.second*controls.row(x.first).transpose();
      s.robust=std::min(1.,robustScale/std::max(1e-12,(fitted-s.p).norm()));
    }
  }
  out.controlPoints.resize(controlCount);
  for(int i=0;i<controlCount;++i)for(int j=0;j<3;++j)out.controlPoints[i][j]=controls(i,j)+mean[j];
  double squared=0;
  for(const auto &row:rows)for(std::size_t j=0;j<row.points.size();++j) {
    const double error=(value(out,row.u,row.arc[j])-row.points[j]).norm();
    squared+=error*error;out.maxError=std::max(out.maxError,error);++out.sampleCount;
  }
  out.rmsError=std::sqrt(squared/out.sampleCount);
  out.firstFrame=sections[rows.front().section].frameIndex;out.lastFrame=sections[rows.back().section].frameIndex;
  out.firstStation=sections[rows.front().section].station;out.lastStation=sections[rows.back().section].station;
  return out.valid() && out.rmsError<=2*robustScale;
}
}
bool BSplineSurfacePatch::valid() const {
  if(!validKnots(knotsU,countU,degreeU) || !validKnots(knotsV,countV,degreeV) ||
     controlPoints.size()!=std::size_t(countU)*countV)return false;
  for(const auto &p:controlPoints)for(double x:p)if(!std::isfinite(x))return false;
  return true;
}
bool BSplineSurfacePatch::evaluate(double u,double v,Eigen::Vector3d &point,Eigen::Vector3d *normal) const {
  if(!valid() || !std::isfinite(u) || !std::isfinite(v) || u<0 || u>1 || v<0 || v>1)return false;
  const auto p=value(*this,u,v);
  Eigen::Vector3d n=Eigen::Vector3d::Zero();
  if(normal) {
    const auto bu=basis(knotsU,3,u),bv=basis(knotsV,3,v);
    const auto du=derivative(knotsU,3,u),dv=derivative(knotsV,3,v);
    Eigen::Vector3d tu=Eigen::Vector3d::Zero(),tv=Eigen::Vector3d::Zero();
    // Subtract an origin to avoid cancellation with large world coordinates.
    const auto &origin=controlPoints.front();
    for(int i=0;i<countU;++i)for(int j=0;j<countV;++j) {
      const auto &cp=controlPoints[i*countV+j];
      const Eigen::Vector3d d(cp[0]-origin[0],cp[1]-origin[1],cp[2]-origin[2]);
      tu+=du[i]*bv[j]*d;tv+=bu[i]*dv[j]*d;
    }
    n=tu.cross(tv);
    if(!n.allFinite() || n.norm()<=1e-12)return false;
    n.normalize();
  }
  if(!p.allFinite())return false;
  point=p;if(normal)*normal=n;return true;
}
bool BSplineSurfaceModel::evaluate(std::size_t patch,double u,double v,Eigen::Vector3d &point,Eigen::Vector3d *normal) const {
  return patch<patches.size() && patches[patch].evaluate(u,v,point,normal);
}
BSplineSurfaceModel fitSurfaceModel(const CloudT &source,
    const std::vector<std::vector<std::size_t>> &groups,const std::vector<SurfaceSection> &sections,
    const std::vector<unsigned char> &labels,unsigned char face,const DetectorConfig &config) {
  BSplineSurfaceModel model;std::vector<Row> rows;
  std::function<void(std::vector<Row>)> fitPart;
  fitPart=[&](std::vector<Row> part) {
    BSplineSurfacePatch patch;
    if(fitRun(part,sections,config,patch))model.patches.push_back(std::move(patch));
    else if(patch.valid() && part.size()>=7) {
      // A rapidly changing boundary can defeat one long chart. Subdivide at
      // a measured section, sharing that section but imposing no artificial
      // continuity across independently fitted patches.
      const auto middle=part.size()/2;
      fitPart(std::vector<Row>(part.begin(),part.begin()+middle+1));
      fitPart(std::vector<Row>(part.begin()+middle,part.end()));
    } else ++model.rejectedPatchCount;
  };
  auto flush=[&] {
    if(rows.empty())return;
    fitPart(std::move(rows));
    rows.clear();
  };
  for(std::size_t i=0;i<sections.size();++i) {
    Row row;row.section=i;
    for(auto index:groups[i])if(labels[index]==face)row.points.push_back(toVector(source[index]).cast<double>());
    if(!rows.empty()) {
      const auto &previous=sections[rows.back().section],&current=sections[i];
      const bool gap=current.frameIndex>=0 ? current.frameIndex!=previous.frameIndex+1 :
          current.station-previous.station>1.5*config.sliceStep;
      if(gap)flush();
    }
    if(row.points.size()<4){flush();continue;}
    rows.push_back(std::move(row));
  }
  flush();
  model.message=model.valid() ? (model.rejectedPatchCount?"partial_support":"ok") : "insufficient_support_or_fit_error";
  return model;
}
}
