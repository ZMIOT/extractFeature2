#include "ProfileCurvature.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}

void robustSplineRegression() {
  // A repeated point exercises zero-length support; an isolated spike changes
  // the robust weights between passes. Reference values are from main 880a55a.
  std::vector<Eigen::Vector2f> points;
  for(int i=0;i<17;++i) {
    const float x=.25f*i;
    points.emplace_back(x,.1f*x*x+(i%2 ? .02f : -.02f)+(i==8 ? 1.5f : 0.f));
    if(i==4)points.back()=points[3];
  }
  const auto fitted=weld::fitParametricProfileSpline(points,.5f,.1f);
  require(fitted.size()==69,"Robust spline sample count changed");
  struct Reference {int index;float x,y,tx,ty,curvature;};
  const Reference expected[] = {
    {0,-.02776932716f,-.01823799126f,.9936991334f,.1120806932f,-.08412594348f},
    {10,.9910771847f,.07649752498f,.9920762181f,.1256375909f,.207397759f},
    {20,1.840291381f,.5702598691f,.4450773597f,.8954921365f,.9843879938f},
    {30,2.019850731f,1.714074612f,-.04906757548f,.9987954497f,.4847694039f},
    {68,3.998784542f,1.584982395f,.8329480886f,.5533511043f,-.1360488236f}
  };
  for(const auto &r:expected) {
    const auto &s=fitted[r.index];
    require(s.valid,"Robust spline sample became invalid");
    require((s.point-Eigen::Vector2f(r.x,r.y)).norm()<2.e-6f &&
            (s.tangent-Eigen::Vector2f(r.tx,r.ty)).norm()<2.e-6f &&
            std::fabs(s.curvature-r.curvature)<2.e-6f &&
            std::fabs(s.rms-.05729499459f)<2.e-6f,
            "Robust spline geometry or residual changed");
  }
}

void degenerateInputs() {
  std::vector<Eigen::Vector2f> points(12,Eigen::Vector2f::Zero());
  require(weld::fitParametricProfileSpline(points,1,.1f).empty(),"Zero-length profile was fitted");
  for(std::size_t i=0;i<points.size();++i)points[i].x()=float(i);
  require(weld::fitParametricProfileSpline(points,0,.1f).empty(),"Zero knot spacing was accepted");
  require(weld::fitParametricProfileSpline(points,1,0).empty(),"Zero sample step was accepted");
  points[5].y()=std::numeric_limits<float>::quiet_NaN();
  require(weld::fitParametricProfileSpline(points,1,.1f).empty(),"Nonfinite profile was fitted");
  points.resize(6);
  require(weld::fitParametricProfileSpline(points,1,.1f).empty(),"Insufficient support was fitted");
}
}

int main() {
  try {
    robustSplineRegression();
    degenerateInputs();
    std::cout<<"Robust profile spline regression passed\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
