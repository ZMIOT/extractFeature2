#include <WeldSeam.h>
#include <iostream>
#include <cmath>

int main() {
  std::string error;
  // In a caller's application, replace this synthetic groove with measured XYZ.
  // Use the library's unmodified defaults; no configuration file is required.
  weld::CloudT::Ptr input(new weld::CloudT);
  for (int x=0; x<80; ++x) for (int j=-250; j<=250; ++j) {
    const float y = .1f*j;
    const float z = std::min(15.f, std::max(0.f, 2*(std::fabs(y)-8.f)));
    input->push_back(weld::PointT(.8f*x, y, z));
  }
  weld::SeamClouds output;
  if (!weld::extractSeams(input, output, weld::defaultDetectorConfig(), &error)) {
    std::cerr << error << '\n'; return 3;
  }
  // Measured clouds and fitted surface functions; storage/display belong to caller.
  std::cout << "left=" << output.left->size() << " right=" << output.right->size()
            << " center=" << output.center->size()
            << " left_surface=" << output.leftSurface->size()
            << " right_surface=" << output.rightSurface->size()
            << " weld_surface=" << output.weldSurface->size() << '\n';
  if (output.leftSurface->empty() || output.rightSurface->empty() || output.weldSurface->empty()) {
    std::cerr << "Insufficient measured support for all three surfaces\n";
    return 4;
  }
  for(const auto *model:{&output.leftSurfaceModel,&output.weldSurfaceModel,&output.rightSurfaceModel}) {
    if(!model->valid()) {std::cerr<<"Surface fitting unavailable: "<<model->message<<'\n';return 5;}
    Eigen::Vector3d point,normal;
    if(!model->evaluate(0,.5,.5,point,&normal))return 6;
    std::cout<<"patches="<<model->patches.size()<<" S(0.5,0.5)="<<point.transpose()
             <<" normal="<<normal.transpose()<<" RMS="<<model->patches[0].rmsError<<'\n';
  }
  return 0;
}
