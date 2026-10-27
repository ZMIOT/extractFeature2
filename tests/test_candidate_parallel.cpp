#include "CandidateExecutor.h"
#include "SeamDetector.h"
#include "AlgorithmParameters.h"
#include <cmath>
#include <future>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool ok, const char *message) { if (!ok) throw std::runtime_error(message); }
void threads(const char *value) {
#ifdef _WIN32
  _putenv_s("WELD_SEAM_CANDIDATE_THREADS", value);
#else
  if (*value) setenv("WELD_SEAM_CANDIDATE_THREADS", value, 1);
  else unsetenv("WELD_SEAM_CANDIDATE_THREADS");
#endif
}
struct RestoreEnvironment {
  std::string previous;
  RestoreEnvironment() { if (const auto *p=std::getenv("WELD_SEAM_CANDIDATE_THREADS")) previous=p; }
  ~RestoreEnvironment() { threads(previous.c_str()); }
};
template<class T> void value(std::ostream &out, const T &v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}
template<class T> void values(std::ostream &out, const std::vector<T> &v) {
  value(out,v.size()); for(const auto &x:v) value(out,x);
}
void model(std::ostream &out, const weld::BSplineSurfaceModel &m) {
  value(out,m.patches.size());value(out,m.rejectedPatchCount);out<<m.message<<'\0';
  for(const auto &p:m.patches) {
    value(out,p.degreeU);value(out,p.degreeV);value(out,p.countU);value(out,p.countV);
    values(out,p.knotsU);values(out,p.knotsV);values(out,p.controlPoints);
    value(out,p.firstFrame);value(out,p.lastFrame);value(out,p.firstStation);value(out,p.lastStation);
    value(out,p.sampleCount);value(out,p.rmsError);value(out,p.maxError);
  }
}
std::string snapshot(const weld::SeamResult &r) {
  std::ostringstream out(std::ios::out|std::ios::binary);
  for(const auto &c:{r.leftSeam,r.rightSeam,r.centerSeam,r.leftSurface,r.weldSurface,r.rightSurface,r.mainBodyCloud,r.filteredCloud}) {
    value(out,c->size());for(const auto &p:*c){value(out,p.x);value(out,p.y);value(out,p.z);}
  }
  for(const auto *m:{&r.leftSurfaceModel,&r.weldSurfaceModel,&r.rightSurfaceModel}) model(out,*m);
  for(const auto *v:{&r.leftFrameIndices,&r.rightFrameIndices,&r.centerFrameIndices})values(out,*v);
  for(const auto *v:{&r.leftSegmentStarts,&r.rightSegmentStarts,&r.centerSegmentStarts})values(out,*v);
  for(const auto *v:{&r.leftPointConfidence,&r.rightPointConfidence,&r.centerPointConfidence})values(out,*v);
  for(const auto *v:{&r.leftPointQuality,&r.rightPointQuality,&r.centerPointQuality}) {
    value(out,v->size());for(const auto &p:*v) {
      value(out,p.sourceConfidence);value(out,p.wallRms);value(out,p.floorRms);value(out,p.supportDistance);
      value(out,p.trackedWall);value(out,p.trackedOffset);
      for(int i=0;i<2;++i)value(out,p.trackedNormal[i]);
      value(out,p.geometryCorrection);value(out,p.rtsCorrection);value(out,p.finalCurveCorrection);
      value(out,p.innovationNis);value(out,p.positionSigma);value(out,p.flags);
      value(out,p.wallAccepted);value(out,p.floorAccepted);value(out,p.wallOffset);value(out,p.floorOffset);
      for(int i=0;i<2;++i)value(out,p.wallNormal[i]);for(int i=0;i<3;++i)value(out,p.floorNormal[i]);
    }
  }
  value(out,r.confidence);value(out,r.leftConfidence);value(out,r.rightConfidence);
  value(out,r.leftHistoryUsageRatio);value(out,r.rightHistoryUsageRatio);
  value(out,r.endpointTrimmedStart);value(out,r.endpointTrimmedEnd);value(out,r.candidateCorrectionRejected);
  value(out,r.inputSamplingPoints);value(out,r.inputSamplingRemovedPoints);value(out,r.removedBodyPoints);
  out<<r.summary;
  return out.str();
}
void poolTests() {
  weld::detail::CandidateExecutor pool(4);
  const auto caller=std::this_thread::get_id();
  pool.run(0,[&](std::size_t){throw std::runtime_error("Unexpected empty work");});
  pool.run(3,[&](std::size_t){require(std::this_thread::get_id()==caller,"Small batch dispatched workers");});
  std::vector<std::atomic<int>> hits(4096);
  for(int repeat=0;repeat<20;++repeat) {
    for(auto &h:hits)h.store(0);
    pool.run(hits.size(),[&](std::size_t i){++hits[i];});
    for(const auto &h:hits)require(h.load()==1,"Pool lost or duplicated a task");
  }
  std::atomic<bool> workerSeen{false};bool caught=false;
  try {
    pool.run(128,[&](std::size_t) {
      if(std::this_thread::get_id()!=caller) {workerSeen=true;throw std::runtime_error("worker failure");}
      const auto start=std::chrono::steady_clock::now();
      while(!workerSeen && std::chrono::steady_clock::now()-start<std::chrono::seconds(5))std::this_thread::yield();
      require(workerSeen,"Workers did not start");
    });
  } catch(const std::runtime_error &e) {caught=std::string(e.what())=="worker failure";}
  require(caught,"Worker exception did not reach the caller");
  for(auto &h:hits)h.store(0);
  pool.run(hits.size(),[&](std::size_t i){++hits[i];});
  for(const auto &h:hits)require(h.load()==1,"Pool was not reusable after exception");
}
void detectionTests(bool bead,weld::ProfileLocalizationMethod method) {
  weld::CloudT::Ptr input(new weld::CloudT);std::vector<weld::ScanFrameRange> frames;
  for(int f=0;f<24;++f) {
    if(f==12)continue;
    const auto begin=input->size();
    for(int j=-300;j<=300;++j) {
      const float y=.1f*j,z=std::min(24.f,std::max(0.f,2*(std::fabs(y)-12)))+
          (bead?2.6f*std::exp(-.5f*std::pow(y/2,2)):0);
      input->push_back({.6f*f+.25f*y+.1f*z,y,z});
    }
    frames.push_back({begin,input->size()-begin,100+f});
  }
  auto config=weld::AlgorithmParameters().get(method);
  config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
  config.trimUnstableEndpoints=false;
  const weld::SeamDetector detector;
  threads("1");
  const auto reference=detector.detect(input,config,nullptr,nullptr,nullptr,&frames,true);
  require(reference.ok(),"Parallel test fixture failed detection");
  const auto expected=snapshot(reference);
  for(const auto *count:{"2","4","8"}) {
    threads(count);
    require(snapshot(detector.detect(input,config,nullptr,nullptr,nullptr,&frames,true))==expected,
            "Candidate parallelism changed output geometry, models, ordering or quality");
  }
  threads("4");
  auto first=std::async(std::launch::async,[&]{return snapshot(detector.detect(input,config,nullptr,nullptr,nullptr,&frames,true));});
  auto second=std::async(std::launch::async,[&]{return snapshot(detector.detect(input,config,nullptr,nullptr,nullptr,&frames,true));});
  require(first.get()==expected && second.get()==expected,"Concurrent detection calls interfered");
  weld::DebugTrace trace;const auto caller=std::this_thread::get_id();int callbacks=0;
  trace.onStep=[&](const weld::DebugStep &){require(std::this_thread::get_id()==caller,"Debug callback moved threads");++callbacks;return false;};
  config.debugFlag=true;
  require(snapshot(detector.detect(input,config,nullptr,nullptr,&trace,&frames,true))==expected && callbacks>0,
          "Serial debug fallback changed results or callbacks");
  std::cout<<"bead="<<bead<<" method="<<int(method)<<" serial/2/4/8/concurrent/debug identical\n";
}
}
int main() try {
  RestoreEnvironment restore;
  poolTests();
  for(bool bead:{false,true})for(auto method:{weld::ProfileLocalizationMethod::CurvatureCandidates,weld::ProfileLocalizationMethod::CurveLineIntersection})
    detectionTests(bead,method);
  return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
