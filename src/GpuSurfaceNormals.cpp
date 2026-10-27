#include "GpuSurfaceNormals.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace weld {
namespace detail {
namespace {

bool equalsIgnoreCase(const char *value, const char *expected) {
  if (!value || !expected) return false;
  while (*value && *expected) {
    if (std::tolower(static_cast<unsigned char>(*value)) !=
        std::tolower(static_cast<unsigned char>(*expected))) return false;
    ++value; ++expected;
  }
  return *value == *expected;
}

bool gpuDisabled() {
  const char *value = std::getenv("WELD_SEAM_GPU");
  return value && (equalsIgnoreCase(value, "0") ||
                   equalsIgnoreCase(value, "off") ||
                   equalsIgnoreCase(value, "false") ||
                   equalsIgnoreCase(value, "cpu"));
}

bool gpuForced() {
  const char *value = std::getenv("WELD_SEAM_GPU");
  return value && (equalsIgnoreCase(value, "1") ||
                   equalsIgnoreCase(value, "on") ||
                   equalsIgnoreCase(value, "true") ||
                   equalsIgnoreCase(value, "opencl"));
}

bool diagnosticsEnabled() {
  const char *value = std::getenv("WELD_SEAM_GPU_DIAGNOSTICS");
  return value && !equalsIgnoreCase(value, "0") &&
         !equalsIgnoreCase(value, "off") &&
         !equalsIgnoreCase(value, "false");
}

#ifdef _WIN32

using cl_uchar = std::uint8_t;
using cl_int = std::int32_t;
using cl_uint = std::uint32_t;
using cl_ulong = std::uint64_t;
using cl_bool = cl_uint;
using cl_bitfield = cl_ulong;
using cl_device_type = cl_bitfield;
using cl_mem_flags = cl_bitfield;
using cl_platform_info = cl_uint;
using cl_device_info = cl_uint;
using cl_program_build_info = cl_uint;
using cl_context_properties = std::intptr_t;
struct _cl_platform_id; using cl_platform_id = _cl_platform_id *;
struct _cl_device_id; using cl_device_id = _cl_device_id *;
struct _cl_context; using cl_context = _cl_context *;
struct _cl_command_queue; using cl_command_queue = _cl_command_queue *;
struct _cl_mem; using cl_mem = _cl_mem *;
struct _cl_program; using cl_program = _cl_program *;
struct _cl_kernel; using cl_kernel = _cl_kernel *;
struct _cl_event; using cl_event = _cl_event *;

constexpr cl_int CL_SUCCESS = 0;
constexpr cl_bool CL_TRUE = 1;
constexpr cl_device_type CL_DEVICE_TYPE_GPU = 1u << 2;
constexpr cl_platform_info CL_PLATFORM_NAME = 0x0902;
constexpr cl_device_info CL_DEVICE_MAX_COMPUTE_UNITS = 0x1002;
constexpr cl_device_info CL_DEVICE_AVAILABLE = 0x1027;
constexpr cl_device_info CL_DEVICE_COMPILER_AVAILABLE = 0x1028;
constexpr cl_device_info CL_DEVICE_NAME = 0x102B;
constexpr cl_device_info CL_DEVICE_VENDOR = 0x102C;
constexpr cl_device_info CL_DEVICE_HOST_UNIFIED_MEMORY = 0x1035;
constexpr cl_program_build_info CL_PROGRAM_BUILD_LOG = 0x1183;
constexpr cl_mem_flags CL_MEM_READ_ONLY = 1u << 2;
constexpr cl_mem_flags CL_MEM_WRITE_ONLY = 1u << 1;
constexpr cl_mem_flags CL_MEM_COPY_HOST_PTR = 1u << 5;

#ifndef CL_API_CALL
#define CL_API_CALL __stdcall
#endif

struct Float2 { float x, y; };

class OpenClRuntime {
public:
  OpenClRuntime() { initialize(); }
  OpenClRuntime(const OpenClRuntime &) = delete;
  OpenClRuntime &operator=(const OpenClRuntime &) = delete;
  ~OpenClRuntime() {
    if (kernel_ && releaseKernel_) releaseKernel_(kernel_);
    if (program_ && releaseProgram_) releaseProgram_(program_);
    if (queue_ && releaseCommandQueue_) releaseCommandQueue_(queue_);
    if (context_ && releaseContext_) releaseContext_(context_);
    if (module_) FreeLibrary(module_);
  }

  bool available() const { return ready_; }
  const std::string &deviceName() const { return deviceName_; }

  bool compute(const std::vector<Float2> &points,
               const std::vector<cl_int> &starts,
               const std::vector<cl_int> &ends, float radiusSquared,
               std::vector<Float2> &normals,
               std::vector<cl_uchar> &supported) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!ready_ || points.empty() || starts.size() != points.size() ||
        ends.size() != points.size()) return false;
    normals.assign(points.size(), Float2{0,0});
    supported.assign(points.size(), 0);
    if (points.size() > std::size_t(std::numeric_limits<cl_int>::max())) return false;
    cl_int error = CL_SUCCESS;
    cl_mem pointBuffer = createBuffer_(context_,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        points.size()*sizeof(Float2),const_cast<Float2 *>(points.data()),&error);
    if (error != CL_SUCCESS || !pointBuffer) return false;
    cl_mem startBuffer = createBuffer_(context_,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        starts.size()*sizeof(cl_int),const_cast<cl_int *>(starts.data()),&error);
    cl_mem endBuffer = createBuffer_(context_,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        ends.size()*sizeof(cl_int),const_cast<cl_int *>(ends.data()),&error);
    cl_mem normalBuffer = createBuffer_(context_,CL_MEM_WRITE_ONLY,
        normals.size()*sizeof(Float2),nullptr,&error);
    cl_mem supportedBuffer = createBuffer_(context_,CL_MEM_WRITE_ONLY,
        supported.size()*sizeof(cl_uchar),nullptr,&error);
    auto release = [&] {
      if (supportedBuffer) releaseMemObject_(supportedBuffer);
      if (normalBuffer) releaseMemObject_(normalBuffer);
      if (endBuffer) releaseMemObject_(endBuffer);
      if (startBuffer) releaseMemObject_(startBuffer);
      if (pointBuffer) releaseMemObject_(pointBuffer);
    };
    if (!startBuffer || !endBuffer || !normalBuffer || !supportedBuffer ||
        error != CL_SUCCESS) { release(); return false; }
    const cl_int count = static_cast<cl_int>(points.size());
    error = setKernelArg_(kernel_,0,sizeof(pointBuffer),&pointBuffer);
    error |= setKernelArg_(kernel_,1,sizeof(startBuffer),&startBuffer);
    error |= setKernelArg_(kernel_,2,sizeof(endBuffer),&endBuffer);
    error |= setKernelArg_(kernel_,3,sizeof(normalBuffer),&normalBuffer);
    error |= setKernelArg_(kernel_,4,sizeof(supportedBuffer),&supportedBuffer);
    error |= setKernelArg_(kernel_,5,sizeof(radiusSquared),&radiusSquared);
    error |= setKernelArg_(kernel_,6,sizeof(count),&count);
    const std::size_t global = points.size();
    if (error == CL_SUCCESS)
      error=enqueueNDRangeKernel_(queue_,kernel_,1,nullptr,&global,nullptr,0,nullptr,nullptr);
    if (error == CL_SUCCESS) error=finish_(queue_);
    if (error == CL_SUCCESS)
      error=enqueueReadBuffer_(queue_,normalBuffer,CL_TRUE,0,
          normals.size()*sizeof(Float2),normals.data(),0,nullptr,nullptr);
    if (error == CL_SUCCESS)
      error=enqueueReadBuffer_(queue_,supportedBuffer,CL_TRUE,0,
          supported.size()*sizeof(cl_uchar),supported.data(),0,nullptr,nullptr);
    release();
    return error == CL_SUCCESS;
  }

private:
  template<class Function> bool load(Function &function,const char *name) {
    function=reinterpret_cast<Function>(GetProcAddress(module_,name));
    return function!=nullptr;
  }
  template<class T> bool deviceValue(cl_device_id device,cl_device_info key,T &value) {
    std::size_t size=0;
    if(getDeviceInfo_(device,key,0,nullptr,&size)!=CL_SUCCESS || size!=sizeof(T))return false;
    return getDeviceInfo_(device,key,sizeof(T),&value,nullptr)==CL_SUCCESS;
  }
  std::string platformString(cl_platform_id platform,cl_platform_info key) {
    std::size_t size=0;if(getPlatformInfo_(platform,key,0,nullptr,&size)!=CL_SUCCESS||!size)return {};
    std::vector<char> value(size);if(getPlatformInfo_(platform,key,size,value.data(),nullptr)!=CL_SUCCESS)return {};
    return std::string(value.data());
  }
  std::string deviceString(cl_device_id device,cl_device_info key) {
    std::size_t size=0;if(getDeviceInfo_(device,key,0,nullptr,&size)!=CL_SUCCESS||!size)return {};
    std::vector<char> value(size);if(getDeviceInfo_(device,key,size,value.data(),nullptr)!=CL_SUCCESS)return {};
    return std::string(value.data());
  }
  static bool containsIgnoreCase(std::string value,std::string wanted) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char c){return char(std::tolower(c));});
    std::transform(wanted.begin(),wanted.end(),wanted.begin(),[](unsigned char c){return char(std::tolower(c));});
    return value.find(wanted)!=std::string::npos;
  }

  void initialize() {
    module_=LoadLibraryW(L"OpenCL.dll");if(!module_)return;
    if(!load(getPlatformIDs_,"clGetPlatformIDs")||!load(getPlatformInfo_,"clGetPlatformInfo")||
       !load(getDeviceIDs_,"clGetDeviceIDs")||!load(getDeviceInfo_,"clGetDeviceInfo")||
       !load(createContext_,"clCreateContext")||!load(createCommandQueue_,"clCreateCommandQueue")||
       !load(createProgramWithSource_,"clCreateProgramWithSource")||!load(buildProgram_,"clBuildProgram")||
       !load(getProgramBuildInfo_,"clGetProgramBuildInfo")||!load(createKernel_,"clCreateKernel")||
       !load(createBuffer_,"clCreateBuffer")||!load(setKernelArg_,"clSetKernelArg")||
       !load(enqueueNDRangeKernel_,"clEnqueueNDRangeKernel")||!load(enqueueReadBuffer_,"clEnqueueReadBuffer")||
       !load(finish_,"clFinish")||!load(releaseMemObject_,"clReleaseMemObject")||
       !load(releaseKernel_,"clReleaseKernel")||!load(releaseProgram_,"clReleaseProgram")||
       !load(releaseCommandQueue_,"clReleaseCommandQueue")||!load(releaseContext_,"clReleaseContext"))return;
    cl_uint platformCount=0;
    if(getPlatformIDs_(0,nullptr,&platformCount)!=CL_SUCCESS||!platformCount)return;
    std::vector<cl_platform_id> platforms(platformCount);
    if(getPlatformIDs_(platformCount,platforms.data(),nullptr)!=CL_SUCCESS)return;
    struct Candidate {cl_platform_id platform=nullptr;cl_device_id device=nullptr;
      std::string name,vendor,platformName;cl_bool unified=CL_TRUE;cl_uint units=0;};
    std::vector<Candidate> candidates;
    for(auto platform:platforms) {
      cl_uint deviceCount=0;
      if(getDeviceIDs_(platform,CL_DEVICE_TYPE_GPU,0,nullptr,&deviceCount)!=CL_SUCCESS||!deviceCount)continue;
      std::vector<cl_device_id> devices(deviceCount);
      if(getDeviceIDs_(platform,CL_DEVICE_TYPE_GPU,deviceCount,devices.data(),nullptr)!=CL_SUCCESS)continue;
      for(auto device:devices) {
        cl_bool available=0,compiler=0;
        if(!deviceValue(device,CL_DEVICE_AVAILABLE,available)||!available||
           !deviceValue(device,CL_DEVICE_COMPILER_AVAILABLE,compiler)||!compiler)continue;
        Candidate candidate;candidate.platform=platform;candidate.device=device;
        candidate.name=deviceString(device,CL_DEVICE_NAME);candidate.vendor=deviceString(device,CL_DEVICE_VENDOR);
        candidate.platformName=platformString(platform,CL_PLATFORM_NAME);
        deviceValue(device,CL_DEVICE_HOST_UNIFIED_MEMORY,candidate.unified);
        deviceValue(device,CL_DEVICE_MAX_COMPUTE_UNITS,candidate.units);
        candidates.push_back(std::move(candidate));
      }
    }
    if(candidates.empty())return;
    const char *selection=std::getenv("WELD_SEAM_OPENCL_DEVICE");
    auto acceptable=[&](const Candidate &c){return !selection||!*selection||
      containsIgnoreCase(c.name,selection)||containsIgnoreCase(c.vendor,selection)||
      containsIgnoreCase(c.platformName,selection);};
    auto best=candidates.end();
    for(auto it=candidates.begin();it!=candidates.end();++it)if(acceptable(*it))
      if(best==candidates.end()||(best->unified&&!it->unified)||
         (best->unified==it->unified&&it->units>best->units))best=it;
    if(best==candidates.end())return;
    cl_int error=CL_SUCCESS;
    context_=createContext_(nullptr,1,&best->device,nullptr,nullptr,&error);
    if(error!=CL_SUCCESS||!context_)return;
    queue_=createCommandQueue_(context_,best->device,0,&error);
    if(error!=CL_SUCCESS||!queue_)return;
    static const char *source=R"CLC(
#define MAX_NEIGHBORS 64
__kernel void profile_normals(__global const float2 *points,
    __global const int *starts,__global const int *ends,
    __global float2 *normals,__global uchar *supported,
    const float radius2,const int count) {
  const int i=(int)get_global_id(0);if(i>=count)return;
  normals[i]=(float2)(0.0f,0.0f);supported[i]=(uchar)0;
  float nearestDistance[MAX_NEIGHBORS];int nearestIndex[MAX_NEIGHBORS];
  int found=0,furthest=0;float furthestDistance=-1.0f;const float2 p=points[i];
  for(int k=starts[i];k<ends[i];++k) {
    const float2 d=points[k]-p;const float distance=dot(d,d);if(distance>radius2)continue;
    if(found<MAX_NEIGHBORS) {
      nearestDistance[found]=distance;nearestIndex[found]=k;
      if(distance>furthestDistance||(distance==furthestDistance&&k>nearestIndex[furthest]))
        {furthest=found;furthestDistance=distance;}++found;
    } else if(distance<furthestDistance||(distance==furthestDistance&&k<nearestIndex[furthest])) {
      nearestDistance[furthest]=distance;nearestIndex[furthest]=k;furthest=0;furthestDistance=nearestDistance[0];
      for(int n=1;n<MAX_NEIGHBORS;++n)
        if(nearestDistance[n]>furthestDistance||
           (nearestDistance[n]==furthestDistance&&nearestIndex[n]>nearestIndex[furthest]))
          {furthest=n;furthestDistance=nearestDistance[n];}
    }
  }
  if(found<5)return;supported[i]=(uchar)1;
  float2 center=(float2)(0.0f,0.0f);for(int n=0;n<found;++n)center+=points[nearestIndex[n]];
  center/=convert_float(found);float xx=0.0f,xy=0.0f,yy=0.0f;
  for(int n=0;n<found;++n){const float2 d=points[nearestIndex[n]]-center;xx+=d.x*d.x;xy+=d.x*d.y;yy+=d.y*d.y;}
  const float inv=1.0f/convert_float(found);xx*=inv;xy*=inv;yy*=inv;
  const float discriminant=hypot(xx-yy,2.0f*xy);
  const float small=0.5f*(xx+yy-discriminant),large=0.5f*(xx+yy+discriminant);
  if(!isfinite(small)||!isfinite(large)||large<1.0e-5f||small>0.15f*(small+large))return;
  const float2 first=(float2)(-xy,xx-small),second=(float2)(yy-small,-xy);
  const float2 normal=dot(first,first)>=dot(second,second)?first:second;
  const float length=hypot(normal.x,normal.y);if(!isfinite(length)||length<=1.0e-12f)return;
  normals[i]=normal/length;
}
)CLC";
    program_=createProgramWithSource_(context_,1,&source,nullptr,&error);
    if(error!=CL_SUCCESS||!program_)return;
    error=buildProgram_(program_,1,&best->device,"-cl-std=CL1.2",nullptr,nullptr);
    if(error!=CL_SUCCESS) {
      if(diagnosticsEnabled()) {
        std::size_t length=0;getProgramBuildInfo_(program_,best->device,CL_PROGRAM_BUILD_LOG,0,nullptr,&length);
        std::vector<char> log(length?length:1);getProgramBuildInfo_(program_,best->device,CL_PROGRAM_BUILD_LOG,log.size(),log.data(),nullptr);
        std::fprintf(stderr,"WeldSeam OpenCL build failed: %s\n",log.data());
      }
      return;
    }
    kernel_=createKernel_(program_,"profile_normals",&error);if(error!=CL_SUCCESS||!kernel_)return;
    device_=best->device;deviceName_=best->name;ready_=true;
  }

  using ContextCallback=void (CL_API_CALL *)(const char *,const void *,std::size_t,void *);
  using ProgramCallback=void (CL_API_CALL *)(cl_program,void *);
  using GetPlatformIDs=cl_int (CL_API_CALL *)(cl_uint,cl_platform_id *,cl_uint *);
  using GetPlatformInfo=cl_int (CL_API_CALL *)(cl_platform_id,cl_platform_info,std::size_t,void *,std::size_t *);
  using GetDeviceIDs=cl_int (CL_API_CALL *)(cl_platform_id,cl_device_type,cl_uint,cl_device_id *,cl_uint *);
  using GetDeviceInfo=cl_int (CL_API_CALL *)(cl_device_id,cl_device_info,std::size_t,void *,std::size_t *);
  using CreateContext=cl_context (CL_API_CALL *)(const cl_context_properties *,cl_uint,const cl_device_id *,ContextCallback,void *,cl_int *);
  using CreateCommandQueue=cl_command_queue (CL_API_CALL *)(cl_context,cl_device_id,cl_bitfield,cl_int *);
  using CreateProgramWithSource=cl_program (CL_API_CALL *)(cl_context,cl_uint,const char **,const std::size_t *,cl_int *);
  using BuildProgram=cl_int (CL_API_CALL *)(cl_program,cl_uint,const cl_device_id *,const char *,ProgramCallback,void *);
  using GetProgramBuildInfo=cl_int (CL_API_CALL *)(cl_program,cl_device_id,cl_program_build_info,std::size_t,void *,std::size_t *);
  using CreateKernel=cl_kernel (CL_API_CALL *)(cl_program,const char *,cl_int *);
  using CreateBuffer=cl_mem (CL_API_CALL *)(cl_context,cl_mem_flags,std::size_t,void *,cl_int *);
  using SetKernelArg=cl_int (CL_API_CALL *)(cl_kernel,cl_uint,std::size_t,const void *);
  using EnqueueNDRangeKernel=cl_int (CL_API_CALL *)(cl_command_queue,cl_kernel,cl_uint,const std::size_t *,const std::size_t *,const std::size_t *,cl_uint,const cl_event *,cl_event *);
  using EnqueueReadBuffer=cl_int (CL_API_CALL *)(cl_command_queue,cl_mem,cl_bool,std::size_t,std::size_t,void *,cl_uint,const cl_event *,cl_event *);
  using Finish=cl_int (CL_API_CALL *)(cl_command_queue);
  using ReleaseMemObject=cl_int (CL_API_CALL *)(cl_mem);using ReleaseKernel=cl_int (CL_API_CALL *)(cl_kernel);
  using ReleaseProgram=cl_int (CL_API_CALL *)(cl_program);using ReleaseCommandQueue=cl_int (CL_API_CALL *)(cl_command_queue);
  using ReleaseContext=cl_int (CL_API_CALL *)(cl_context);
  HMODULE module_=nullptr;cl_device_id device_=nullptr;cl_context context_=nullptr;
  cl_command_queue queue_=nullptr;cl_program program_=nullptr;cl_kernel kernel_=nullptr;
  bool ready_=false;std::string deviceName_;std::mutex mutex_;
  GetPlatformIDs getPlatformIDs_=nullptr;GetPlatformInfo getPlatformInfo_=nullptr;
  GetDeviceIDs getDeviceIDs_=nullptr;GetDeviceInfo getDeviceInfo_=nullptr;
  CreateContext createContext_=nullptr;CreateCommandQueue createCommandQueue_=nullptr;
  CreateProgramWithSource createProgramWithSource_=nullptr;BuildProgram buildProgram_=nullptr;
  GetProgramBuildInfo getProgramBuildInfo_=nullptr;CreateKernel createKernel_=nullptr;
  CreateBuffer createBuffer_=nullptr;SetKernelArg setKernelArg_=nullptr;
  EnqueueNDRangeKernel enqueueNDRangeKernel_=nullptr;EnqueueReadBuffer enqueueReadBuffer_=nullptr;
  Finish finish_=nullptr;ReleaseMemObject releaseMemObject_=nullptr;ReleaseKernel releaseKernel_=nullptr;
  ReleaseProgram releaseProgram_=nullptr;ReleaseCommandQueue releaseCommandQueue_=nullptr;ReleaseContext releaseContext_=nullptr;
};

OpenClRuntime &runtime(){static OpenClRuntime value;return value;}
#endif
} // namespace

bool computeSurfaceNormalsGpu(
    const CloudT &source,const std::vector<std::vector<std::size_t>> &groups,
    int uAxis,int vAxis,float radius,GpuSurfaceNormals &result) {
  result=GpuSurfaceNormals();
  if(gpuDisabled()||!std::isfinite(radius)||radius<=0||uAxis<0||uAxis>2||
     vAxis<0||vAxis>2||uAxis==vAxis)return false;
  std::size_t count=0;
  for(const auto &group:groups) {
    if(group.size()>std::size_t(std::numeric_limits<std::int32_t>::max())||
       count>std::size_t(std::numeric_limits<std::int32_t>::max())-group.size())return false;
    count+=group.size();
  }
  if(!count||(!gpuForced()&&count<20000))return false;
#ifdef _WIN32
  auto &opencl=runtime();if(!opencl.available())return false;
  std::vector<Float2> points;points.reserve(count);std::vector<cl_int> starts,ends;
  starts.reserve(count);ends.reserve(count);
  for(const auto &group:groups) {
    const auto begin=static_cast<cl_int>(points.size());
    for(auto index:group) {if(index>=source.size())return false;const auto &p=source[index];points.push_back({p.data[uAxis],p.data[vAxis]});}
    const auto end=static_cast<cl_int>(points.size());
    for(std::size_t i=0;i<group.size();++i){starts.push_back(begin);ends.push_back(end);}
  }
  std::vector<Float2> normals;std::vector<cl_uchar> supported;
  const auto start=std::chrono::steady_clock::now();
  if(!opencl.compute(points,starts,ends,radius*radius,normals,supported))return false;
  result.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
  result.device=opencl.deviceName();result.normals.resize(groups.size());result.supported.resize(groups.size());
  std::size_t offset=0;
  for(std::size_t i=0;i<groups.size();++i)for(std::size_t j=0;j<groups[i].size();++j,++offset) {
    result.normals[i].emplace_back(normals[offset].x,normals[offset].y);result.supported[i].push_back(supported[offset]);
  }
  if(diagnosticsEnabled())std::fprintf(stderr,
      "WeldSeam OpenCL surface normals: device=%s points=%zu time_ms=%.3f\n",
      result.device.c_str(),count,result.milliseconds);
  return true;
#else
  (void)source;(void)groups;(void)count;return false;
#endif
}

} // namespace detail
} // namespace weld
