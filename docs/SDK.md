# WeldSeam 焊缝提取库

WeldSeam 是一个 C++ 点云算法库，输入原始 XYZ 点云，返回左侧壁、右侧壁、焊道面三个点云，并保留焊缝轨迹输出。交付内容为 DLL、链接库和头文件，不包含独立 EXE，不依赖 Qt 或 VTK。点云文件的读取、保存和显示由调用方处理。

调用时使用连续曲率算法的内置默认参数，无需读取或修改 JSON 配置文件。剖面分箱 `profileBinSize` 默认 `0.8f`，两套算法一致。

`dev` 分支的左右轨迹独立选择路径、平滑和裁剪，允许点数及起止帧不同；单侧有效时也可返回成功。调用方应分别使用各侧帧号和质量数组，不能再用相同数组下标配对左右点。`sideEndpointTrimmedStart/End[0/1]` 分别记录左右裁剪点数；旧 `endpointTrimmedStart/End` 保留为左侧计数。结果结构新增字段，须使用配套头文件和 DLL 重新编译。单侧轨迹有效且侧壁有实测支撑时，也可从该侧壁脚向内提取连通底面；范围由当前截面的实测轮廓决定，遇到明显缺口或陡立边界停止，不生成缺失的另一侧轨迹或侧壁。若两侧站位不一致，统一宽度及单/双道转换跳过，保留独立轨迹并在结果说明中注明。

当前 SDK 版本为 **5.0.0**。新增入口采样参数 `inputVoxelSize` 和采样点数统计，参数和检测结果结构体布局发生变化，调用工程必须使用配套的新头文件、LIB 和 DLL 重新编译。保留左、中、右三个 B 样条曲面模型。

默认 `limitSurfaceEndsToTracks=true`：仅左右侧壁分别受各自轨迹终点限制，底面首尾独立按实测支撑保留，不随较短侧壁截断。曲面模型使用各自保留的实测点拟合。左右终点沿用独立端部变化判断，不使用固定帧号或孔洞位置。设为 `false` 可使侧壁也保留轨迹终点之外仍有测量支持的部分，不改变底面。被裁剪截面的边界通过本帧的大面侧壁和实测壁脚重新确认，底面沿连续测量支撑提取。端部只有单侧有效轨迹时，从有侧壁拟合支撑或可靠跟踪法向的壁脚向内追踪底面，不再按双侧公共范围删除。追踪遇到明显缺口或陡立边界停止，不生成缺失侧壁或轨迹。双侧覆盖区内部短暂丢失一侧壁脚时，仍仅在相邻有效截面之间按站位插值边界，跨度不超过 `max(4.5*sliceStep, profileCurvatureRadius)`，防止扩展到外部平板。分离的面及缺失帧不补齐。调用方应分别使用各输出的实际范围。新增配置字段须配套重新编译调用工程。

SDK 3.0 的浅槽端部面提取修复已包含在本版中。不能仅替换旧版 DLL 升级到 SDK 5。

`surfaceRimRepairLength` 为侧壁外缘局部凹口重建最大跨度，默认 `24`，设为 `0` 关闭，单位与输入一致。正值仅影响左右 B 样条模型：以两端相容的实测外缘及可靠壁脚建立低权重补面约束，不改底面模型、轨迹、实测面点云、实测样本计数或误差统计口径。只使用原始外缘作支撑，不递归扩张；首尾、显式无效边界、跨度超限或两侧形状不相容的缺口不补。最大横向补面长度还受 `min(surfaceRimRepairLength/2, 2*profileCurvatureRadius)` 限制。真实凹口和扫描缺失可能不可区分，需保留真实凹口时调用方应设为 `0`。该字段新增到 `DetectorConfig`，调用工程须使用配套头文件、LIB 和 DLL 重新编译。

侧壁下部弯曲而偏离大面拟合平面时，可由当前截面的可靠壁脚法向与已识别侧壁双端约束，沿连续实测点恢复短过渡带。搜索长度限于 `max(4*max(0.5, profileBinSize), profileCurvatureRadius)`，明显测量缺口停止；未能连到已有侧壁时不恢复。只输出原始点，保留底面和已有面分类，新增点不作为再次扩张的种子。无可靠壁脚、被拒绝的边界及分离平行面不使用此恢复。

端部裁剪前会验证疑似偏离的尾段是否仍有连续的大面侧壁与实测底面接触支撑。若前段能够接上，且原候选持续偏离真实壁脚，则逐侧恢复实测壁脚后重新判断终点；真实弯折仍受原几何阈值限制。短暂无效接触点省略，持续缺口终止延伸。恢复点带 `QualityMeasuredWallFoot` 标志，`supportDistance=0`，坐标来自该帧实测点，旧候选的局部拟合约束不再标为有效；各侧质量、帧号和轨迹数组保持对齐。调试输出 `measuredWallFeet` 记录供验证的壁脚点（纵向为检测站位）。

入口采样默认 `config.inputVoxelSize=0.5f`，设为 `0` 关闭；单位与输入点云一致。提供 `frames` 时逐帧采样，保留帧号和每帧首尾点；每个体素选择最接近中心的实测点，输出仍为原始坐标。主体筛选使用采样点；面提取与拟合恢复主体保留体素内的原始密集点。有帧号时，轨迹定位和局部精修也恢复这些原始轮廓。关闭主体筛选时，面提取使用完整原始点云。帧站位、帧间距和扫描平面始终从原始帧计算，避免采样密度改变纵向尺度。无帧号时，轨迹检测仍使用采样点云，面提取使用恢复后的密集实测点。启用时跳过旧 `voxelSize` 二次降采样。该参数不同于仅用于主体分析的 `hullSurfaceVoxelSize`；较大体素仍可能改变主体范围，不能保证各档轨迹完全相同。比较轨迹应按原始帧号对齐，而不是直接按数组下标。

方向默认自动判断：`autoLongitudinalAxis=true` 在主体提取后做 PCA，按最大绝对分量选择 X/Y/Z 延伸轴。有原始帧时，优先使用保留主体的帧中心 PCA 推进方向，避免短扫描的横向宽度及倾斜扫描平面造成误判；无帧时使用主体点云主方向。`autoOpeningDirection=true` 对多个截面做鲁棒凹凸拟合及投票，选择槽口朝向。主方向不明确或槽口证据不足时，分别使用 `longitudinalAxis`、`interiorOpeningDirection` 备用值。需要固定方向时必须把对应自动开关设为 `false`。此功能保持输入坐标，仍使用现有轴向截面定义，不提供任意 PCA 坐标旋转。新增配置/结果字段后，调用方必须配套更新头文件与 DLL 并重新编译。

## 交付目录

将整个 `install` 目录交给调用方，主要内容如下：

```text
install/
├─ bin/                    WeldSeam.dll 及运行所需 DLL
│  └─ config/              随库附带的配置文件，默认调用无需使用
├─ include/WeldSeam/       公共头文件
├─ lib/
│  ├─ WeldSeam.lib         Windows 导入库
│  └─ cmake/WeldSeam/      CMake 集成文件
└─ README.md               接入说明
```

## 开发环境

当前交付版本使用 Windows x64、MSVC v142、Release `/MD`、C++14、PCL 1.11.0 和 Eigen 3.3.7。调用工程需使用匹配的编译环境及依赖版本。

编译时需要 PCL 及其配套开发依赖；算法库使用 PCL 的 `common`、`filters` 模块，不需要 VTK 可视化模块。

## 调用示例

在已有点云处理代码中包含头文件，直接传入点云即可：

```cpp
#include <WeldSeam.h>
#include <string>

bool extractWeldClouds(const weld::CloudT::ConstPtr& input,
                      weld::SeamClouds& output,
                      std::string& error)
{
    const bool generateSurfaceModels = true;
    return weld::extractSeams(input, output, generateSurfaceModels,
                              weld::defaultDetectorConfig(), &error);
}
```

第三个参数控制是否生成曲面模型：传 `false`（或 `0`）时只生成轨迹和三个实测面点云，
跳过耗时的 B 样条拟合，此时三个 `SurfaceModel` 均为空且 `message == "not_fitted"`；
传 `true`（或 `1`）时点云和曲面模型都生成。不传该参数的旧调用方式保持原行为，即生成曲面模型。

`weld::CloudT` 即 `pcl::PointCloud<pcl::PointXYZ>`。返回点云如下：

| 输出成员 | 内容 |
| --- | --- |
| `output.leftSurface` | 左侧壁实测点云 |
| `output.rightSurface` | 右侧壁实测点云 |
| `output.weldSurface` | 左右侧壁之间的焊道面实测点云 |
| `output.leftSurfaceModel` | 左侧壁 B 样条曲面模型 |
| `output.weldSurfaceModel` | 中间焊道 B 样条曲面模型 |
| `output.rightSurfaceModel` | 右侧壁 B 样条曲面模型 |
| `output.left` / `output.right` | 默认双道模式的左右焊缝轨迹 |
| `output.center` | 单道模式的中心轨迹，默认双道模式下为空 |

算法沿长度分段拟合左右侧壁大面，提取贴面的实测点，并将两个侧壁面之间的剩余主体点提取为焊道面。三个面保留输入点的原始坐标和顺序，彼此不重复使用同一个输入点；左右按横向坐标区分，默认 Y 较小的一侧为左。只提取有效检测区段内有几何支撑的面，不跨缺帧外推，也不包含外部上沿。默认输入单位为毫米。

浅槽端部已有轨迹、但侧壁不足以拟合大面时，改在当前帧壁脚附近拟合有界的短过渡面，并提取两侧壁脚之间的实测底面。短面必须有原始点支撑，且相对局部底面向外升起；不会把平坦底板强行分成三个面，也不会复制邻帧点补齐。这样有测量支撑的浅槽端部可以与轨迹保持相同的帧范围，无需缩小剖面分箱或截短轨迹。

当左右实测壁脚轨迹均存在、但其中一侧的侧壁点太少而无法可靠拟合时，该侧壁允许为空，底面仍按本截面的两个实测壁脚提取；至少需要另一侧有可靠侧壁支持，因此不会仅凭平板上的两个给定点虚构槽底或缺失侧壁。

调用成功返回 `true`；三个面支撑不足时可能部分为空，需要完整三面时检查相应点云。调用失败返回 `false`，清空全部输出，通过 `error` 获取原因。

## B 样条曲面输出

三个模型在实测面点云分类后独立拟合。使用 Eigen 稀疏最小二乘、三次张量积 B 样条、二阶差分平滑和三次鲁棒重加权，不增加 PCL surface / VTK 运行依赖。原始六个点云输出保持原义。

```cpp
Eigen::Vector3d point, normal;
const auto& model = output.weldSurfaceModel;
for (std::size_t patch = 0; patch < model.patches.size(); ++patch) {
    if (model.evaluate(patch, 0.5, 0.5, point, &normal)) {
        // point 为原始坐标系中的 XYZ；normal 为参数方向决定的单位法向。
    }
}
```

每片函数为 `S(u,v) = sum_i sum_j N_i,3(u) N_j,3(v) P_ij`，参数域是 `[0,1] × [0,1]`。`u` 按扫描顺序的截面重心累计弦长归一化；`v` 按截面轮廓累计弦长归一化。两侧壁用截面内的主方向排序，支持竖直侧壁。法向符号随参数方向确定，不保证朝向工件外侧。

`patches` 中每片保存完整节点向量 `knotsU/knotsV`（端点重复四次）、控制网格大小 `countU/countV`、三维控制点 `controlPoints[i*countV+j]`、起止帧/站位及 `rmsError/maxError`。误差是在分配的扫描/轮廓参数处相对实测点的欧氏残差，不是最近点到曲面距离。拟合在原始坐标单位下计算；空间切片模式 `firstFrame/lastFrame=-1`，站位是检测器工作站位，帧模式下不一定等于原始 X。

同一面内部的缺口按实际站位距离判断：两端间距不超过 `max(4.5*sliceStep, profileCurvatureRadius)` 时，用两端实测轮廓建立低权重插值约束，保持曲面连续（帧模式下 `sliceStep` 使用原始帧间距）；截面内部通过等弦长插值补充拟合样本。更长缺口、显式无效的面边界仍会分片，不向首尾外推。每片至少需要四个有效实测截面、每截面四个不同点及二维几何支撑。误差 RMS 超过 `2*max(0.01, profileWallTolerance, hullSurfaceDistance)` 的长区段会在实测截面处继续细分；无法满足误差条件的短片不输出。共享边界截面的细分片固定为相同边界曲线，保证位置连续（C0），不保证法向连续。`rejectedPatchCount` 记录最终被拒绝的非空子区段数，`message` 为 `ok`、`partial_support`、`insufficient_support_or_fit_error`，未执行拟合时为 `not_fitted`。`valid()` 只表示至少有一个可用片，不表示完整覆盖整个面。曲面拟合失败不使原有焊缝检测失败。求值越界、参数非有限或模型结构无效时返回 `false` 且不修改输出；请求法向时退化点也返回 `false`。

拟合保留区段内全部有效截面，每截面使用 48～192 个等弦长样本平衡点密度；控制网格最多 32×12。跨截面插值样本的初始权重为实测轮廓样本的 1/4；误差和 `sampleCount` 仅在有效截面的实测点上统计。面提取先恢复主体保留体素内的原始密集点，再根据同截面、同一面的固定轮廓邻域恢复因局部法向不稳定而漏选的实测点。合成约束只用于拟合，不混入三个实测面点云。适用于截面可按一个主方向排序的槽壁及底面；明显回折、多层重叠轮廓需要先拆分。侧壁局部缺点造成的壁脚边缘凹缺，可使用可靠壁脚和同侧邻近截面建立低权重补面约束；补面距离受 `max(4.5*sliceStep, profileCurvatureRadius)` 限制，较弯曲的缺口还须由两端邻近实测壁脚确认。优先沿可靠壁脚参数化，减少缺点引起的边缘波动；不跨到另一侧壁，不补无可靠边界的长缺口。补面只影响拟合模型，原始面点云不增加合成点，`sampleCount` 和误差仍仅统计实测点。三个不同面之间不强制严格共边，不跨长缺口连网格；每片的 `v=0/1` 是拟合的测量边界，可有拟合误差。

GUI 默认显示蓝、橙、绿三色 B 样条三角网格。“B 样条实体曲面”开关控制网格，“叠加实测面点云”用于对照；各面的颜色开关同时控制该面的网格和点云。未生成模型的面回退显示实测点云。

GUI/CLI 保存结果时增加三个文件：`*_left_surface_model.json`、`*_weld_surface_model.json`、`*_right_surface_model.json`。每个文件保存独立模型、参数约定和误差；模型不可用时仍覆盖写入 `valid:false` 与空 `patches`，避免残留旧模型。`*_params.json` 的 `surfaceModelFiles` 给出三个文件名。SDK 本身不执行文件 I/O，调用方可直接保存上述公开模型字段。

如果不需要错误信息，可以直接调用：

```cpp
weld::SeamClouds output;
if (weld::extractSeams(input, output)) {
    const auto& leftWall = output.leftSurface;
    const auto& rightWall = output.rightSurface;
    const auto& weldFace = output.weldSurface;
    if (!leftWall->empty() && !rightWall->empty() && !weldFace->empty()) {
        // 在这里使用三个面的 XYZ 点云。
    }
}
```

## 集成到现有工程

在调用工程的 `CMakeLists.txt` 中，定义自己的目标后添加：

```cmake
find_package(WeldSeam 5 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE WeldSeam::WeldSeam)
```

将 `your_target` 替换为已有的程序或库目标名称。配置工程时，将 `CMAKE_PREFIX_PATH` 指向本 SDK 的 `install` 目录，`PCL_DIR` 指向本机 PCL 的 `cmake` 目录。CMake 会传递所需头文件目录和链接依赖。

使用 Visual Studio 手动集成时，将 `install/include/WeldSeam` 加入包含目录，将 `install/lib` 加入库目录并链接 `WeldSeam.lib`，同时配置匹配版本的 PCL 及配套依赖。

## 运行部署

将 `install/bin` 下的全部 DLL 放到调用程序所在目录，或将 `install/bin` 加入运行环境的 `PATH`。默认调用不需要复制或加载 `config` 目录。

主体提取中的法向估计和逐点归属查询会根据点数自动使用最多 4 个 CPU 工作线程，
小点云仍串行执行。截面粗/细轮廓复用排序，
纵向跨度不足的分箱跳过无效斜率配对；轮廓 B 样条取样只计算当前跨度内 4 个非零基函数，
不再用整列零基函数参与矩阵乘法。连续截面跟踪保持原顺序，算法阈值与输出精度不变。

正常提取还会在每次调用内复用固定线程池，并行计算同一截面的独立候选（含槽底样条拟合、
区域生长和几何验证）。默认总参与线程数不超过 4（包含调用线程，并受硬件线程数限制）；
少于 16 个候选的批次串行执行，结果严格按原候选顺序合并。前后帧侧壁跟踪和最终路径选择仍串行。
传入 `DebugTrace` 时，候选阶段回退到串行，调试回调始终在调用线程执行。
需要限核或性能对照时，可在启动程序前设置环境变量 `WELD_SEAM_CANDIDATE_THREADS=1..8`；
`1` 仅关闭候选阶段并行，不关闭主体提取已有的并行。未设置或值无效时使用默认值。
不要在提取过程中修改环境变量。线程预算按每次调用计算；调用方并发处理多份点云时应相应调小，
避免线程过量。无需修改公共参数结构，也不会修改进程级 Eigen 线程设置。

面点云提取的二维邻域法向计算在批量点数不少于 20000 时自动尝试 OpenCL GPU。
运行时动态加载系统 `OpenCL.dll`，优先选择非共享内存的独立 GPU；没有 OpenCL 平台、
设备不可用、内核编译/执行失败或批量过小时自动走原 CPU 路径，因此 SDK 不新增必须链接的
OpenCL 导入库，也不要求调用接口增加参数。成功使用时 `stageMilliseconds` 会增加
`gpu_surface_normals` 项；`surface_extraction` 仍包含完整的分组、分类和曲面模型拟合时间。
同一结果还会分别记录 `surface_cloud_classification` 和 `surface_model_fitting`，便于区分面点云
生成与最终三个曲面模型拟合；关闭曲面模型时不会产生后一项。
首次使用某个设备时驱动需要编译 OpenCL 内核，首轮可能出现额外延迟。

必须在第一次提取前设置以下可选环境变量，运行中修改不保证生效：

- `WELD_SEAM_GPU=0|off|cpu`：关闭 GPU，强制使用 CPU；
- `WELD_SEAM_GPU=1|on|opencl`：即使小批量也尝试 GPU，失败仍自动回退 CPU；
- `WELD_SEAM_OPENCL_DEVICE=<名称片段>`：指定设备，如 `NVIDIA` 或 `Intel`；
- `WELD_SEAM_GPU_DIAGNOSTICS=1`：向标准错误输出设备、点数和 GPU 子阶段耗时。

未设置 `WELD_SEAM_GPU` 时为自动模式。当前 GPU 范围仅包含面点云局部法向；主体提取、
截面候选/B 样条、帧间跟踪、最终轨迹平滑和曲面模型拟合仍使用 CPU。调用方并发提取时，
OpenCL 队列内部串行复用，CPU 其余阶段不受影响。

## 从源码构建

已有 SDK 的调用方可跳过本节。在源码根目录执行以下命令，构建并安装 SDK：

```powershell
cmake --preset sdk
cmake --build --preset sdk-release
ctest --preset sdk-release
cmake --install build/sdk --config Release --prefix install
```

构建前请按本机环境调整 `CMakePresets.json` 中的 PCL 路径。SDK 预设关闭 GUI 和命令行应用；安装结果位于 `install` 目录。
