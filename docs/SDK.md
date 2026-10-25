# WeldSeam 焊缝提取库

WeldSeam 是一个 C++ 点云算法库，输入原始 XYZ 点云，返回左侧壁、右侧壁、焊道面三个点云，并保留焊缝轨迹输出。交付内容为 DLL、链接库和头文件，不包含独立 EXE，不依赖 Qt 或 VTK。点云文件的读取、保存和显示由调用方处理。

调用时使用连续曲率算法的内置默认参数，无需读取或修改 JSON 配置文件。

当前 SDK 版本为 **4.0.0**。新增左、中、右三个 B 样条曲面模型，返回结构体布局发生变化，调用工程必须使用配套的新头文件、LIB 和 DLL 重新编译。

SDK 3.0 的浅槽端部面提取修复已包含在本版中。此前仅替换 DLL 的升级方式不适用于 SDK 4。

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
    return weld::extractSeams(input, output,
                              weld::defaultDetectorConfig(), &error);
}
```

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

缺帧、缺少某侧面点或空间站位缺口会分片。每片至少需要四个连续有效截面、每截面四个不同点及二维几何支撑。误差 RMS 超过 `2*max(0.01, profileWallTolerance, hullSurfaceDistance)` 的长区段会在实测截面处继续细分；无法满足误差条件的短片不输出。细分片可以共享一个边界截面，但不保证片间严格连续。`rejectedPatchCount` 记录最终被拒绝的非空子区段数，`message` 为 `ok`、`partial_support`、`insufficient_support_or_fit_error`，未执行拟合时为 `not_fitted`。`valid()` 只表示至少有一个可用片，不表示完整覆盖整个面。曲面拟合失败不使原有焊缝检测失败。求值越界、参数非有限或模型结构无效时返回 `false` 且不修改输出；请求法向时退化点也返回 `false`。

拟合以每片最多 256 个截面、每截面 48 个等弦长样本平衡点密度；控制网格最多 32×12，误差在全部保留实测点上统计。适用于截面可按一个主方向排序的槽壁及底面；明显回折、多层重叠轮廓需要先拆分。拟合仅覆盖当前提取出的面区域，不恢复被上游筛除的点，也不强制三面严格共边。不同片不可跨缺口连网格；每片的 `v=0/1` 是拟合的测量边界，可有拟合误差。

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
find_package(WeldSeam 4 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE WeldSeam::WeldSeam)
```

将 `your_target` 替换为已有的程序或库目标名称。配置工程时，将 `CMAKE_PREFIX_PATH` 指向本 SDK 的 `install` 目录，`PCL_DIR` 指向本机 PCL 的 `cmake` 目录。CMake 会传递所需头文件目录和链接依赖。

使用 Visual Studio 手动集成时，将 `install/include/WeldSeam` 加入包含目录，将 `install/lib` 加入库目录并链接 `WeldSeam.lib`，同时配置匹配版本的 PCL 及配套依赖。

## 运行部署

将 `install/bin` 下的全部 DLL 放到调用程序所在目录，或将 `install/bin` 加入运行环境的 `PATH`。默认调用不需要复制或加载 `config` 目录。

## 从源码构建

已有 SDK 的调用方可跳过本节。在源码根目录执行以下命令，构建并安装 SDK：

```powershell
cmake --preset sdk
cmake --build --preset sdk-release
ctest --preset sdk-release
cmake --install build/sdk --config Release --prefix install
```

构建前请按本机环境调整 `CMakePresets.json` 中的 PCL 路径。SDK 预设关闭 GUI 和命令行应用；安装结果位于 `install` 目录。
