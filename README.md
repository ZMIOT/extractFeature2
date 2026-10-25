# WeldSeam 焊缝提取库

WeldSeam 是一个 C++ 点云算法库，输入原始 XYZ 点云，返回左侧壁、右侧壁、焊道面三个点云，并保留焊缝轨迹输出。交付内容为 DLL、链接库和头文件，不包含独立 EXE，不依赖 Qt 或 VTK。点云文件的读取、保存和显示由调用方处理。

调用时使用连续曲率算法的内置默认参数，无需读取或修改 JSON 配置文件。

当前 SDK 版本为 **5.0.0**。新增入口实测点采样参数 `inputVoxelSize`，参数和检测统计结构体布局发生变化，调用工程必须使用配套的新头文件、LIB 和 DLL 重新编译。保留 `leftSurfaceModel`、`weldSurfaceModel`、`rightSurfaceModel` 三个 B 样条曲面输出，可求坐标和法向。用法、分片与有效范围见 [SDK 曲面说明](docs/SDK.md#b-样条曲面输出)。

SDK 3.0 的浅槽端部面提取修复已包含在本版中。不能仅替换旧版 DLL 升级到 SDK 5。GUI/CLI 导出 `*_left_surface_model.json`、`*_weld_surface_model.json`、`*_right_surface_model.json`，保存完整 B 样条函数。

`inputVoxelSize=0` 关闭入口采样，也是默认值。GUI 在“关键参数 → 入口降采样体素（0 关闭）”手动输入 size，支持最多六位小数，范围与配置一致为 0..1000000000，例如 `0.05`、`0.075`、`0.1`；单位与点云一致。每帧独立选择最靠近体素中心的原始测量点，保留每帧首尾点、原始坐标、帧号与点顺序，不生成体素均值点。主体提取、候选生成、局部精修及三个表面均使用采样后的点云；启用时跳过旧的 `voxelSize` 二次降采样。采样可能影响轨迹和曲面拟合，应对照未采样结果评估后选值。CLI 使用 `--input-voxel 0.075`，JSON 配置使用 `"inputVoxelSize": 0.075`；示例数值不代表已验证的精度保证。

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
| `output.left` / `output.right` | 默认双道模式的左右焊缝轨迹 |
| `output.center` | 单道模式的中心轨迹，默认双道模式下为空 |

算法沿长度分段拟合左右侧壁大面，提取贴面的实测点，并将两个侧壁面之间的剩余主体点提取为焊道面。三个面保留输入点的原始坐标和顺序，彼此不重复使用同一个输入点；左右按横向坐标区分，默认 Y 较小的一侧为左。只提取有效检测区段内有几何支撑的面，不跨缺帧外推，也不包含外部上沿。默认输入单位为毫米。

浅槽端部已有轨迹、但侧壁不足以拟合大面时，改在当前帧壁脚附近拟合有界的短过渡面，并提取两侧壁脚之间的实测底面。短面必须有原始点支撑，且相对局部底面向外升起；不会把平坦底板强行分成三个面，也不会复制邻帧点补齐。这样有测量支撑的浅槽端部可以与轨迹保持相同的帧范围，无需缩小剖面分箱或截短轨迹。

调用成功返回 `true`；三个面支撑不足时可能部分为空，需要完整三面时检查相应点云。调用失败返回 `false`，清空全部输出，通过 `error` 获取原因。

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

## 从源码构建

已有 SDK 的调用方可跳过本节。在源码根目录执行以下命令，构建并安装 SDK：

```powershell
cmake --preset sdk
cmake --build --preset sdk-release
ctest --preset sdk-release
cmake --install build/sdk --config Release --prefix install
```

构建前请按本机环境调整 `CMakePresets.json` 中的 PCL 路径。SDK 预设关闭 GUI 和命令行应用；安装结果位于 `install` 目录。

本仓库只保存源码、构建配置及必要的 SDK 说明和许可证，不包含点云、编译产物或安装包。部分回归测试依赖本地 `tests/data/*.txt` 实测数据；运行完整测试前需恢复这些数据，数据未随仓库上传。
