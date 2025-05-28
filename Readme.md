# Readme

## 功能说明

​	根据输入的口扫模型.obj文件，通过在多个平面（目前只支持在3个平面上计算）上计算轴嵴连线，生成阶梯状包围盒（目前只能生成两个侧面），用以约束基台形变。

​	效果示例：

​				

​						<img src="E:\MylabProjects\ConvexHull\convex\convex1\轴嵴连线文档.assets\image-20250410154726920-1744271247387-1.png" alt="image-20250410154726920" style="zoom: 80%;" />



<img src="E:\MylabProjects\ConvexHull\convex\convex1\轴嵴连线文档.assets\image-20250410154743419.png" alt="image-20250410154743419" style="zoom:50%;" />

<img src="E:\MylabProjects\ConvexHull\convex\convex1\轴嵴连线文档.assets\image-20250410154943646.png" alt="image-20250410154943646" style="zoom:80%;" />



## 接口说明

convex.h提供了两个接口：

`bool processTeethData(const std::string& input_obj);`

`bool processTeethData(const std::string& input_obj, float percent1, float percent2, float percent3);`

参数说明：

* .obj 文件。包含上半口或下半口扫描模型的点云数据，Z 轴已对齐咬合面法线朝向。标签文件(.json格式)需与.obj文件在同一根目录下，遵循 FDI 命名规范（例如 11、12、21 等）。

* 三个浮点数，指定切割平面的位置。可指定3个0到1之间的浮点数，用于指定三个平面的位置，轴嵴连线的计算将在这三个平面上进行。默认为1，0.7，0.5。数值代表牙齿高度在平面以上的部分占整颗牙齿高度的比例。当该平面为1时，切割平面在牙齿底部，为0.5时，切割平面位于牙齿中间。



输出：

- .obj格式的文件。在输入的.obj文件目录下，输出另一个.obj格式文件，包含轴嵴连线求出的阶梯状包围盒（目前能求出包围盒的两个侧面）



## 环境和依赖

- **CGAL**：提供二维凸包计算功能（convex_hull_2）。已经放在dependency文件夹中。
- **nlohmann/json**：用于解析 .json 文件和生成输出 JSON 数据。已经放在dependency文件夹中。
- **tinyspline:** 用于样条曲线拟合，头文件已经放在工程目录下，实现放在工程中一起编译(parson.c，tinyspline.c，tinysplinecxx.cxx)。
- 使用C++14



## 注意事项

* 舌侧的轴嵴连线准确性依赖于牙齿分割准确性（颊侧轴嵴连线不受影响）。牙齿分割不准确可能导致问题。

  



------

## 凸包算法求轴嵴连线思路

### 问题定义

轴嵴线点是指牙齿在舌侧或颊侧的突出点，轴嵴横向连线是通过连接不同牙齿的轴嵴线点形成的曲线，用于描述牙齿排列规律和几何约束条件。程序的目标是分别生成舌侧段和颊侧段的轴嵴横向连线。

### 凸包算法

凸包是包含给定点集的最小凸多边形，能够有效描述点集的外围边界。在本程序中，凸包用于近似牙齿点云的外部轮廓，从而提取轴嵴线点。

### 算法思路简介

1. 数据预处理

   - 将牙齿的 XYZ 坐标投影到 XY 平面（忽略 Z 轴），形成二维点集。
   - 根据 .json 文件的标签，将点分组到 pointsMap，并滤除标签为 0 的点。
   - 重构标签，使其从 1 开始连续编号。

2. 旋转点集

   - 对每颗牙齿的点集，计算其质心并绕质心旋转 180°。
   - 原始点集存储在 pointsMap，旋转后的点集存储在 rotatedPointsMap。

3. 计算凸包

   - **颊侧凸包（外侧）**：对 pointsMap 中的点集使用 CGAL::convex_hull_2 计算凸包，得到颊侧轴嵴线点。
   - **舌侧凸包（内侧）**：对 rotatedPointsMap 计算凸包，结果旋转回原始位置，得到舌侧轴嵴线点。

4. 优化凸包点

   - **首尾牙齿单独处理**：通过 computeSideConvexHull 计算两端牙齿的局部凸包。
   - **补充缺失点**：调用 OptimizeConvexHull，检查哪些牙齿未被凸包覆盖，在这些牙齿的点集中找到与相邻凸包点连线最近的点，作为补充点。
   - **点序排列**：使用 sortPointsAlongCurve 将凸包点按轴嵴连线方向顺序排列，便于后续拟合。

5. 样条曲线拟合

   - 样条曲线（如 B 样条或贝塞尔曲线）拟合优化后的凸包点，生成平滑的轴嵴横向连线。

### 具体实现细节

- **颊侧凸包（外侧）**：直接对 pointsMap 计算凸包，排除首尾牙齿干扰，结果存储在 outerTeethHull。
- **舌侧凸包（内侧）**：对 rotatedPointsMap 计算凸包，若显示原始状态，则将结果旋转回原位，存储在 innerTeethHull。
- **优化（OptimizeConvexHull）**：通过 getClosestPoint 找到最近点，确保每颗牙齿在凸包中有代表点。

