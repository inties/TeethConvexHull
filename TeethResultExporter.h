#pragma once
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
namespace TeethConvex {
    struct Point2D {
        double x, y;
        int label;
    };

    class TeethResultExporter {
    public:
        TeethResultExporter() = default;
        ~TeethResultExporter() = default;

        // 在每次迭代结束后，调用此接口收集数据
        // 注意：样条曲线直接作为参数传入，不再由本类计算
        void addIterationData(
            const std::map<int, std::vector<Point_2>>& pointsMap,
            const std::map<int, std::pair<double, double>>& zRanges,
            const std::vector<Point_2>& innerHull,
            const std::vector<Point_2>& outerHull,
            const std::map<int, std::vector<Point_2>>& innerHullMap,
            const std::map<int, std::vector<Point_2>>& outerHullMap,
            const std::vector<Point_2>& innerBSpline,
            const std::vector<Point_2>& outerBSpline
        );

        // 所有迭代完成后，写 JSON
        void exportJSON(const std::string& jsonPath);

        // 所有迭代完成后，写 OBJ
        void exportOBJ(const std::string& objPath);
        void exportOBJ_2(const std::string& objPath);
		void exportSingleSplineOBJ(const std::string& objPath, bool isInner);

    private:
        // 内部存储每次迭代的数据
        struct IterationData {
            std::map<int, std::vector<Point_2>> pointsMap;       // 二维点集
            std::map<int, std::pair<double, double>> zRanges;    // minZ/maxZ
            std::vector<Point_2> innerHull, outerHull;           // 凸包点
            std::map<int, std::vector<Point_2>> innerHullMap, outerHullMap; // 分组凸包点
            std::vector<Point_2> innerBSpline, outerBSpline;     // 样条曲线点
            std::unordered_map<std::string, int> inPointLabel;   // 内侧spline点→label
            std::unordered_map<std::string, int> outPointLabel;  // 外侧spline点→label
        };
        std::vector<IterationData> m_iterations;

        // 质心相关成员变量
        Point_2 m_innerCentroid;  // 内侧质心，只计算一次
        Point_2 m_outerCentroid;  // 外侧质心，只计算一次
        bool m_centroidsComputed = false;  // 标记质心是否已计算

        // 将 {x,y} 字符串化，用于快速 lookup
        std::string makeKey(const Point_2& p) const {
            return std::to_string(p.x()) + "," + std::to_string(p.y());
        }

        // 重新分布样条曲线点，使其均匀分布
        void reAssignSpline(
            std::vector<Point_2>& innerBSpline,
            std::vector<Point_2>& outerBSpline
        );

        // 计算点集的质心
        Point_2 calculateCentroid(const std::vector<Point_2>& points) const;

        // 从质心发射射线并求与样条曲线的交点
        std::vector<Point_2> getIntersectionPoints(
            const Point_2& centroid,
            const std::vector<Point_2>& splinePoints,
            double angleStep = 2.0  // 默认2度间隔
        ) const;

        // 计算射线与样条曲线所有线段的交点，返回最近的交点
        bool raySplineIntersection(
            const Point_2& rayOrigin,
            double rayAngle,
            const std::vector<Point_2>& splinePoints,
            Point_2& intersection
        ) const;

        // 计算射线与单个线段的交点
        bool raySegmentIntersection(
            const Point_2& rayOrigin,
            double rayAngle,
            const Point_2& segStart,
            const Point_2& segEnd,
            Point_2& intersection
        ) const;

        // 角度转换为方向向量
        std::pair<double, double> angleToDirection(double angleInDegrees) const;

        // 从 splineMap（标签→样条点）建立 {key→label} 查找表
        std::unordered_map<std::string, int> createPointToLabelMap(
            const std::map<int, std::vector<Point_2>>& splineMap
        ) const;

        // 从 hullMap 和 splinePoints 计算 splineMap
        void assignLabelsToSplinePoints(
            const std::vector<Point_2>& splinePoints,
            const std::map<int, std::vector<Point_2>>& hullMap,
            std::map<int, std::vector<Point_2>>& splineMap
        );

        // 辅助把 CGAL::Point_2 + label 组合成 Point2D
        std::vector<Point2D> convertToPoint2D(
            const std::vector<Point_2>& splinePoints,
            const std::unordered_map<std::string, int>& pointLabelMap
        ) const;
    };
}