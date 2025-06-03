#pragma once
#include <vector>
#include <map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include"header.h"
#include"DebugJson.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

namespace TeethConvex {
    class TeethCurveFitter {
    public:
        TeethCurveFitter(const std::map<int, std::vector<Point_2>>& pointsMap);
        void setPoints(const std::map<int, std::vector<Point_2>>& pointsMap);

        void computeConvexHulls();
        void computeSideConvexHull(const std::map<int, std::vector<Point_2>>& pointsMap,
            std::vector<Point_2>& firstTeethHull,
            std::vector<Point_2>& lastTeethHull);
        void reAssignMap(std::map<int, std::vector<Point_2>>& resultPointMap,
            const std::vector<Point_2>& result, const std::map<int, std::vector<Point_2>>& pointsMap);
        void OptimizeConvexHull(std::map<int, std::vector<Point_2>>& resultPointMap,
            const std::map<int, std::vector<Point_2>>& pointsMap);

        std::vector<Point_2> fitBSpline(const std::vector<Point_2>& points, size_t num_samples = 100);
        const std::vector<Point_2>& getInnerBSpline() const;
        const std::vector<Point_2>& getOuterBSpline() const;

        // 添加新的getter方法
        const std::vector<Point_2>& getInnerHull() const;
        const std::vector<Point_2>& getOuterHull() const;
        const std::map<int, std::vector<Point_2>>& getInnerHullMap() const;
        const std::map<int, std::vector<Point_2>>& getOuterHullMap() const;

        // 其他getter...

        static enum class returnTeeth {
            all,
            pre,
            current,
            next
        };
    private:
        std::vector<Point_2> points;
        std::map<int, std::vector<Point_2>> pointsMap;
        std::map<int, std::vector<Point_2>> rotatedPointsMap;
        std::vector<Point_2> innerhull, outerhull;
        std::map<int, std::vector<Point_2>> innerhullMap;
        std::map<int, std::vector<Point_2>> outerhullMap;
        std::vector<Point_2> innerBSplinePoints;
        std::vector<Point_2> outerBSplinePoints;
        int findEndPoint_low(const std::vector<Point_2>& points);
        int findEndPoint_high(const std::vector<Point_2>& points);
        Point_2 getClosestPoint(const std::vector<Point_2>& points, const Point_2& linePoint1, const Point_2& linePoint2);
        std::vector<Point_2> sortPointsAlongCurve(const std::vector<Point_2>& points, bool curve_right2left, bool low_high);
        float pointToLineDistance(Point_2 point, Point_2 linePoint1, Point_2 linePoint2);
        Point_2 CalculateCentroid(const std::vector<Point_2>& targetPoints);
        std::vector<Point_2> Rotate(const std::vector<Point_2>& targetPoints, Point_2 rotateCenter);
        vector<Point_2> filtSideTeeth(const std::vector<Point_2>& teethHull);
        // 其他私有成员和辅助函数
    };
}