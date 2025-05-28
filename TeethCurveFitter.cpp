#include "TeethCurveFitter.h"
#include <CGAL/convex_hull_2.h>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <iostream>
#include "tinysplinecxx.h"
using namespace tinyspline;
namespace TeethConvex {
    inline Point_2 CalculateCentroid(std::vector<Point_2>& targetPoints) {
        double sumX = 0, sumY = 0;
        for (const auto& p : targetPoints) {
            sumX += p.x();
            sumY += p.y();
        }
        return Point_2(sumX / targetPoints.size(), sumY / targetPoints.size());
    }
    inline  std::vector<Point_2> Rotate(const std::vector<Point_2>& targetPoints, Point_2 rotateCenter) {
        std::vector<Point_2> resultPoints;
        for (auto& p : targetPoints) {
            Point_2 temp(-p.x() + 2 * rotateCenter.x(), -p.y() + 2 * rotateCenter.y());
            resultPoints.push_back(temp);
        }
        return resultPoints;
    }
    template<typename Iterator>
    static std::vector<Point_2> computeUnitConvex(Iterator it, TeethCurveFitter::returnTeeth op) {
        auto prevIt = std::prev(it);
        auto nextIt = std::next(it);
        int pointsNum = prevIt->second.size() + it->second.size() + nextIt->second.size();
        Point_2* points = new Point_2[pointsNum];
        int index = 0;

        for (auto point : prevIt->second) points[index++] = point;
        for (auto point : it->second) points[index++] = point;
        for (auto point : nextIt->second) points[index++] = point;

        std::vector<Point_2> temp, TeethHull;
        CGAL::convex_hull_2(points, points + pointsNum, std::back_inserter(temp));

        switch (op) {
        case TeethCurveFitter::returnTeeth::pre:
            for (auto point : temp) {
                if (std::find(prevIt->second.begin(), prevIt->second.end(), point) != prevIt->second.end())
                    TeethHull.push_back(point);
            }
            break;
        case TeethCurveFitter::returnTeeth::next:
            for (auto point : temp) {
                if (std::find(nextIt->second.begin(), nextIt->second.end(), point) != nextIt->second.end())
                    TeethHull.push_back(point);
            }
            break;
        default:
            break;
        }
        delete[] points;
        return TeethHull;
    }
    // 构造函数
    TeethCurveFitter::TeethCurveFitter(const std::map<int, std::vector<Point_2>>& pointsMap)
        : pointsMap(pointsMap) {
        // 初始化旋转点集
        for (const auto& pair : pointsMap) {
            Point_2 rotateCenter = CalculateCentroid(pair.second);
            rotatedPointsMap[pair.first] = Rotate(pair.second, rotateCenter);
        }
    }

    void TeethCurveFitter::setPoints(const std::map<int, std::vector<Point_2>>& pointsMap) {
        points.clear();
        for (const auto& pair : pointsMap) {
            points.insert(points.end(), pair.second.begin(), pair.second.end());
        }
    }

    // 计算首尾牙齿的局部凸包
    void TeethCurveFitter::computeSideConvexHull(const std::map<int, std::vector<Point_2>>& pointsMap,
        std::vector<Point_2>& firstTeethHull,
        std::vector<Point_2>& lastTeethHull) {
        auto secondIter = std::next(pointsMap.begin());
        auto penultimateIter = std::prev(pointsMap.end(), 2);
        firstTeethHull = computeUnitConvex(secondIter, returnTeeth::pre);
        lastTeethHull = computeUnitConvex(penultimateIter, returnTeeth::next);
    }

    // 计算凸包主流程
    void TeethCurveFitter::computeConvexHulls() {
        std::vector<Point_2> result;
        std::vector<Point_2> innerTeethHull, outerTeethHull;
        std::vector<Point_2> firstTeethHull, lastTeethHull;
        std::map<int, std::vector<Point_2>> resultPointMap;

        // 计算首尾牙齿的局部凸包
        computeSideConvexHull(pointsMap, firstTeethHull, lastTeethHull);

        // 计算外侧凸包（颊侧）
        result.clear();
        resultPointMap.clear();
        setPoints(pointsMap);
        std::vector<Point_2> tempResult;
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(tempResult));
        OptimizeConvexHull(resultPointMap, tempResult, pointsMap);
        Point_2 firstTeeth, lastTeeth;
        for (auto pair : resultPointMap) {
            result.insert(result.end(), pair.second.begin(), pair.second.end());
            if (pair.first == 1) {
                float sum_x = 0, sum_y = 0;
                for (auto point : pair.second) {
                    sum_x += point.x();
                    sum_y += point.y();
                }
                firstTeeth = Point_2(sum_x / pair.second.size(), sum_y / pair.second.size());
            }
            if (pair.first == static_cast<int>(resultPointMap.size())) {
                float sum_x = 0, sum_y = 0;
                for (auto point : pair.second) {
                    sum_x += point.x();
                    sum_y += point.y();
                }
                lastTeeth = Point_2(sum_x / pair.second.size(), sum_y / pair.second.size());
            }
        }
        outerhullMap = resultPointMap;
        outerTeethHull = result;
        tempResult.clear();

        // 计算内侧凸包（舌侧）
        result.clear();
        resultPointMap.clear();
        setPoints(rotatedPointsMap);
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result));
        OptimizeConvexHull(resultPointMap, result, rotatedPointsMap);
        result.clear();
        // 使用外侧凸包的首尾点覆盖内侧凸包的首尾点
        resultPointMap[1] = std::vector<Point_2>{ firstTeeth };
        resultPointMap[resultPointMap.size()] = std::vector<Point_2>{ lastTeeth };
        // 逐牙旋转回去
        for (auto pair : resultPointMap) {
            Point_2 rotatedCenter = CalculateCentroid(rotatedPointsMap[pair.first]);
            std::vector<Point_2> temp = Rotate(pair.second, rotatedCenter);
            result.insert(result.end(), temp.begin(), temp.end());
        }
        innerhullMap = resultPointMap;
        innerTeethHull = result;
        // 排序凸包点
        std::vector<Point_2> newInnerHull = sortPointsAlongCurve(innerTeethHull, false, false);
        std::vector<Point_2> newOuterHull = sortPointsAlongCurve(outerTeethHull, false, false);
        innerhull = newInnerHull;
        outerhull = newOuterHull;
        // 拟合B样条曲线
        innerBSplinePoints = fitBSpline(newInnerHull, 100);
        outerBSplinePoints = fitBSpline(newOuterHull, 100);
    }
    std::vector<Point_2> TeethCurveFitter::fitBSpline(const std::vector<Point_2>& points, size_t num_samples) {
        if (points.size() < 3) {
            std::cerr << "点数不足，无法进行B样条拟合!" << std::endl;
            return {};
        }

        size_t degree = 3;
        BSpline spline(points.size(), 2, degree);
        std::vector<real> ctrlp;
        for (const auto& point : points) {
            ctrlp.push_back(point.x());
            ctrlp.push_back(point.y());
        }

        spline.setControlPoints(ctrlp);

        std::vector<Point_2> splinePoints;
        for (size_t i = 0; i < num_samples; ++i) {
            real t = i / static_cast<real>(num_samples - 1);
            auto result = spline.eval(t).result();
            splinePoints.emplace_back(result[0], result[1]);
        }
        return splinePoints;
    }
    // 凸包优化
    void TeethCurveFitter::OptimizeConvexHull(std::map<int, std::vector<Point_2>>& resultPointMap,
        const std::vector<Point_2>& result,
        const std::map<int, std::vector<Point_2>>& pointsMap) {
        std::vector<bool> teethWithConvexhull(pointsMap.size(), false);
        teethWithConvexhull[0] = true;
        teethWithConvexhull.back() = true;
        for (Point_2 point : result) {
            for (auto pair : pointsMap) {
                auto it = std::find(pair.second.begin(), pair.second.end(), point);
                if (it != pair.second.end()) {
                    teethWithConvexhull[pair.first - 1] = true;
                    resultPointMap[pair.first].push_back(point);
                    break;
                }
            }
        }
        for (int i = 0; i < teethWithConvexhull.size(); i++) {
            if (!teethWithConvexhull[i]) {
                int j = i, k = i;
                while (!teethWithConvexhull[--j]);
                while (!teethWithConvexhull[++k]);
                auto firsthull_lastpoint = resultPointMap[j + 1].end() - 1;
                auto nexthull_firstpoint = resultPointMap[k + 1].begin();
                Point_2 closePoint = getClosestPoint(pointsMap.at(i + 1), *firsthull_lastpoint, *nexthull_firstpoint);
                resultPointMap[i + 1] = { closePoint };
            }
        }
    }

    // 排序点
    std::vector<Point_2> TeethCurveFitter::sortPointsAlongCurve(const std::vector<Point_2>& points, bool curve_right2left, bool low_high) {
        if (points.empty()) return {};
        std::vector<int> sortedIndices;
        std::vector<Point_2> sortedPoints;
        int startIdx = low_high ? findEndPoint_low(points) : findEndPoint_high(points);
        Point_2 firstPoint = points[startIdx];
        sortedIndices.push_back(startIdx);
        sortedPoints.push_back(firstPoint);
        while (sortedIndices.size() < points.size()) {
            int lastIdx = sortedIndices.back();
            double minDist = std::numeric_limits<double>::max();
            int nextIdx = -1;
            for (size_t i = 0; i < points.size(); ++i) {
                if (std::find(sortedIndices.begin(), sortedIndices.end(), i) == sortedIndices.end()) {
                    double dist = (points[lastIdx].x() - points[i].x()) * (points[lastIdx].x() - points[i].x()) +
                        (points[lastIdx].y() - points[i].y()) * (points[lastIdx].y() - points[i].y());
                    if (dist < minDist) {
                        minDist = dist;
                        nextIdx = i;
                    }
                }
            }
            if (nextIdx == -1) break;
            sortedIndices.push_back(nextIdx);
            sortedPoints.push_back(points[nextIdx]);
        }
        Point_2 lastPoint = sortedPoints.back();
        if ((lastPoint.x() > firstPoint.x() && curve_right2left) || (lastPoint.x() < firstPoint.x() && !curve_right2left))
            std::reverse(sortedPoints.begin(), sortedPoints.end());
        return sortedPoints;
    }



    int TeethCurveFitter::findEndPoint_low(const std::vector<Point_2>& points) {
        double minY = points[0].y();
        int minYIdx = 0;
        for (size_t i = 0; i < points.size(); i++) {
            if (points[i].y() < minY) {
                minYIdx = i;
                minY = points[i].y();
            }
        }
        return minYIdx;
    }

    int TeethCurveFitter::findEndPoint_high(const std::vector<Point_2>& points) {
        double maxY = points[0].y();
        int maxYIdx = 0;
        for (size_t i = 0; i < points.size(); i++) {
            if (points[i].y() > maxY) {
                maxYIdx = i;
                maxY = points[i].y();
            }
        }
        return maxYIdx;
    }

    float TeethCurveFitter::pointToLineDistance(Point_2 point, Point_2 linePoint1, Point_2 linePoint2) {
        double dx = linePoint1.x() - linePoint2.x();
        double dy = linePoint1.y() - linePoint2.y();
        double px = point.x() - linePoint2.x();
        double py = point.y() - linePoint2.y();
        double norm = std::sqrt(dx * dx + dy * dy);
        if (norm == 0) return 0.0f;
        return std::fabs((dx * py - dy * px) / norm);
    }

    Point_2 TeethCurveFitter::getClosestPoint(const std::vector<Point_2>& points, const Point_2& linePoint1, const Point_2& linePoint2) {
        float min_distance = pointToLineDistance(points[0], linePoint1, linePoint2);
        Point_2 closePoint = points[0];
        for (const Point_2& point : points) {
            float distance = pointToLineDistance(point, linePoint1, linePoint2);
            if (min_distance > distance) {
                min_distance = distance;
                closePoint = point;
            }
        }
        return closePoint;
    }

    Point_2 TeethCurveFitter::CalculateCentroid(const std::vector<Point_2>& targetPoints) {
        double sumX = 0, sumY = 0;
        for (const auto& p : targetPoints) {
            sumX += p.x();
            sumY += p.y();
        }
        return Point_2(sumX / targetPoints.size(), sumY / targetPoints.size());
    }

    std::vector<Point_2> TeethCurveFitter::Rotate(const std::vector<Point_2>& targetPoints, Point_2 rotateCenter) {
        std::vector<Point_2> resultPoints;
        for (auto& p : targetPoints) {
            Point_2 temp(-p.x() + 2 * rotateCenter.x(), -p.y() + 2 * rotateCenter.y());
            resultPoints.push_back(temp);
        }
        return resultPoints;
    }

    const std::vector<Point_2>& TeethCurveFitter::getInnerBSpline() const {
        return innerBSplinePoints;
    }
    const std::vector<Point_2>& TeethCurveFitter::getOuterBSpline() const {
        return outerBSplinePoints;
    }

    const std::vector<Point_2>& TeethCurveFitter::getInnerHull() const {
        return innerhull;
    }

    const std::vector<Point_2>& TeethCurveFitter::getOuterHull() const {
        return outerhull;
    }

    const std::map<int, std::vector<Point_2>>& TeethCurveFitter::getInnerHullMap() const {
        return innerhullMap;
    }

    const std::map<int, std::vector<Point_2>>& TeethCurveFitter::getOuterHullMap() const {
        return outerhullMap;
    }
}
