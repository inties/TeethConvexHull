#include "TeethResultExporter.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace TeethConvex {
    void TeethResultExporter::addIterationData(
        const std::map<int, std::vector<Point_2>>& pointsMap,
        const std::map<int, std::pair<double, double>>& zRanges,
        const std::vector<Point_2>& innerHull,
        const std::vector<Point_2>& outerHull,
        const std::map<int, std::vector<Point_2>>& innerHullMap,
        const std::map<int, std::vector<Point_2>>& outerHullMap,
        const std::vector<Point_2>& innerBSpline,
        const std::vector<Point_2>& outerBSpline
    ) {
        // 创建样条曲线的副本用于重新分布
        std::vector<Point_2> redistributedInnerBSpline = innerBSpline;
        std::vector<Point_2> redistributedOuterBSpline = outerBSpline;
        
        // 重新分布样条曲线点，使其均匀分布
        reAssignSpline(redistributedInnerBSpline, redistributedOuterBSpline);

        IterationData data;
        data.pointsMap = pointsMap;
        data.zRanges = zRanges;
        data.innerHull = innerHull;
        data.outerHull = outerHull;
        data.innerHullMap = innerHullMap;
        data.outerHullMap = outerHullMap;
        data.innerBSpline = redistributedInnerBSpline;
        data.outerBSpline = redistributedOuterBSpline;
        //data.innerBSpline = innerBSpline;
        //data.outerBSpline = outerBSpline;
        // 为样条点分配标签
        std::map<int, std::vector<Point_2>> innerSplineMap, outerSplineMap;
        assignLabelsToSplinePoints(redistributedInnerBSpline, innerHullMap, innerSplineMap);
        assignLabelsToSplinePoints(redistributedOuterBSpline, outerHullMap, outerSplineMap);

        // 建立点到标签的映射
        data.inPointLabel = createPointToLabelMap(innerSplineMap);
        data.outPointLabel = createPointToLabelMap(outerSplineMap);

        m_iterations.push_back(std::move(data));
    }

    void TeethResultExporter::exportJSON(const std::string& jsonPath) {
        nlohmann::json json;

        for (size_t i = 0; i < m_iterations.size(); ++i) {
            const auto& iter = m_iterations[i];
            std::string suffix = std::to_string(i);

            // 写入pointsMap
            nlohmann::json pointsMapJson = nlohmann::json::object();
            for (const auto& pair : iter.pointsMap) {
                nlohmann::json pointsArray = nlohmann::json::array();
                for (const auto& point : pair.second) {
                    pointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
                }
                pointsMapJson[std::to_string(pair.first)] = pointsArray;
            }
            json["pointsMap_" + suffix] = pointsMapJson;

            // 写入innerHull
            nlohmann::json innerHullJson = nlohmann::json::array();
            for (const auto& point : iter.innerHull) {
                innerHullJson.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            json["hull_in_" + suffix] = innerHullJson;

            // 写入outerHull
            nlohmann::json outerHullJson = nlohmann::json::array();
            for (const auto& point : iter.outerHull) {
                outerHullJson.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            json["hull_out_" + suffix] = outerHullJson;

            // 写入内侧样条曲线（带标签）
            nlohmann::json innerBSplineJson = nlohmann::json::array();
            for (const auto& point : iter.innerBSpline) {
                std::string key = makeKey(point);
                int label = -1; // 默认标签
                if (iter.inPointLabel.find(key) != iter.inPointLabel.end()) {
                    label = iter.inPointLabel.at(key);
                }
                innerBSplineJson.push_back({
                    {"x", point.x()},
                    {"y", point.y()},
                    {"label", label}
                    });
            }
            json["spline_in_" + suffix] = innerBSplineJson;

            // 写入外侧样条曲线（带标签）
            nlohmann::json outerBSplineJson = nlohmann::json::array();
            for (const auto& point : iter.outerBSpline) {
                std::string key = makeKey(point);
                int label = -1; // 默认标签
                if (iter.outPointLabel.find(key) != iter.outPointLabel.end()) {
                    label = iter.outPointLabel.at(key);
                }
                outerBSplineJson.push_back({
                    {"x", point.x()},
                    {"y", point.y()},
                    {"label", label}
                    });
            }
            json["spline_out_" + suffix] = outerBSplineJson;

            // 写入z范围
            nlohmann::json zRangesJson = nlohmann::json::object();
            for (const auto& pair : iter.zRanges) {
                zRangesJson[std::to_string(pair.first)] = {
                    {"min", pair.second.first},
                    {"max", pair.second.second}
                };
            }
            json["zRanges_" + suffix] = zRangesJson;
        }

        // 写入文件
        std::ofstream file(jsonPath);
        if (!file.is_open()) {
            std::cerr << "无法创建 JSON 文件: " << jsonPath << std::endl;
            return;
        }
        file << json.dump(4); // 4是缩进空格数
        file.close();

        std::cout << "success export JSON file: " << jsonPath << std::endl;
    }

	void TeethResultExporter::exportOBJ_2(const std::string& objPath) {// 生成两个OBJ文件，一个内侧样条曲线，一个外侧样条曲线
        if (m_iterations.empty()) {
            std::cerr << "no iteration data, cannot generate OBJ file" << std::endl;
            return;
        }

		std::string basePath = objPath.substr(0, objPath.find_last_of("."));
		std::string innerObjPath = basePath + "_innerTeethConvex.obj";
		std::string outerObjPath = basePath + "_outerTeethConvex.obj";
        exportSingleSplineOBJ(innerObjPath, true, 1, false);//写入内侧样条曲线
        exportSingleSplineOBJ(outerObjPath, false, 1, false);//写入外侧样条曲线

        std::cout << "success generate OBJ files: " << innerObjPath << "inner.obj and " << outerObjPath << "outer.obj" << std::endl;
    }

	void TeethResultExporter::exportOBJ(const std::string& objPath) {// 生成一个OBJ文件，包含内侧和外侧样条曲线
        if (m_iterations.empty()) {
            std::cerr << "no iteration data, cannot generate OBJ file" << std::endl;
            return;
        }

        
		int startIdx = exportSingleSplineOBJ(objPath, true, 1, false);//写入内侧样条曲线
		exportSingleSplineOBJ(objPath, false, startIdx, true);//追加外侧样条曲线
        
        std::cout << "success generate OBJ files: " << objPath<< std::endl;
    }

    int TeethResultExporter::exportSingleSplineOBJ(const std::string& objPath, bool isInner, int startIdx,bool append) {
        std::ofstream file(objPath, append ? std::ios::app : std::ios::out); // 支持追加
        if (!file.is_open()) {
            std::cerr << "cannot create or open OBJ file: " << objPath << std::endl;
            return -1;
        }

        // 收集所有顶点
        std::vector<std::vector<Point_3>> allIterationVertices;
        int totalVertexCount = 0;

        // 遍历所有迭代生成顶点
        for (size_t i = 0; i < m_iterations.size(); ++i) {
            const auto& curIter = m_iterations[i];
            
            // 选择内侧或外侧样条曲线
            const auto& splinePoints = isInner ? curIter.innerBSpline : curIter.outerBSpline;
            const auto& pointLabelMap = isInner ? curIter.inPointLabel : curIter.outPointLabel;

            std::vector<Point_3> iterationVertices;
            
            for (const auto& point : splinePoints) {
                std::string key = makeKey(point);
                int label = -1;

                // 获取点的标签
                if (pointLabelMap.find(key) != pointLabelMap.end()) {
                    label = pointLabelMap.at(key);
                } else {
                    continue; // 跳过没有标签的点
                }

                // 获取对应标签的最小z值
                double minZ = 0.0;
                if (curIter.zRanges.find(label) != curIter.zRanges.end()) {
                    minZ = curIter.zRanges.at(label).first;
                }

                // 创建3D点
                iterationVertices.emplace_back(point.x(), point.y(), minZ);
            }
            
            allIterationVertices.push_back(iterationVertices);
            totalVertexCount += iterationVertices.size();
        }

        // 写入所有顶点
        for (const auto& iterationVertices : allIterationVertices) {
            for (const auto& vertex : iterationVertices) {
                file << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
            }
        }

        // 计算每次迭代顶点的起始索引（OBJ文件中顶点索引从1开始）
        std::vector<int> iterationStartIndices;
        int currentIndex = startIdx;
        for (const auto& iterationVertices : allIterationVertices) {
            iterationStartIndices.push_back(currentIndex);
            currentIndex += iterationVertices.size();
        }

        // 生成连接关系
        // 1. 同一迭代内的连接（相邻点之间）
        for (size_t i = 0; i < allIterationVertices.size(); ++i) {
            const auto& iterationVertices = allIterationVertices[i];
            int startIndex = iterationStartIndices[i];
            
            for (size_t j = 0; j < iterationVertices.size(); ++j) {
                if (j < iterationVertices.size() - 1) {
                    // 连接到下一个点
                    file << "l " << (startIndex + j) << " " << (startIndex + j + 1) << std::endl;
                }
            }
        }

        // 2. 相邻迭代间的连接（对应位置的点）
        for (size_t i = 0; i < allIterationVertices.size() - 1; ++i) {
            const auto& currentIterVertices = allIterationVertices[i];
            const auto& nextIterVertices = allIterationVertices[i + 1];
            
            int currentStartIndex = iterationStartIndices[i];
            int nextStartIndex = iterationStartIndices[i + 1];
            
            // 连接对应位置的点（取较小的长度以避免越界）
            size_t minSize = std::min(currentIterVertices.size(), nextIterVertices.size());
            for (size_t j = 0; j < minSize; ++j) {
                file << "l " << (currentStartIndex + j) << " " << (nextStartIndex + j) << std::endl;
            }
        }

        file.close();
        std::cout << "generated " << (isInner ? "inner" : "outer") << " spline OBJ with " 
                  << totalVertexCount << " vertices" << std::endl;
		return currentIndex; // 返回总顶点数，供后续处理使用
    }

    std::unordered_map<std::string, int> TeethResultExporter::createPointToLabelMap(
        const std::map<int, std::vector<Point_2>>& splineMap) const {
        std::unordered_map<std::string, int> pointToLabelMap;
        for (const auto& pair : splineMap) {
            int label = pair.first;
            for (const auto& point : pair.second) {
                std::string key = makeKey(point);
                pointToLabelMap[key] = label;
            }
        }
        return pointToLabelMap;
    }

    void TeethResultExporter::assignLabelsToSplinePoints(
        const std::vector<Point_2>& splinePoints,
        const std::map<int, std::vector<Point_2>>& hullMap,
        std::map<int, std::vector<Point_2>>& splineMap) {
        if (splinePoints.empty() || hullMap.empty()) {
            return;
        }

        // 获取所有标签并排序
        std::vector<int> labels;
        for (const auto& pair : hullMap) {
            labels.push_back(pair.first);
        }
        std::sort(labels.begin(), labels.end());

        // 为每个B样条曲线上的点分配标签
        int prevLabel = -1; // 上一个点的标签

        for (size_t i = 0; i < splinePoints.size(); ++i) {
            const Point_2& splinePoint = splinePoints[i];
            int bestLabel = -1;
            double minDistance = std::numeric_limits<double>::max();

            // 确定搜索范围
            std::vector<int> searchLabels;
            if (prevLabel == -1) {
                // 第一个点，搜索所有标签
                searchLabels = labels;
            }
            else {
                // 根据优化策略，只搜索当前标签和下一个标签
                searchLabels.push_back(prevLabel);
                if (prevLabel + 1 <= labels.back()) {
                    searchLabels.push_back(prevLabel + 1);
                }
            }

            // 在确定的标签范围内查找最近点
            for (int label : searchLabels) {
                const auto& hullPoints = hullMap.at(label);
                for (const auto& hullPoint : hullPoints) {
                    double distance = CGAL::squared_distance(splinePoint, hullPoint);
                    if (distance < minDistance) {
                        minDistance = distance;
                        bestLabel = label;
                    }
                }
            }

            // 如果在优化范围内没有找到足够近的点，则搜索所有标签
            if (bestLabel == -1 || minDistance > 1.0) { // 可以调整阈值
                for (int label : labels) {
                    if (std::find(searchLabels.begin(), searchLabels.end(), label) != searchLabels.end()) {
                        continue; // 跳过已经搜索过的标签
                    }

                    const auto& hullPoints = hullMap.at(label);
                    for (const auto& hullPoint : hullPoints) {
                        double distance = CGAL::squared_distance(splinePoint, hullPoint);
                        if (distance < minDistance) {
                            minDistance = distance;
                            bestLabel = label;
                        }
                    }
                }
            }

            // 将点添加到对应的标签组
            if (bestLabel != -1) {
                splineMap[bestLabel].push_back(splinePoint);
                prevLabel = bestLabel;
            }
        }
    }

    std::vector<Point2D> TeethResultExporter::convertToPoint2D(
        const std::vector<Point_2>& splinePoints,
        const std::unordered_map<std::string, int>& pointLabelMap) const {
        std::vector<Point2D> result;
        result.reserve(splinePoints.size());

        for (const auto& point : splinePoints) {
            std::string key = makeKey(point);
            int label = -1;

            auto it = pointLabelMap.find(key);
            if (it != pointLabelMap.end()) {
                label = it->second;
            }

            result.push_back({ point.x(), point.y(), label });
        }

        return result;
    }

    void TeethResultExporter::reAssignSpline(
        std::vector<Point_2>& innerBSpline,
        std::vector<Point_2>& outerBSpline
    ) {
        // 如果样条曲线为空，直接返回
        if (innerBSpline.empty() || outerBSpline.empty()) {
            return;
        }

        // 如果是第一次迭代，计算质心
        if (!m_centroidsComputed) {
            m_innerCentroid = calculateCentroid(innerBSpline);
            m_outerCentroid = calculateCentroid(outerBSpline);
            m_centroidsComputed = true;
        }

        // 重新分布内侧样条曲线点
        std::vector<Point_2> newInnerBSpline = getIntersectionPoints(m_innerCentroid, innerBSpline);
        if (!newInnerBSpline.empty()) {
            innerBSpline = std::move(newInnerBSpline);
        }

        // 重新分布外侧样条曲线点
        std::vector<Point_2> newOuterBSpline = getIntersectionPoints(m_outerCentroid, outerBSpline);
        if (!newOuterBSpline.empty()) {
            outerBSpline = std::move(newOuterBSpline);
        }
    }

    Point_2 TeethResultExporter::calculateCentroid(const std::vector<Point_2>& points) const {
        if (points.empty()) {
            return Point_2(0, 0);
        }

        double sumX = 0.0, sumY = 0.0;
        for (const auto& point : points) {
            sumX += point.x();
            sumY += point.y();
        }

        return Point_2(sumX / points.size(), sumY / points.size());
    }

    std::vector<Point_2> TeethResultExporter::getIntersectionPoints(
        const Point_2& centroid,
        const std::vector<Point_2>& splinePoints,
        double angleStep
    ) const {
        std::vector<Point_2> intersectionPoints;
        
        // 从0度到360度，每angleStep度发射一条射线
        for (double angle = 90.0; angle < 360.0+90.0; angle += angleStep) {
            Point_2 intersection;
            if (raySplineIntersection(centroid, angle, splinePoints, intersection)) {
                intersectionPoints.push_back(intersection);
            }
        }

        return intersectionPoints;
    }

    bool TeethResultExporter::raySplineIntersection(
        const Point_2& rayOrigin,
        double rayAngle,
        const std::vector<Point_2>& splinePoints,
        Point_2& intersection
    ) const {
        if (splinePoints.size() < 2) {
            return false;
        }

        double minDistance = std::numeric_limits<double>::max();
        Point_2 closestIntersection;
        bool found = false;

        // 遍历样条曲线的所有线段
        for (size_t i = 0; i < splinePoints.size() - 1; ++i) {
            Point_2 segmentIntersection;
            if (raySegmentIntersection(rayOrigin, rayAngle, splinePoints[i], splinePoints[i + 1], segmentIntersection)) {
                // 计算交点到射线起点的距离
                // double distance = CGAL::squared_distance(rayOrigin, segmentIntersection);
                double distance = pow(rayOrigin.x() - segmentIntersection.x(), 2) + pow(rayOrigin.y() - segmentIntersection.y(), 2);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestIntersection = segmentIntersection;
                    found = true;
                }
            }
        }

        if (found) {
            intersection = closestIntersection;
        }
        return found;
    }

    bool TeethResultExporter::raySegmentIntersection(
        const Point_2& rayOrigin,
        double rayAngle,
        const Point_2& segStart,
        const Point_2& segEnd,
        Point_2& intersection
    ) const {
        // 将角度转换为方向向量
        auto direction = angleToDirection(rayAngle);
        double rayDx = direction.first;
        double rayDy = direction.second;

        // 射线参数方程: P = rayOrigin + t * (rayDx, rayDy), t >= 0
        // 线段参数方程: Q = segStart + s * (segEnd - segStart), 0 <= s <= 1

        double segDx = segEnd.x() - segStart.x();
        double segDy = segEnd.y() - segStart.y();

        // 求解线性方程组
        // rayOrigin.x + t * rayDx = segStart.x + s * segDx
        // rayOrigin.y + t * rayDy = segStart.y + s * segDy

        double denominator = rayDx * segDy - rayDy * segDx;
        
        // 如果分母为0，说明射线和线段平行
        if (std::abs(denominator) < 1e-10) {
            return false;
        }

        double dx = segStart.x() - rayOrigin.x();
        double dy = segStart.y() - rayOrigin.y();

        double t = (dx * segDy - dy * segDx) / denominator;
        double s = (dx * rayDy - dy * rayDx) / denominator;

        // 检查参数范围
        if (t >= 0 && s >= 0 && s <= 1) {
            // 计算交点
            intersection = Point_2(
                rayOrigin.x() + t * rayDx,
                rayOrigin.y() + t * rayDy
            );
            return true;
        }

        return false;
    }

    std::pair<double, double> TeethResultExporter::angleToDirection(double angleInDegrees) const {
        double radians = angleInDegrees * M_PI / 180.0;
        return std::make_pair(std::cos(radians), std::sin(radians));
    }
}