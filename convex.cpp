#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <nlohmann/json.hpp>
#include "tinysplinecxx.h"
#include <sys/stat.h>
#include <direct.h>
#include <unordered_map>
#include"file_utils.h"
#include "json2obj.h"




//// 获取文件名部分（不含扩展名）
//std::string getBaseFileName(const std::string& filePath) {
//    std::string fileName = filePath;
//    size_t pos = fileName.find_last_of("/\\");
//    if (pos != std::string::npos) {
//        fileName = fileName.substr(pos + 1);
//    }
//    pos = fileName.find_last_of(".");
//    if (pos != std::string::npos) {
//        return fileName.substr(0, pos);
//    }
//    return fileName;
//}

enum class returnTeeth {
    all,
    pre,
    current,
    next
};

using json = nlohmann::json;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
using std::vector;
using namespace tinyspline;

// 假设Vector2类已经定义，如果没有需要添加
struct Vector2 {
    double x, y;
    Vector2(double x_, double y_) : x(x_), y(y_) {}
    Vector2 getVertical() const { return Vector2(-y, x); }
    void normalize() {
        double len = std::sqrt(x * x + y * y);
        if (len > 0) { x /= len; y /= len; }
    }
    double operator*(const Vector2& other) const { return x * other.x + y * other.y; }
};

class TeethProcessor {
private:
    // 原始点集
    std::map<int, std::vector<Point_3>> originPointsMap;
    
    // 当前迭代使用的点集
    Point_2* points = nullptr;
    int num_points = 0;
    std::vector<Point_2> innerhull, outerhull;
    std::map<int, std::vector<Point_2>> pointsMap;
    std::map<int, std::vector<Point_2>> rotatedPointsMap;
    std::map<int, std::vector<Point_2>> innerhullMap;
    std::map<int, std::vector<Point_2>> outerhullMap;

    std::vector<Point_2> innerBSplinePoints;
    std::vector<Point_2> outerBSplinePoints;
    std::map<int, std::vector<Point_2>> innerSplineMap;
    std::map<int, std::vector<Point_2>> outerSplineMap;
    bool isRotated = false;
    
    // 每个标签的z值范围
    std::map<int, std::pair<double, double>> zRanges;
    std::map<int, float> min_Zs;

    // 存储所有迭代的样条曲线点和z值
    struct SplineData {
        std::vector<Point_2> innerSpline;
        std::vector<Point_2> outerSpline;
        std::map<int, float> minZ;
        std::unordered_map<std::string, int> innerPointToLabelMap;
        std::unordered_map<std::string, int> outerPointToLabelMap;
    };
    std::vector<SplineData> allSplineData;

    // 初始化原始点集
    void initializeOriginPoints(const std::string& obj_filename) {
        size_t last_dot = obj_filename.find_last_of('.');
        if (last_dot == std::string::npos) {
            throw std::runtime_error("文件名中没有点");
        }
        std::string base_name = obj_filename.substr(0, last_dot);
        std::string json_filename = base_name + ".json";

        if (!fileExists(obj_filename) || !fileExists(json_filename)) {
            std::cerr << "输入文件不存在: " << obj_filename << " 或 " << json_filename << std::endl;
            return;
        }

        // 读取OBJ文件
        std::ifstream obj_file(obj_filename);
        std::string line;
        std::vector<Point_3> allPoints;
        while (std::getline(obj_file, line)) {
            if (line.substr(0, 2) == "v ") {
                std::istringstream iss(line.substr(2));
                double x, y, z;
                iss >> x >> y >> z;
                allPoints.emplace_back(x, y, z);
            }
        }

        // 读取JSON标签文件
        std::ifstream json_file(json_filename);
        json jsonData;
        json_file >> jsonData;
        std::vector<int> labels = jsonData["labels"];

        // 临时存储原始标签的点
        std::map<int, std::vector<Point_3>> tempPointsMap;
        
        // 按标签分组并计算z值范围
        for (size_t i = 0; i < allPoints.size(); ++i) {
            int label = labels[i];
            if (label == 0) continue;

            double z = allPoints[i].z();
            if (zRanges.find(label) == zRanges.end()) {
                zRanges[label] = {z, z};
            } else {
                auto& range = zRanges[label];
                range.first = std::min(range.first, z);
                range.second = std::max(range.second, z);
            }
            
            tempPointsMap[label].push_back(allPoints[i]);
        }

        // 重构标签
        std::map<int, int> keyMap;
        std::vector<int> group_1, group_2;

        for (const auto& pair : tempPointsMap) {
            if (pair.first < 20 || pair.first > 40) group_1.push_back(pair.first);
            else group_2.push_back(pair.first);
        }

        std::sort(group_1.begin(), group_1.end(), std::greater<int>());
        std::sort(group_2.begin(), group_2.end());
        int newkey = 1;
        for (int key : group_1) keyMap[key] = newkey++;
        for (int key : group_2) keyMap[key] = newkey++;

        // 重构zRanges
        std::map<int, std::pair<double, double>> newZRanges;
        for (const auto& pair : zRanges) {
            if (keyMap.find(pair.first) != keyMap.end()) {
                newZRanges[keyMap[pair.first]] = pair.second;
            }
        }
        zRanges = newZRanges;

        // 重构originPointsMap
        for (const auto& pair : tempPointsMap) {
            if (keyMap.find(pair.first) != keyMap.end()) {
                originPointsMap[keyMap[pair.first]] = pair.second;
            }
        }
    }

    // 根据percent筛选点
    void filterPointsByZRange(float percent) {
        pointsMap.clear();
        rotatedPointsMap.clear();
        min_Zs.clear();

        for (const auto& pair : originPointsMap) {
            int label = pair.first;
            const auto& range = zRanges[label];
            double min_z = range.first;
            double max_z = range.second;
            double threshold = max_z - percent * (max_z - min_z);
            // 记录调整后的min_z
            min_Zs[label] = static_cast<float>(min_z + (1.0 - percent) * (max_z - min_z));
            // 筛选点并转换为2D
            for (const auto& point : pair.second) {
                if (point.z() >= threshold && point.z() <= max_z) {
                    pointsMap[label].emplace_back(point.x(), point.y());
                }
            }
        }
        // 关键：每轮都用本轮的pointsMap[label]，以其质心为中心旋转
        for (const auto& pair : pointsMap) {
            int label = pair.first;
            const auto& pts = pair.second;
            if (!pts.empty()) {
                Point_2 rotateCenter = CalculateCentroid(pts);
                rotatedPointsMap[label] = Rotate(pts, rotateCenter);
            }
        }
    }

    // 计算3D点的质心
    Point_2 CalculateCentroid3D(const std::vector<Point_3>& points) {
        double sumX = 0, sumY = 0;
        for (const auto& p : points) {
            sumX += p.x();
            sumY += p.y();
        }
        return Point_2(sumX / points.size(), sumY / points.size());
    }

    // 3D点旋转
    std::vector<Point_3> Rotate3D(const std::vector<Point_3>& points, Point_2 center) {
        std::vector<Point_3> result;
        for (const auto& p : points) {
            result.emplace_back(-p.x() + 2 * center.x(), -p.y() + 2 * center.y(), p.z());
        }
        return result;
    }

    // 设置点集
    void setPoints(const std::map<int, std::vector<Point_2>>& pointsMap) {
        delete[] points;
        num_points = 0;
        for (const auto& pair : pointsMap) {
            num_points += pair.second.size();
        }
        points = new Point_2[num_points];
        int index = 0;
        for (const auto& pair : pointsMap) {
            for (const auto& point : pair.second) {
                points[index++] = point;
            }
        }
    }

    // 加载输入文件
    void loadPoints(const std::string& obj_filename, float percent) {
        size_t last_dot = obj_filename.find_last_of('.');
        if (last_dot == std::string::npos) {
            throw std::runtime_error("文件名中没有点");
        }
        std::string base_name = obj_filename.substr(0, last_dot);
        std::string json_filename = base_name + ".json";

        if (!fileExists(obj_filename) || !fileExists(json_filename)) {
            std::cerr << "输入文件不存在: " << obj_filename << " 或 " << json_filename << std::endl;
            return;
        }

        pointsMap.clear();
        rotatedPointsMap.clear();
        innerhull.clear();
        outerhull.clear();

        read_obj_and_json_files(obj_filename, json_filename, pointsMap, percent);

        GetRotatePoints();
    }

    // B样条拟合
    std::vector<Point_2> fitBSpline(const std::vector<Point_2>& points, size_t num_samples = 100) {
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

    // 优化凸包
    void OptimizeConvexHull(std::map<int, vector<Point_2>>& resultPointMap, vector<Point_2>& result, std::map<int, vector<Point_2>>& pointsMap) {
        vector<bool> teethWithConvexhull(pointsMap.size(), false);
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
                Point_2 closePoint = getClosestPoint(pointsMap[i + 1], *firsthull_lastpoint, *nexthull_firstpoint);
                resultPointMap[i + 1] = { closePoint };
            }
        }
    }

    // 计算凸包
    void computeConvexHull() {
        std::vector<Point_2> result;
        std::vector<Point_2> innerTeethHull, outerTeethHull;
        std::vector<Point_2> firstTeethHull, lastTeethHull;
        std::map<int, std::vector<Point_2>> resultPointMap;

        // 计算首尾牙齿的局部凸包
        computeSideConvexHull(pointsMap, firstTeethHull, lastTeethHull);

        // 计算外侧凸包（颊侧）
        result.clear();
        resultPointMap.clear();
        vector<Point_2> tempResult;
        setPoints(pointsMap);
        CGAL::convex_hull_2(points, points + num_points, std::back_inserter(tempResult));
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
        CGAL::convex_hull_2(points, points + num_points, std::back_inserter(result));
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
        std::cout << "newInnerHull.size():" << innerTeethHull.size() << std::endl;
        std::cout << "newOuterHull.size():" << outerTeethHull.size() << std::endl;
        // 排序凸包点
        std::vector<Point_2> newInnerHull = sortPointsAlongCurve(innerTeethHull, false, false);
        std::vector<Point_2> newOuterHull = sortPointsAlongCurve(outerTeethHull, false, false);
        std::cout << "newInnerHull.size():" << newInnerHull.size() << std::endl;
        std::cout << "newOuterHull.size():" << newOuterHull.size() << std::endl;
        innerhull = newInnerHull;
        outerhull = newOuterHull;
        // 拟合B样条曲线
        innerBSplinePoints = fitBSpline(newInnerHull, 100);
        outerBSplinePoints = fitBSpline(newOuterHull, 100);

        // 保存数据

    }

    // 重置处理器状态
    

    // 保存数据到文件
    void saveData(const std::string& path, int iteration = 0) {
        json j;

        // 如果文件已存在且不是第一次迭代，则读取现有文件
        if (iteration > 0 && fileExists(path)) {
            std::ifstream existing_file(path);
            if (existing_file.is_open()) {
                existing_file >> j;
                existing_file.close();
            }
        }

        // 只在第一次迭代时保存pointsMap、hull_in和hull_out
        if (iteration == 0) {
            for (const auto& pair : pointsMap) {
                json pointsArray = json::array();
                for (const auto& point : pair.second) {
                    pointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
                }
                j["pointsMap"][std::to_string(pair.first)] = pointsArray;
            }

            json hullArray_in = json::array();
            for (const auto& point : innerhull) {
                hullArray_in.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            j["hull_in"] = hullArray_in;

            json hullArray_out = json::array();
            for (const auto& point : outerhull) {
                hullArray_out.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            j["hull_out"] = hullArray_out;

            // 保存innerhullMap和outerhullMap
            for (const auto& pair : innerhullMap) {
                json hullPointsArray = json::array();
                for (const auto& point : pair.second) {
                    hullPointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
                }
                j["innerhullMap"][std::to_string(pair.first)] = hullPointsArray;
            }

            for (const auto& pair : outerhullMap) {
                json hullPointsArray = json::array();
                for (const auto& point : pair.second) {
                    hullPointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
                }
                j["outerhullMap"][std::to_string(pair.first)] = hullPointsArray;
            }
        }

        // 创建点坐标到标签的映射
        auto innerPointToLabelMap = createPointToLabelMap(innerSplineMap);
        auto outerPointToLabelMap = createPointToLabelMap(outerSplineMap);

        // 每次迭代都保存spline_in和spline_out，使用不同的键名
        std::string spline_in_key = iteration == 0 ? "spline_in" : "spline_in_" + std::to_string(iteration);
        std::string spline_out_key = iteration == 0 ? "spline_out" : "spline_out_" + std::to_string(iteration);

        // 为spline_in添加标签
        json splineArray_in = json::array();
        for (const auto& point : innerBSplinePoints) {
            // 使用点坐标作为键查找标签
            std::string key = std::to_string(point.x()) + "," + std::to_string(point.y());
            int label = -1;

            if (innerPointToLabelMap.find(key) != innerPointToLabelMap.end()) {
                label = innerPointToLabelMap[key];
            }

            // 使用标签作为第三个坐标
            splineArray_in.push_back({ {"x", point.x()}, {"y", point.y()}, {"label", label} });
        }
        j[spline_in_key] = splineArray_in;

        // 为spline_out添加标签
        json splineArray_out = json::array();
        for (const auto& point : outerBSplinePoints) {
            // 使用点坐标作为键查找标签
            std::string key = std::to_string(point.x()) + "," + std::to_string(point.y());
            int label = -1;

            if (outerPointToLabelMap.find(key) != outerPointToLabelMap.end()) {
                label = outerPointToLabelMap[key];
            }

            // 使用标签作为第三个坐标
            splineArray_out.push_back({ {"x", point.x()}, {"y", point.y()}, {"label", label} });
        }
        j[spline_out_key] = splineArray_out;

        // 保存innerSplineMap和outerSplineMap
        std::string innerSplineMap_key = iteration == 0 ? "innerSplineMap" : "innerSplineMap_" + std::to_string(iteration);
        std::string outerSplineMap_key = iteration == 0 ? "outerSplineMap" : "outerSplineMap_" + std::to_string(iteration);

        json innerSplineMapJson = json::object();
        for (const auto& pair : innerSplineMap) {
            json splinePointsArray = json::array();
            for (const auto& point : pair.second) {
                splinePointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            innerSplineMapJson[std::to_string(pair.first)] = splinePointsArray;
        }
        j[innerSplineMap_key] = innerSplineMapJson;

        json outerSplineMapJson = json::object();
        for (const auto& pair : outerSplineMap) {
            json splinePointsArray = json::array();
            for (const auto& point : pair.second) {
                splinePointsArray.push_back({ {"x", point.x()}, {"y", point.y()} });
            }
            outerSplineMapJson[std::to_string(pair.first)] = splinePointsArray;
        }
        j[outerSplineMap_key] = outerSplineMapJson;

        // 保存min_Zs数据，使用不同的键名
        std::string min_z_key = iteration == 0 ? "min_z" : "min_z_" + std::to_string(iteration);
        json minZArray = json::object();
        for (const auto& pair : min_Zs) {
            minZArray[std::to_string(pair.first)] = pair.second;
        }
        j[min_z_key] = minZArray;

        std::ofstream file(path);
        if (!file.is_open()) {
            std::cerr << "无法创建/打开文件: " << path << std::endl;
            return;
        }
        file << j.dump(4);
        file.close();
        std::cout << "数据已成功保存至 " << path << " (迭代 " << iteration << ")" << std::endl;
    }

    // 读取OBJ和JSON文件
    void read_obj_and_json_files(
        const std::string& obj_filename,
        const std::string& json_filename,
        std::map<int, std::vector<Point_2>>& TargetpointsMap,
        float percent)
    {
        std::map<int, std::vector<Point_2>> pointsMap;
        std::vector<Point_3> allPoints;
        min_Zs.clear(); // 清空之前的min_Zs数据

        // 读取OBJ文件，包含z值
        std::ifstream obj_file(obj_filename);
        std::string line;
        while (std::getline(obj_file, line)) {
            if (line.substr(0, 2) == "v ") {
                std::istringstream iss(line.substr(2));
                double x, y, z;
                iss >> x >> y >> z;
                allPoints.emplace_back(x, y, z);
            }
        }
        std::cout << json_filename << std::endl;
        // 读取JSON标签文件
        std::ifstream json_file(json_filename);
        json jsonData;
        json_file >> jsonData;
        std::vector<int> labels = jsonData["labels"];

        // 按标签分组并计算每个标签的z值范围
        std::map<int, std::pair<double, double>> labelZRanges; // 存储每个标签的min_z和max_z

        // 首先计算每个标签的z值范围，排除标签为0的点
        for (size_t i = 0; i < allPoints.size(); ++i) {
            int label = labels[i];

            // 跳过标签为0的点
            if (label == 0) {
                continue;
            }

            double z = allPoints[i].z();

            // 更新该标签的z值范围
            if (labelZRanges.find(label) == labelZRanges.end()) {
                labelZRanges[label] = { z, z };
            }
            else {
                auto& range = labelZRanges[label];
                range.first = std::min(range.first, z);  // min_z
                range.second = std::max(range.second, z); // max_z
            }
        }

        // 记录每个标签的最小Z值到min_Zs
        for (const auto& pair : labelZRanges) {
            // 根据percent调整min_z值
            double min_z = pair.second.first;
            double max_z = pair.second.second;
            double adjusted_min_z = min_z + (1.0 - percent) * (max_z - min_z);
            min_Zs[pair.first] = static_cast<float>(adjusted_min_z);
        }

        // 筛选点，只保留z值在指定范围内的点
        for (size_t i = 0; i < allPoints.size(); ++i) {
            int label = labels[i];

            // 跳过标签为0的点
            if (label == 0) {
                continue;
            }

            double z = allPoints[i].z();

            auto& range = labelZRanges[label];
            double min_z = range.first;
            double max_z = range.second;
            double threshold = max_z - percent * (max_z - min_z);

            if (z >= threshold && z <= max_z) {
                pointsMap[label].emplace_back(Point_2(allPoints[i].x(), allPoints[i].y()));
            }
        }
        rebuildPointsMap(pointsMap, TargetpointsMap);
        return;
    }

    // 重构PointsMap键值
    void rebuildPointsMap(const std::map<int, std::vector<Point_2>>& pointsMap, std::map<int, std::vector<Point_2>>& TargetpointsMap) {

        std::map<int, int> keyMap;
        std::vector<int> group_1, group_2;

        for (const auto& pair : pointsMap) {
            if (pair.first < 20 || pair.first > 40) group_1.push_back(pair.first);
            else group_2.push_back(pair.first);
        }

        std::sort(group_1.begin(), group_1.end(), std::greater<int>());
        std::sort(group_2.begin(), group_2.end());
        int newkey = 1;
        for (int key : group_1) keyMap[key] = newkey++;
        for (int key : group_2) keyMap[key] = newkey++;
        for (const auto& pair : pointsMap) TargetpointsMap[keyMap[pair.first]] = pair.second;

        // 重构min_Zs的标签
        std::map<int, float> newMinZs;
        for (const auto& pair : min_Zs) {
            if (keyMap.find(pair.first) != keyMap.end()) {
                newMinZs[keyMap[pair.first]] = pair.second;
            }
        }
        min_Zs = newMinZs;

        return;
    }

    // 计算旋转点
    void GetRotatePoints() {
        for (auto pair : pointsMap) {
            Point_2 rotateCenter = CalculateCentroid(pair.second);
            rotatedPointsMap[pair.first] = Rotate(pair.second, rotateCenter);
        }
    }

    // 计算质心
    Point_2 CalculateCentroid(const vector<Point_2>& targetPoints) {
        double sumX = 0, sumY = 0;
        for (const auto& p : targetPoints) {
            sumX += p.x();
            sumY += p.y();
        }
        return Point_2(sumX / targetPoints.size(), sumY / targetPoints.size());
    }

    // 旋转180度
    vector<Point_2> Rotate(const vector<Point_2>& targetPoints, Point_2 rotateCenter) {
        vector<Point_2> resultPoints;
        for (auto& p : targetPoints) {
            Point_2 temp(-p.x() + 2 * rotateCenter.x(), -p.y() + 2 * rotateCenter.y());
            resultPoints.push_back(temp);
        }
        return resultPoints;
    }

    // 计算两端牙齿的凸包
    void computeSideConvexHull(std::map<int, std::vector<Point_2>>& pointsMap, vector<Point_2>& firstTeethHull, vector<Point_2>& lastTeethHull) {
        auto secondIter = std::next(pointsMap.begin());
        auto penultimateIter = std::prev(pointsMap.end(), 2);
        firstTeethHull = computeUnitConvex(secondIter, returnTeeth::pre);
        lastTeethHull = computeUnitConvex(penultimateIter, returnTeeth::next);
    }

    // 计算单元凸包
    template<typename Iterator>
    std::vector<Point_2> computeUnitConvex(Iterator it, returnTeeth op) {
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
        case returnTeeth::pre:
            for (auto point : temp) {
                if (std::find(prevIt->second.begin(), prevIt->second.end(), point) != prevIt->second.end())
                    TeethHull.push_back(point);
            }
            break;
        case returnTeeth::next:
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

    // 点到直线的距离
    float pointToLineDistance(Point_2 point, Point_2 linePoint1, Point_2 linePoint2) {
        Vector2 line(linePoint1.x() - linePoint2.x(), linePoint1.y() - linePoint2.y());
        Vector2 pointToLine(point.x() - linePoint2.x(), point.y() - linePoint2.y());
        Vector2 linenormal = line.getVertical();
        linenormal.normalize();
        return fabs(linenormal * pointToLine);
    }

    // 获取最近点
    Point_2 getClosestPoint(const vector<Point_2>& points, const Point_2& linePoint1, const Point_2& linePoint2) {
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

    // 沿曲线排序点
    vector<Point_2> sortPointsAlongCurve(const std::vector<Point_2>& points, bool curve_right2left, bool low_high = true) {
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

    // 找到Y坐标最小点索引
    int findEndPoint_low(const vector<Point_2>& points) {
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

    // 找到Y坐标最大点索引
    int findEndPoint_high(const vector<Point_2>& points) {
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



    // 为B样条曲线上的点分配标签
    void assignSplineLabels() {
        // 清空之前的标签映射
        innerSplineMap.clear();
        outerSplineMap.clear();

        // 为内侧B样条曲线上的点分配标签
        assignLabelsToSplinePoints(innerBSplinePoints, innerhullMap, innerSplineMap);

        // 为外侧B样条曲线上的点分配标签
        assignLabelsToSplinePoints(outerBSplinePoints, outerhullMap, outerSplineMap);
    }

    // 为B样条曲线上的点分配标签的辅助函数
    void assignLabelsToSplinePoints(const std::vector<Point_2>& splinePoints,
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

    // 创建点坐标到标签的映射
    std::unordered_map<std::string, int> createPointToLabelMap(const std::map<int, std::vector<Point_2>>& splineMap) {
        std::unordered_map<std::string, int> pointToLabelMap;
        for (const auto& pair : splineMap) {
            int label = pair.first;
            for (const auto& point : pair.second) {
                std::string key = std::to_string(point.x()) + "," + std::to_string(point.y());
                pointToLabelMap[key] = label;
            }
        }
        return pointToLabelMap;
    }

    // 直接生成.obj文件
    void generateObjFile(const std::string& outputPath) {
        std::ofstream file(outputPath);
        if (!file.is_open()) {
            std::cerr << "无法创建文件: " << outputPath << std::endl;
            return;
        }

        int vertexCount = 0;
        std::vector<std::vector<Point_3>> vertexGroups;

        // 处理每一组样条曲线数据
        for (size_t i = 0; i < allSplineData.size(); ++i) {
            auto& data = allSplineData[i];
            std::vector<Point_3> vertices;

            // 处理内侧样条曲线
            for (const auto& point : data.innerSpline) {
                std::string key = std::to_string(point.x()) + "," + std::to_string(point.y());
                auto it = data.innerPointToLabelMap.find(key);
                if (it != data.innerPointToLabelMap.end()) {
                    int label = it->second;
                    float z1 = data.minZ[label];
                    float z2 = (i < allSplineData.size() - 1) ? allSplineData[i + 1].minZ[label] : z1 + 5.0f;
                    vertices.emplace_back(point.x(), point.y(), z1);
                    vertices.emplace_back(point.x(), point.y(), z2);
                }
            }
            vertexGroups.push_back(vertices);

            // 处理外侧样条曲线
            vertices.clear();
            for (const auto& point : data.outerSpline) {
                std::string key = std::to_string(point.x()) + "," + std::to_string(point.y());
                auto it = data.outerPointToLabelMap.find(key);
                if (it != data.outerPointToLabelMap.end()) {
                    int label = it->second;
                    float z1 = data.minZ[label];
                    float z2 = (i < allSplineData.size() - 1) ? allSplineData[i + 1].minZ[label] : z1 + 5.0f;
                    vertices.emplace_back(point.x(), point.y(), z1);
                    vertices.emplace_back(point.x(), point.y(), z2);
                }
            }
            vertexGroups.push_back(vertices);
        }

        // 写入所有顶点
        for (const auto& group : vertexGroups) {
            for (const auto& v : group) {
                file << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
                vertexCount++;
            }
        }

        // 写入线段
        std::vector<int> offsets = { 0 };
        for (size_t i = 0; i < vertexGroups.size() - 1; ++i) {
            offsets.push_back(offsets.back() + vertexGroups[i].size());
        }

        for (size_t i = 0; i < vertexGroups.size(); ++i) {
            const auto& group = vertexGroups[i];
            int groupSize = group.size() / 2;
            int offset = offsets[i];

            for (int j = 0; j < groupSize; ++j) {
                // 垂直线段
                int v1 = offset + j * 2 + 1;
                int v2 = offset + j * 2 + 2;
                file << "l " << v1 << " " << v2 << std::endl;

                // 水平线段
                if (j < groupSize - 1) {
                    v2 = offset + j * 2 + 2;
                    int v4 = offset + (j + 1) * 2 + 2;
                    file << "l " << v2 << " " << v4 << std::endl;
                }
            }
        }

        std::cout << "已生成.obj文件: " << outputPath << std::endl;
        std::cout << "包含 " << vertexCount << " 个顶点" << std::endl;
    }

   

public:
    TeethProcessor() = default;
    ~TeethProcessor() { delete[] points; }
    // 主处理函数
    void process(const std::string& input_obj, const std::string& output_json, float percent, int iteration = 0) {
        // 如果是第一次迭代，初始化原始点集
        if (iteration == 0) {
            initializeOriginPoints(input_obj);
        }
        
        // 每次迭代都根据percent筛选点
        filterPointsByZRange(percent);

        // 检查筛选后的点集是否为空
        if (pointsMap.empty()) {
            std::cerr << "筛选后的点集为空，请检查percent值或输入数据。" << std::endl;
            return;
        }
        
        // 计算凸包
        computeConvexHull();

        // 为B样条曲线上的点分配标签
        assignSplineLabels();

        // 保存数据
        saveData(output_json, iteration);

        // 在每次迭代后存储样条曲线数据
        SplineData data;
        data.innerSpline = innerBSplinePoints;
        data.outerSpline = outerBSplinePoints;
        data.minZ = min_Zs;
        data.innerPointToLabelMap = createPointToLabelMap(innerSplineMap);
        data.outerPointToLabelMap = createPointToLabelMap(outerSplineMap);
        allSplineData.push_back(data);

        // 如果是最后一次迭代，生成.obj文件
        // if (iteration == 2) {  // 假设总共有3次迭代
        //     std::cout << allSplineData.size() << std::endl;
        //     std::cout << allSplineData[0].innerSpline.size() << std::endl;
        //     std::cout << allSplineData[1].innerSpline.size() << std::endl;
        //     std::cout << allSplineData[2].innerSpline.size() << std::endl;
        //     std::string output_obj = output_json.substr(0, output_json.find_last_of('.')) + ".obj";
        //     generateObjFile(output_obj);
        // }
    }
    void reset() {
        delete[] points;
        points = nullptr;
        num_points = 0;
        innerhull.clear();
        innerhullMap.clear();
        outerhull.clear();
        outerhullMap.clear();
        pointsMap.clear();
        rotatedPointsMap.clear();
        innerBSplinePoints.clear();
        outerBSplinePoints.clear();
        innerSplineMap.clear();
        outerSplineMap.clear();
        isRotated = false;
        min_Zs.clear(); // 重置min_Zs
       
    }
};

// 处理牙齿数据的函数接口（带默认percent值）


// 处理牙齿数据的函数接口（带三个percent值）
bool processTeethData(const std::string& input_obj, float percent1, float percent2, float percent3) {
    // 检查输入文件是否存在
    if (!fileExists(input_obj)) {
        std::cerr << "输入文件不存在: " << input_obj << std::endl;
        return false;
    }

    // 验证percent值
    if (percent1 <= 0.0f || percent1 > 1.0f ||
        percent2 <= 0.0f || percent2 > 1.0f ||
        percent3 <= 0.0f || percent3 > 1.0f) {
        std::cerr << "错误：percent值必须在0到1之间" << std::endl;
        return false;
    }

    // 验证percent值是否递减
    if (percent1 < percent2 || percent2 < percent3) {
        std::cerr << "错误：percent值必须递减 (percent1 >= percent2 >= percent3)" << std::endl;
        return false;
    }

    // 获取输入文件的目录和基本文件名
    std::string input_dir = getDirectory(input_obj);
    std::string base_name = getBaseFileName(input_obj);
    std::string json_path = input_dir + "/" + base_name + ".json";

    // 设置输出JSON文件路径
    std::string output_json = input_dir + "/" + base_name + "_output.json";

    // 创建TeethProcessor对象
    TeethProcessor processor;

    // 处理每个percent值
    std::vector<float> percents = { percent1, percent2, percent3 };
    for (int i = 0; i < percents.size(); i++) {
        std::cout << "处理 percent = " << percents[i] << " (迭代 " << i << ")" << std::endl;

        // 处理当前percent值
        processor.process(input_obj, output_json, percents[i], i);

        // 如果不是最后一次迭代，重置处理器状态
        if (i < percents.size() - 1) {
            processor.reset();
        }
    }

    // 生成.obj文件
    std::string output_obj = input_dir + "/" + base_name + "_output.obj";
    std::cout << "正在生成.obj文件: " << output_obj << std::endl;

    if (jsonToObj(output_json, output_obj)) {
        std::cout << "成功生成.obj文件" << std::endl;
        return true;
    }
    else {
        std::cerr << "生成.obj文件失败" << std::endl;
        return false;
    }
}
bool processTeethData(const std::string& input_obj) {
    return processTeethData(input_obj, 1.0f, 0.7f, 0.5f);
}