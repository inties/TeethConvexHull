#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>
#include <sys/stat.h>
#include <direct.h>
#include"json2obj.h"
using json = nlohmann::json;









// 从JSON文件加载数据
void loadData(const std::string& filePath,
    std::map<int, std::vector<Point2D>>& pointsMap,
    std::vector<Point2D>& hullIn,
    std::vector<Point2D>& hullOut,
    std::vector<Point2D>& splineIn,
    std::vector<Point2D>& splineOut,
    std::vector<Point2D>& splineIn1,
    std::vector<Point2D>& splineOut1,
    std::vector<Point2D>& splineIn2,
    std::vector<Point2D>& splineOut2,
    std::map<int, double>& minZ,
    std::map<int, double>& minZ1,
    std::map<int, double>& minZ2) {

    // 打开JSON文件
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filePath << std::endl;
        return;
    }

    // 解析JSON
    json data;
    file >> data;

    // 读取pointsMap
    if (data.contains("pointsMap")) {
        for (const auto& item : data["pointsMap"].items()) {
            int label = std::stoi(item.key());
            std::vector<Point2D> points;
            for (const auto& point : item.value()) {
                points.emplace_back(point["x"], point["y"]);
            }
            pointsMap[label] = points;
        }
    }

    // 读取hull_in
    if (data.contains("hull_in")) {
        for (const auto& point : data["hull_in"]) {
            hullIn.emplace_back(point["x"], point["y"]);
        }
    }

    // 读取hull_out
    if (data.contains("hull_out")) {
        for (const auto& point : data["hull_out"]) {
            hullOut.emplace_back(point["x"], point["y"]);
        }
    }

    // 读取spline_in
    if (data.contains("spline_in")) {
        for (const auto& point : data["spline_in"]) {
            splineIn.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取spline_out
    if (data.contains("spline_out")) {
        for (const auto& point : data["spline_out"]) {
            splineOut.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取spline_in_1
    if (data.contains("spline_in_1")) {
        for (const auto& point : data["spline_in_1"]) {
            splineIn1.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取spline_out_1
    if (data.contains("spline_out_1")) {
        for (const auto& point : data["spline_out_1"]) {
            splineOut1.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取spline_in_2
    if (data.contains("spline_in_2")) {
        for (const auto& point : data["spline_in_2"]) {
            splineIn2.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取spline_out_2
    if (data.contains("spline_out_2")) {
        for (const auto& point : data["spline_out_2"]) {
            splineOut2.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // 读取min_z
    if (data.contains("min_z")) {
        for (const auto& item : data["min_z"].items()) {
            minZ[std::stoi(item.key())] = item.value();
        }
    }

    // 读取min_z_1
    if (data.contains("min_z_1")) {
        for (const auto& item : data["min_z_1"].items()) {
            minZ1[std::stoi(item.key())] = item.value();
        }
    }

    // 读取min_z_2
    if (data.contains("min_z_2")) {
        for (const auto& item : data["min_z_2"].items()) {
            minZ2[std::stoi(item.key())] = item.value();
        }
    }

    std::cout << "成功加载数据" << std::endl;
}

// 生成.obj文件
void generateObj(const std::vector<Point2D>& splineIn,
    const std::vector<Point2D>& splineOut,
    const std::vector<Point2D>& splineIn1,
    const std::vector<Point2D>& splineOut1,
    const std::vector<Point2D>& splineIn2,
    const std::vector<Point2D>& splineOut2,
    const std::map<int, double>& minZ,
    const std::map<int, double>& minZ1,
    const std::map<int, double>& minZ2,
    const std::string& outputPath,
    double zMin,
    double zMax) {

    std::ofstream file(outputPath);
    if (!file.is_open()) {
        std::cerr << "无法创建文件: " << outputPath << std::endl;
        return;
    }

    int vertexCount = 0;
    std::vector<std::vector<Point3D>> vertexGroups;

    // 处理spline_in (使用minZ和minZ1)
    std::vector<Point3D> vertices;
    for (const auto& point : splineIn) {
        int label = point.label;
        double z1 = (minZ.find(label) != minZ.end()) ? minZ.at(label) : zMin;
        double z2 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 处理spline_out (使用minZ和minZ1)
    vertices.clear();
    for (const auto& point : splineOut) {
        int label = point.label;
        double z1 = (minZ.find(label) != minZ.end()) ? minZ.at(label) : zMin;
        double z2 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 处理spline_in_1 (使用minZ1和minZ2)
    vertices.clear();
    for (const auto& point : splineIn1) {
        int label = point.label;
        double z1 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMin;
        double z2 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 处理spline_out_1 (使用minZ1和minZ2)
    vertices.clear();
    for (const auto& point : splineOut1) {
        int label = point.label;
        double z1 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMin;
        double z2 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 处理spline_in_2 (使用minZ2和zMax)
    vertices.clear();
    for (const auto& point : splineIn2) {
        int label = point.label;
        double z1 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMin;
        double z2 = z1 + 5.0;  // 顶部(没有更高层)
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 处理spline_out_2 (使用minZ2和zMax)
    vertices.clear();
    for (const auto& point : splineOut2) {
        int label = point.label;
        double z1 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMin;
        double z2 = z1 + 5.0;  // 顶部(没有更高层)
        vertices.emplace_back(point.x, point.y, z1);  // 底部
        vertices.emplace_back(point.x, point.y, z2);  // 顶部
    }
    vertexGroups.push_back(vertices);

    // 写入所有顶点
    for (const auto& group : vertexGroups) {
        for (const auto& v : group) {
            file << "v " << v.x << " " << v.y << " " << v.z << std::endl;
            vertexCount++;
        }
    }

    // 写入线段 (l)
    // 计算每组顶点的偏移量
    std::vector<int> offsets = { 0 };
    for (size_t i = 0; i < vertexGroups.size() - 1; ++i) {
        offsets.push_back(offsets.back() + vertexGroups[i].size());
    }

    // 为每组顶点生成线段
    for (size_t i = 0; i < vertexGroups.size(); ++i) {
        const auto& group = vertexGroups[i];
        int groupSize = group.size() / 2;  // 每组有底部和顶部顶点
        int offset = offsets[i];

        for (int j = 0; j < groupSize; ++j) {
            // 垂直线段
            int v1 = offset + j * 2 + 1;
            int v2 = offset + j * 2 + 2;
            file << "l " << v1 << " " << v2 << std::endl;

            // 水平线段(如果不是最后一个点)
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

// 将JSON文件转换为OBJ文件
bool jsonToObj(const std::string& jsonPath, const std::string& outputPath) {
    // 检查JSON文件是否存在
    if (!fileExists(jsonPath)) {
        std::cerr << "JSON文件不存在: " << jsonPath << std::endl;
        return false;
    }

    // 准备数据结构
    std::map<int, std::vector<Point2D>> pointsMap;
    std::vector<Point2D> hullIn, hullOut;
    std::vector<Point2D> splineIn, splineOut;
    std::vector<Point2D> splineIn1, splineOut1;
    std::vector<Point2D> splineIn2, splineOut2;
    std::map<int, double> minZ, minZ1, minZ2;

    // 加载数据
    loadData(jsonPath, pointsMap, hullIn, hullOut, splineIn, splineOut,
        splineIn1, splineOut1, splineIn2, splineOut2, minZ, minZ1, minZ2);

    // 生成.obj文件
    generateObj(splineIn, splineOut, splineIn1, splineOut1, splineIn2, splineOut2,
        minZ, minZ1, minZ2, outputPath);

    return true;
}


Point2D::Point2D(double x_, double y_, int label_) : x(x_), y(y_), label(label_) {}

Point3D::Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
