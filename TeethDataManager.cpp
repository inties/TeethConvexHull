#include "TeethDataManager.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <iostream>
#include <CGAL/convex_hull_2.h>
using json = nlohmann::json;
using std::vector;
namespace TeethConvex {
    TeethDataManager::TeethDataManager(const std::string& obj_path, const std::string& json_path)
        : objPath(obj_path), jsonPath(json_path) {
    }
    void TeethDataManager::load() {
        initializeOriginPoints();
        zfilteredRanges = zRanges; 
    }

    void TeethDataManager::initializeOriginPoints() {
        // 读取OBJ文件
        std::ifstream obj_file(objPath);
        if (!obj_file.is_open()) {
            // 错误处理：文件无法打开
            std::cerr << "无法打开 OBJ 文件: " << objPath << std::endl;
            return;
        }
        std::string line;
        std::vector<Point_3> allPoints;
        allPoints.reserve(100000);
        while (std::getline(obj_file, line)) {
            // 避免 substr(0, 2)
            if (line.length() > 1 && line[0] == 'v' && line[1] == ' ') {
                // 找到第一个空格后的位置
                size_t pos = 2; // 从 'v ' 后开始查找
                size_t end_pos;

                double x, y, z;

                // 解析 X
                end_pos = line.find(' ', pos);
                if (end_pos == std::string::npos) continue; // 格式错误
                x = std::strtod(line.c_str() + pos, nullptr); // 使用 strtod 解析

                // 解析 Y
                pos = end_pos + 1;
                end_pos = line.find(' ', pos);
                if (end_pos == std::string::npos) continue; // 格式错误
                y = std::strtod(line.c_str() + pos, nullptr);

                // 解析 Z
                pos = end_pos + 1;
                z = std::strtod(line.c_str() + pos, nullptr); // Z 后面可能没有空格，直接到行尾

                allPoints.emplace_back(x, y, z);
            }
        }
        obj_file.close(); 
        // 读取JSON标签文件
        std::ifstream json_file(jsonPath);
        json jsonData;
        json_file >> jsonData;
        std::vector<int> labels = jsonData["labels"];

        // 按标签分组并计算z值范围
        std::map<int, std::vector<Point_3>> tempPointsMap;
        for (size_t i = 0; i < allPoints.size(); ++i) {
            int label = labels[i];
            if (label == 0) continue;
            double z = allPoints[i].z();
            if (zRanges.find(label) == zRanges.end()) {
                zRanges[label] = { z, z };
            }
            else {
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

        // 初始化zfilteredRanges
        zfilteredRanges = zRanges;
        allPoints.clear();
        labels.clear();
    }



    // 新实现：去掉const，更新zfilteredRanges
    std::map<int, std::vector<Point_2>> TeethDataManager::filterPointsByPercent(float percent) {
        std::map<int, std::vector<Point_2>> pointsMap;  
        zfilteredRanges.clear();
        for (const auto& pair : originPointsMap) {
            int label = pair.first;
            const auto& range = zRanges.at(label);
            double min_z = range.first;
            double max_z = range.second;
            double threshold = max_z - percent * (max_z - min_z);
            double cur_min_z = std::numeric_limits<double>::max();
            double cur_max_z = -std::numeric_limits<double>::max();
            for (const auto& point : pair.second) {
                if (point.z() >= threshold && point.z() <= max_z) {
                    pointsMap[label].emplace_back(point.x(), point.y());
                    if (point.z() < cur_min_z) cur_min_z = point.z();
                    if (point.z() > cur_max_z) cur_max_z = point.z();
                }
            }
            
            if (!pointsMap[label].empty()) {
                zfilteredRanges[label] = { cur_min_z, cur_max_z };
            }
		}
        return pointsMap;
    }
    std::map<int, std::vector<Point_2>> TeethDataManager::calculateConvexByTeeth(std::map<int,vector<Point_2>>&pointsmap) {
		std::map<int, std::vector<Point_2>> result;
        for (auto pair : pointsmap) {
            vector<Point_2>temp;
			vector<Point_2>& points = pair.second;
			CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(temp));
			result[pair.first] = std::move(temp);
        }
        return result;
    }
    const std::map<int, std::vector<Point_3>>& TeethDataManager::getOriginPointsMap() const {
        return originPointsMap;
    }

    const std::map<int, std::pair<double, double>>& TeethDataManager::getZRanges() const {
        return zfilteredRanges;
    }
}
