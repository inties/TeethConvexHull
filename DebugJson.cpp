#include "DebugJson.h"


// 静态成员变量定义
std::string DebugJson::jsonFilePath = "E:\\MylabProjects\\debugPoints.json";

void DebugJson::setJsonPath(const std::string& path) {
    jsonFilePath = path;
}

void DebugJson::writeToJson(const std::map<int, std::vector<Point2DStruct>>& resultMap) {
    writeToJson(resultMap, jsonFilePath);
}

void DebugJson::writeToJson(const std::map<int, std::vector<Point2DStruct>>& resultMap, const std::string& filePath) {
    try {
        nlohmann::json json;
        nlohmann::json pointsmap = nlohmann::json::object();

        for (const auto& pair : resultMap) {
            nlohmann::json points_array = nlohmann::json::array();
            const std::vector<Point2DStruct>& points = pair.second;

            for (const Point2DStruct& point : points) {
                points_array.push_back({
                    {"x", point.x},
                    {"y", point.y}
                    });
            }
            pointsmap[std::to_string(pair.first)] = points_array;
        }

        json["resultpoints"] = pointsmap;
        json["timestamp"] = std::time(nullptr);
        json["point_count"] = resultMap.size();

        // 写入文件
        
        std::ofstream file(filePath);
        if (!file.is_open()) {
            std::cerr << "无法创建JSON文件: " << filePath << std::endl;
            return;
        }

        file << json.dump(4); // 4是缩进空格数，使JSON格式更易读
        file.close();

        std::cout << "成功写入JSON文件: " << filePath << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "写入JSON文件时发生错误: " << e.what() << std::endl;
    }
}
void DebugJson::writeToJson(const std::vector<Point2DStruct>& resultPoints) {
    
    nlohmann::json json;
	nlohmann::json points_array = nlohmann::json::array();
	
	// 写入文件
    std::ifstream inputFile(jsonFilePath);
    if (inputFile.is_open()) {
        try {
            inputFile >> json;
        }
        catch (const nlohmann::json::parse_error& e) {
            std::cerr << "JSON解析错误: " << e.what() << std::endl;
            // 如果解析失败，创建一个新的JSON对象
        }
        inputFile.close();
    }
    for (const Point2DStruct& point : resultPoints) {
        points_array.push_back({
            {"x", point.x},
            {"y", point.y}
            });
    }
    json["debugPoints"] = points_array;
	std::ofstream file(jsonFilePath);
	if (!file.is_open()) {
		std::cerr << "无法创建JSON文件: " << jsonFilePath << std::endl;
		return;
	}
	file << json.dump(4); // 4是缩进空格数，使JSON格式更易读
	file.close();
	std::cout << "成功写入JSON文件: " << jsonFilePath << std::endl;
}


