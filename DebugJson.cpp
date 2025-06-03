#include "DebugJson.h"


// ��̬��Ա��������
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

        // д���ļ�
        
        std::ofstream file(filePath);
        if (!file.is_open()) {
            std::cerr << "�޷�����JSON�ļ�: " << filePath << std::endl;
            return;
        }

        file << json.dump(4); // 4�������ո�����ʹJSON��ʽ���׶�
        file.close();

        std::cout << "�ɹ�д��JSON�ļ�: " << filePath << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "д��JSON�ļ�ʱ��������: " << e.what() << std::endl;
    }
}
void DebugJson::writeToJson(const std::vector<Point2DStruct>& resultPoints) {
    
    nlohmann::json json;
	nlohmann::json points_array = nlohmann::json::array();
	
	// д���ļ�
    std::ifstream inputFile(jsonFilePath);
    if (inputFile.is_open()) {
        try {
            inputFile >> json;
        }
        catch (const nlohmann::json::parse_error& e) {
            std::cerr << "JSON��������: " << e.what() << std::endl;
            // �������ʧ�ܣ�����һ���µ�JSON����
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
		std::cerr << "�޷�����JSON�ļ�: " << jsonFilePath << std::endl;
		return;
	}
	file << json.dump(4); // 4�������ո�����ʹJSON��ʽ���׶�
	file.close();
	std::cout << "�ɹ�д��JSON�ļ�: " << jsonFilePath << std::endl;
}


