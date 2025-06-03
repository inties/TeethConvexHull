#pragma once
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

struct Point2DStruct {
	float x;
	float y;
	Point2DStruct(float x, float y) {
		this->x = x;
		this->y = y;
	}
};
using namespace std;
class DebugJson
{
private:
	static std::string jsonFilePath;

public:
	// ����JSON�ļ�·��
	static void setJsonPath(const std::string& path);

	// д��㼯�����ݵ�JSON�ļ�
	static void writeToJson(const std::map<int, std::vector<Point2DStruct>>& resultMap);
	static void writeToJson(const std::vector<Point2DStruct>& resultPoints);

	// д��㼯�����ݵ�ָ��·����JSON�ļ�
	static void writeToJson(const std::map<int, std::vector<Point2DStruct>>& resultMap, const std::string& filePath);

};
