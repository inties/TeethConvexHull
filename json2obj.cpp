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









// ��JSON�ļ���������
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

    // ��JSON�ļ�
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "�޷����ļ�: " << filePath << std::endl;
        return;
    }

    // ����JSON
    json data;
    file >> data;

    // ��ȡpointsMap
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

    // ��ȡhull_in
    if (data.contains("hull_in")) {
        for (const auto& point : data["hull_in"]) {
            hullIn.emplace_back(point["x"], point["y"]);
        }
    }

    // ��ȡhull_out
    if (data.contains("hull_out")) {
        for (const auto& point : data["hull_out"]) {
            hullOut.emplace_back(point["x"], point["y"]);
        }
    }

    // ��ȡspline_in
    if (data.contains("spline_in")) {
        for (const auto& point : data["spline_in"]) {
            splineIn.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡspline_out
    if (data.contains("spline_out")) {
        for (const auto& point : data["spline_out"]) {
            splineOut.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡspline_in_1
    if (data.contains("spline_in_1")) {
        for (const auto& point : data["spline_in_1"]) {
            splineIn1.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡspline_out_1
    if (data.contains("spline_out_1")) {
        for (const auto& point : data["spline_out_1"]) {
            splineOut1.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡspline_in_2
    if (data.contains("spline_in_2")) {
        for (const auto& point : data["spline_in_2"]) {
            splineIn2.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡspline_out_2
    if (data.contains("spline_out_2")) {
        for (const auto& point : data["spline_out_2"]) {
            splineOut2.emplace_back(point["x"], point["y"], point["label"]);
        }
    }

    // ��ȡmin_z
    if (data.contains("min_z")) {
        for (const auto& item : data["min_z"].items()) {
            minZ[std::stoi(item.key())] = item.value();
        }
    }

    // ��ȡmin_z_1
    if (data.contains("min_z_1")) {
        for (const auto& item : data["min_z_1"].items()) {
            minZ1[std::stoi(item.key())] = item.value();
        }
    }

    // ��ȡmin_z_2
    if (data.contains("min_z_2")) {
        for (const auto& item : data["min_z_2"].items()) {
            minZ2[std::stoi(item.key())] = item.value();
        }
    }

    std::cout << "�ɹ���������" << std::endl;
}

// ����.obj�ļ�
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
        std::cerr << "�޷������ļ�: " << outputPath << std::endl;
        return;
    }

    int vertexCount = 0;
    std::vector<std::vector<Point3D>> vertexGroups;

    // ����spline_in (ʹ��minZ��minZ1)
    std::vector<Point3D> vertices;
    for (const auto& point : splineIn) {
        int label = point.label;
        double z1 = (minZ.find(label) != minZ.end()) ? minZ.at(label) : zMin;
        double z2 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // ����spline_out (ʹ��minZ��minZ1)
    vertices.clear();
    for (const auto& point : splineOut) {
        int label = point.label;
        double z1 = (minZ.find(label) != minZ.end()) ? minZ.at(label) : zMin;
        double z2 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // ����spline_in_1 (ʹ��minZ1��minZ2)
    vertices.clear();
    for (const auto& point : splineIn1) {
        int label = point.label;
        double z1 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMin;
        double z2 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // ����spline_out_1 (ʹ��minZ1��minZ2)
    vertices.clear();
    for (const auto& point : splineOut1) {
        int label = point.label;
        double z1 = (minZ1.find(label) != minZ1.end()) ? minZ1.at(label) : zMin;
        double z2 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMax;
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // ����spline_in_2 (ʹ��minZ2��zMax)
    vertices.clear();
    for (const auto& point : splineIn2) {
        int label = point.label;
        double z1 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMin;
        double z2 = z1 + 5.0;  // ����(û�и��߲�)
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // ����spline_out_2 (ʹ��minZ2��zMax)
    vertices.clear();
    for (const auto& point : splineOut2) {
        int label = point.label;
        double z1 = (minZ2.find(label) != minZ2.end()) ? minZ2.at(label) : zMin;
        double z2 = z1 + 5.0;  // ����(û�и��߲�)
        vertices.emplace_back(point.x, point.y, z1);  // �ײ�
        vertices.emplace_back(point.x, point.y, z2);  // ����
    }
    vertexGroups.push_back(vertices);

    // д�����ж���
    for (const auto& group : vertexGroups) {
        for (const auto& v : group) {
            file << "v " << v.x << " " << v.y << " " << v.z << std::endl;
            vertexCount++;
        }
    }

    // д���߶� (l)
    // ����ÿ�鶥���ƫ����
    std::vector<int> offsets = { 0 };
    for (size_t i = 0; i < vertexGroups.size() - 1; ++i) {
        offsets.push_back(offsets.back() + vertexGroups[i].size());
    }

    // Ϊÿ�鶥�������߶�
    for (size_t i = 0; i < vertexGroups.size(); ++i) {
        const auto& group = vertexGroups[i];
        int groupSize = group.size() / 2;  // ÿ���еײ��Ͷ�������
        int offset = offsets[i];

        for (int j = 0; j < groupSize; ++j) {
            // ��ֱ�߶�
            int v1 = offset + j * 2 + 1;
            int v2 = offset + j * 2 + 2;
            file << "l " << v1 << " " << v2 << std::endl;

            // ˮƽ�߶�(����������һ����)
            if (j < groupSize - 1) {
                v2 = offset + j * 2 + 2;
                int v4 = offset + (j + 1) * 2 + 2;
                file << "l " << v2 << " " << v4 << std::endl;
            }
        }
    }

    std::cout << "������.obj�ļ�: " << outputPath << std::endl;
    std::cout << "���� " << vertexCount << " ������" << std::endl;
}

// ��JSON�ļ�ת��ΪOBJ�ļ�
bool jsonToObj(const std::string& jsonPath, const std::string& outputPath) {
    // ���JSON�ļ��Ƿ����
    if (!fileExists(jsonPath)) {
        std::cerr << "JSON�ļ�������: " << jsonPath << std::endl;
        return false;
    }

    // ׼�����ݽṹ
    std::map<int, std::vector<Point2D>> pointsMap;
    std::vector<Point2D> hullIn, hullOut;
    std::vector<Point2D> splineIn, splineOut;
    std::vector<Point2D> splineIn1, splineOut1;
    std::vector<Point2D> splineIn2, splineOut2;
    std::map<int, double> minZ, minZ1, minZ2;

    // ��������
    loadData(jsonPath, pointsMap, hullIn, hullOut, splineIn, splineOut,
        splineIn1, splineOut1, splineIn2, splineOut2, minZ, minZ1, minZ2);

    // ����.obj�ļ�
    generateObj(splineIn, splineOut, splineIn1, splineOut1, splineIn2, splineOut2,
        minZ, minZ1, minZ2, outputPath);

    return true;
}


Point2D::Point2D(double x_, double y_, int label_) : x(x_), y(y_), label(label_) {}

Point3D::Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
