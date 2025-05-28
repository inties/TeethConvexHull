#ifndef JSON2OBJ_H
#define JSON2OBJ_H

#include <string>
#include <vector>
#include <map>
#include"file_utils.h"

// �����Ľṹ
struct Point2D {
    double x, y;
    int label;

    Point2D(double x_ = 0, double y_ = 0, int label_ = -1);
};

// �����Ľṹ����z���꣩
struct Point3D {
    double x, y, z;

    Point3D(double x_ = 0, double y_ = 0, double z_ = 0);
};



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
    std::map<int, double>& minZ2);

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
    double zMin = -10.0,
    double zMax = 10.0);

// ��JSON�ļ�ת��ΪOBJ�ļ�
bool jsonToObj(const std::string& jsonPath, const std::string& outputPath);

#endif // JSON2OBJ_H#pragma once
