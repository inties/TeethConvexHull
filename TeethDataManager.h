#pragma once
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
namespace TeethConvex {
    class TeethDataManager {
    public:
        TeethDataManager(const std::string& obj_path, const std::string& json_path);
        void load();
        std::map<int, std::vector<Point_2>> filterPointsByPercent(float percent);
        std::map<int, std::vector<Point_2>> calculateConvexByTeeth(std::map<int,std::vector<Point_2>>&pointsmap);
        const std::map<int, std::vector<Point_3>>& getOriginPointsMap() const;
        const std::map<int, std::pair<double, double>>& getZRanges() const;
    private:
		
        std::string objPath, jsonPath;
        std::map<int, std::vector<Point_3>> originPointsMap;
        std::map<int, std::pair<double, double>> zRanges;
        std::map<int, std::pair<double, double>> zfilteredRanges;
        void initializeOriginPoints();
    };
}