#include "TeethPipeline.h"
#include "TeethDataManager.h"
#include "TeethCurveFitter.h"
#include "TeethResultExporter.h"
#include <iostream>
#include <future>
#include <vector>
#include <algorithm>
#include<mutex>

namespace TeethConvex {
    // 定义结果结构体，用于多线程安全处理
    struct IterationResult {
        std::map<int, std::vector<Point_2>> pointsMap;
        std::map<int, std::pair<double, double>> zRanges;
        std::vector<Point_2> innerHull;
        std::vector<Point_2> outerHull;
        std::map<int, std::vector<Point_2>> innerHullMap;
        std::map<int, std::vector<Point_2>> outerHullMap;
        std::vector<Point_2> innerBSpline;
        std::vector<Point_2> outerBSpline;
        size_t index;  // 保存迭代索引，确保结果顺序
    };
    void TeethPipeline::run(
        const std::string& objPath,
        const std::string& jsonPath,
        const std::string& outputObjPath,
        const std::vector<float>& percents
    ) {
        TeethDataManager dataManager(objPath, jsonPath);
        dataManager.load();

        TeethResultExporter exporter;
        std::mutex mutex_;
        // 用于存储每个线程的 future
        std::vector<std::future<IterationResult>> futures;

        for (size_t i = 0; i < percents.size(); ++i) {
            float percent = percents[i];

            // 启动异步任务，返回计算结果而不是直接修改共享数据
            futures.push_back(std::async(std::launch::async, [&mutex_,&dataManager, percent, i]() -> IterationResult {
                {
                    std::lock_guard<std::mutex>lock(mutex_);
                    std::cout << "Pipeline: processing percent = " << percent << " (iteration " << i << ")" << std::endl;
                }
                // 创建结果结构体
                IterationResult result;
                result.index = i;

                {
					std::lock_guard<std::mutex> lock(mutex_);
                    std::map<int,std::vector<Point_2>>pointsmap= dataManager.filterPointsByPercent(percent);
					result.pointsMap = dataManager.calculateConvexByTeeth(pointsmap);
                    result.zRanges = dataManager.getZRanges();
                }
                TeethCurveFitter fitter(result.pointsMap);
                fitter.computeConvexHulls();

                result.innerHull = fitter.getInnerHull();
                result.outerHull = fitter.getOuterHull();
                result.innerHullMap = fitter.getInnerHullMap();
                result.outerHullMap = fitter.getOuterHullMap();
                result.innerBSpline = fitter.getInnerBSpline();
                result.outerBSpline = fitter.getOuterBSpline();

                return result;
                }));
        }

        std::vector<IterationResult> results;
        results.reserve(futures.size());

        // 等待所有线程完成并收集结果
        for (auto& f : futures) {
            try {
                results.push_back(f.get());
            }
            catch (const std::exception& e) {
                std::cerr << "处理异常: " << e.what() << std::endl;
            }
        }


        std::sort(results.begin(), results.end(),
            [](const IterationResult& a, const IterationResult& b) {
                return a.index < b.index;
            });

        //for (const auto& result : results) {
        //    exporter.addIterationData(
        //        result.pointsMap,
        //        result.zRanges,
        //        result.innerHull,
        //        result.outerHull,
        //        result.innerHullMap,
        //        result.outerHullMap,
        //        result.innerBSpline,
        //        result.outerBSpline
        //    );
        //}
        //exporter.exportJSON("E:\\MylabProjects\\debug.json");

        //exporter.exportOBJ(outputObjPath);
        //exporter.exportOBJ_2("E:/MylabProjects/debug");
        std::cout << "Pipeline: Finished processing. Results saved to " << outputObjPath << std::endl;
    }

}