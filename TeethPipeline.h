#pragma once
#include <string>
#include <vector>
namespace TeethConvex {
    class TeethPipeline {
    public:
        /**
         * 运行整个处理管道
         * @param objPath 输入OBJ文件路径
         * @param jsonPath 输入JSON标签文件路径
         * @param outputJsonPath 输出结果JSON文件路径
         * @param outputObjPath 输出结果OBJ文件路径
         * @param percents 处理百分比数组
         */
        static void run(
            const std::string& objPath,
            const std::string& jsonPath,
            const std::string& outputObjPath,
            const std::vector<float>& percents
        );
    };
}