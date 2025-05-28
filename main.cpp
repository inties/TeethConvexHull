//#include <iostream>
//#include <string>
//#include <fstream>
//#include "convex.h"
//#include "file_utils.h"
//
//int main(int argc, char* argv[]) {
//    if (argc < 2) {
//        std::cerr << "用法: " << argv[0] << " <input_obj_file> [percent1] [percent2] [percent3]" << std::endl;
//        std::cerr << "注意: 如果不提供percent参数，将使用默认值: 1.0, 0.7, 0.5" << std::endl;
//        return 1;
//    }
//
//    std::string input_obj_file = argv[1];
//
//    // 检查输入文件是否存在
//    if (!fileExists(input_obj_file)) {
//        std::cerr << "错误：输入文件不存在: " << input_obj_file << std::endl;
//        return 1;
//    }
//
//    // 处理percent参数
//    if (argc == 2) {
//        // 使用默认值
//        std::cout << "使用默认percent值: 1.0, 0.7, 0.5" << std::endl;
//        bool success = processTeethData(input_obj_file);
//
//        if (success) {
//            std::cout << "处理成功完成！" << std::endl;
//            return 0;
//        }
//        else {
//            std::cerr << "处理失败！" << std::endl;
//            return 1;
//        }
//    }
//    else if (argc == 5) {
//        // 使用三个percent值
//        try {
//            float p1 = std::stof(argv[2]);
//            float p2 = std::stof(argv[3]);
//            float p3 = std::stof(argv[4]);
//
//            // 验证percent值
//            if (p1 <= 0.0f || p1 > 1.0f || p2 <= 0.0f || p2 > 1.0f || p3 <= 0.0f || p3 > 1.0f) {
//                std::cerr << "错误：percent值必须在0到1之间" << std::endl;
//                return 1;
//            }
//
//            // 验证percent值是否递减
//            if (p1 < p2 || p2 < p3) {
//                std::cerr << "错误：percent值必须递减 (p1 >= p2 >= p3)" << std::endl;
//                return 1;
//            }
//
//            std::cout << "使用自定义percent值: " << p1 << ", " << p2 << ", " << p3 << std::endl;
//            bool success = processTeethData(input_obj_file, p1, p2, p3);
//
//            if (success) {
//                std::cout << "处理成功完成！" << std::endl;
//                return 0;
//            }
//            else {
//                std::cerr << "处理失败！" << std::endl;
//                return 1;
//            }
//        }
//        catch (const std::exception& e) {
//            std::cerr << "错误：无效的percent参数: " << e.what() << std::endl;
//            return 1;
//        }
//    }
//    else {
//        std::cerr << "错误：必须提供三个percent参数或使用默认值" << std::endl;
//        std::cerr << "用法: " << argv[0] << " <input_obj_file> [percent1] [percent2] [percent3]" << std::endl;
//        return 1;
//    }
//}