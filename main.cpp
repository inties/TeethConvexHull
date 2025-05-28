//#include <iostream>
//#include <string>
//#include <fstream>
//#include "convex.h"
//#include "file_utils.h"
//
//int main(int argc, char* argv[]) {
//    if (argc < 2) {
//        std::cerr << "�÷�: " << argv[0] << " <input_obj_file> [percent1] [percent2] [percent3]" << std::endl;
//        std::cerr << "ע��: ������ṩpercent��������ʹ��Ĭ��ֵ: 1.0, 0.7, 0.5" << std::endl;
//        return 1;
//    }
//
//    std::string input_obj_file = argv[1];
//
//    // ��������ļ��Ƿ����
//    if (!fileExists(input_obj_file)) {
//        std::cerr << "���������ļ�������: " << input_obj_file << std::endl;
//        return 1;
//    }
//
//    // ����percent����
//    if (argc == 2) {
//        // ʹ��Ĭ��ֵ
//        std::cout << "ʹ��Ĭ��percentֵ: 1.0, 0.7, 0.5" << std::endl;
//        bool success = processTeethData(input_obj_file);
//
//        if (success) {
//            std::cout << "����ɹ���ɣ�" << std::endl;
//            return 0;
//        }
//        else {
//            std::cerr << "����ʧ�ܣ�" << std::endl;
//            return 1;
//        }
//    }
//    else if (argc == 5) {
//        // ʹ������percentֵ
//        try {
//            float p1 = std::stof(argv[2]);
//            float p2 = std::stof(argv[3]);
//            float p3 = std::stof(argv[4]);
//
//            // ��֤percentֵ
//            if (p1 <= 0.0f || p1 > 1.0f || p2 <= 0.0f || p2 > 1.0f || p3 <= 0.0f || p3 > 1.0f) {
//                std::cerr << "����percentֵ������0��1֮��" << std::endl;
//                return 1;
//            }
//
//            // ��֤percentֵ�Ƿ�ݼ�
//            if (p1 < p2 || p2 < p3) {
//                std::cerr << "����percentֵ����ݼ� (p1 >= p2 >= p3)" << std::endl;
//                return 1;
//            }
//
//            std::cout << "ʹ���Զ���percentֵ: " << p1 << ", " << p2 << ", " << p3 << std::endl;
//            bool success = processTeethData(input_obj_file, p1, p2, p3);
//
//            if (success) {
//                std::cout << "����ɹ���ɣ�" << std::endl;
//                return 0;
//            }
//            else {
//                std::cerr << "����ʧ�ܣ�" << std::endl;
//                return 1;
//            }
//        }
//        catch (const std::exception& e) {
//            std::cerr << "������Ч��percent����: " << e.what() << std::endl;
//            return 1;
//        }
//    }
//    else {
//        std::cerr << "���󣺱����ṩ����percent������ʹ��Ĭ��ֵ" << std::endl;
//        std::cerr << "�÷�: " << argv[0] << " <input_obj_file> [percent1] [percent2] [percent3]" << std::endl;
//        return 1;
//    }
//}