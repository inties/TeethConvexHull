#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <string>
#include <fstream>

// ����ļ��Ƿ����
inline bool fileExists(const std::string& filePath) {
    std::ifstream file(filePath);
    return file.good();
}

// ��ȡ�ļ�����Ŀ¼
inline std::string getDirectory(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    return (pos != std::string::npos) ? filePath.substr(0, pos) : "";
}

// ��ȡ�ļ�����������չ����
inline std::string getBaseFileName(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    std::string fileName = (pos != std::string::npos) ? filePath.substr(pos + 1) : filePath;
    pos = fileName.find_last_of(".");
    return (pos != std::string::npos) ? fileName.substr(0, pos) : fileName;
}
// ��ȡ�ļ�������
inline std::string getFileName(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    if (pos != std::string::npos) {
        return filePath.substr(pos + 1);
    }
    return filePath;
}
#endif // FILE_UTILS_H 