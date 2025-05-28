#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <string>
#include <fstream>

// 检查文件是否存在
inline bool fileExists(const std::string& filePath) {
    std::ifstream file(filePath);
    return file.good();
}

// 获取文件所在目录
inline std::string getDirectory(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    return (pos != std::string::npos) ? filePath.substr(0, pos) : "";
}

// 获取文件名（不含扩展名）
inline std::string getBaseFileName(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    std::string fileName = (pos != std::string::npos) ? filePath.substr(pos + 1) : filePath;
    pos = fileName.find_last_of(".");
    return (pos != std::string::npos) ? fileName.substr(0, pos) : fileName;
}
// 获取文件名部分
inline std::string getFileName(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    if (pos != std::string::npos) {
        return filePath.substr(pos + 1);
    }
    return filePath;
}
#endif // FILE_UTILS_H 