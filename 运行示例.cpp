//使用演示
#include"TeethPipeline.h"
#include<vector>
using namespace std;
int main() {
	//调用方法1：指定OBJ文件路径，要求obj文件同一目录下有同名的.json文件（牙齿标签），将以默认percents参数处理。输出轴嵴连线obj文件到同一目录下。
	TeethConvex::TeethPipeline::run("E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.obj");
	
	//调用方法2：指定obj文件路径、json文件路径、输出obj文件路径和percents参数，进行处理。
	/*vector<float>percents={1,0.9,0.8,0.7,0.6};
	TeethConvex::TeethPipeline::run("E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.obj", "E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.json",
		"C:/Users/邓力源/Desktop/model/teethTest_lower_output.obj",percents);*/
	
	return 0;
}