//ʹ����ʾ
#include"TeethPipeline.h"
#include<vector>
using namespace std;
int main() {
	//���÷���1��ָ��OBJ�ļ�·����Ҫ��obj�ļ�ͬһĿ¼����ͬ����.json�ļ������ݱ�ǩ��������Ĭ��percents�������������������obj�ļ���ͬһĿ¼�¡�
	TeethConvex::TeethPipeline::run("E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.obj");
	
	//���÷���2��ָ��obj�ļ�·����json�ļ�·�������obj�ļ�·����percents���������д���
	/*vector<float>percents={1,0.9,0.8,0.7,0.6};
	TeethConvex::TeethPipeline::run("E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.obj", "E:\\LabFiles\\check\\Teeth3DS\\data_part_5\\lower\\0166M778\\0166M778_lower.json",
		"C:/Users/����Դ/Desktop/model/teethTest_lower_output.obj",percents);*/
	
	return 0;
}