//ʹ����ʾ
#include"TeethPipeline.h"
#include<vector>
using namespace std;
int main() {
	vector<float>percents={1,0.9,0.8,0.7,0.6,0.5};
	TeethConvex::TeethPipeline::run("C:/Users/����Դ/Desktop/model/teethTest_lower.obj", "C:/Users/����Դ/Desktop/model/teethTest_lower.json",
		"C:/Users/����Դ/Desktop/model/teethTest_lower_output.obj",percents);
	return 0;
}