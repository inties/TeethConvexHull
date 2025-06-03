//使用演示
#include"TeethPipeline.h"
#include<vector>
using namespace std;
int main() {
	vector<float>percents={0.6};
	TeethConvex::TeethPipeline::run("C:\\Users\\邓力源\\Desktop\\model\\yms_upper-pca.obj", "C:\\Users\\邓力源\\Desktop\\model\\yms_upper-pca.json",
		"C:/Users/邓力源/Desktop/model/teethTest_lower_output.obj",percents);
	return 0;
}