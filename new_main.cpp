//ʹ����ʾ
#include"TeethPipeline.h"
#include<vector>
using namespace std;
int main() {
	vector<float>percents={0.6};
	TeethConvex::TeethPipeline::run("C:\\Users\\����Դ\\Desktop\\model\\yms_upper-pca.obj", "C:\\Users\\����Դ\\Desktop\\model\\yms_upper-pca.json",
		"C:/Users/����Դ/Desktop/model/teethTest_lower_output.obj",percents);
	return 0;
}