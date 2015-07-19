// Stock_predict.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "dataload.h"
#include "global_trend.h"
#include "regression.h"
#include "arima.h"

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned days(510), p(2);
	std::ofstream file;
	file.open("result.txt", std::ios_base::out);
	StockData training(198, 244, 55, days);
	training.BatchLoad();
	std::cout << "1.all file loaded!" << std::endl;
	GlobalTrend training_global(training.value);
	//training_global.GetGlobal();
	//std::cout << "2.global fit done!" << std::endl;
	//for(int i = 0; i < days; ++i) {
	//	std::cout << "AR day " << i << std::endl;
	//	arima mymodel;
	//	std::vector<float> array(training_global.global_value[i], training_global.global_value[i] + 55);
	//	mymodel.AR_ols(array, p);
	//	std::vector<float> outcome(mymodel.AR_forecast(array, 25));
	//	file << outcome[24] << std::endl;
	//}
	training_global.get_each_sd();
	file.close();
	return 0;
}

