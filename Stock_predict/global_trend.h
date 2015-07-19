#include "stdafx.h"
#include "dataload.h"
// value(stock_ind, day_ind, time_ind)
// feature(day_ind, time_ind, feature_ind);

#ifndef _GLOBAL_TREND_H_
#define _GLOBAL_TREND_H_

class GlobalTrend
{
public:
	// global_value(day_ind, time_ind)
	float **global_value;
	// global_value(day_ind, stock_ind)
	float **each_sd;
	float **each_var;
	GlobalTrend(stock_values& stock);
	~GlobalTrend();
	void GetGlobal();
	void get_each_sd();
	float CDF_Fit(float *value, const int &length, float &fit_sigma);
	void arma_fit();
private:
	int ndays, ntimes, nstocks;
	float sigma, mu;
	float step_a;
	std::vector<float> x_mu_sigma;     //(x-mu)/sqrt(2*pow(sigma,2))
	std::vector<float> residual;       //residual
	float ***stocks, ***stock_each;
	bool gd_slope(float &slope_mu, float &slope_sigma);
	float erff(const float &x);
};


#endif