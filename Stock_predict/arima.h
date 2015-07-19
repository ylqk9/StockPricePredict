#include "stdafx.h"

class arima
{
private:
	unsigned ARIMA_p, ARIMA_d, ARIMA_q;
	bool     ARIMA_predicted, SD_Var_calculated;
	unsigned diff_order, ar_order, ma_order;
	float    sd, expectation;
	unsigned autocov_order;
	std::vector<float> autocorrelation, autocovariance;
	std::vector<float> diff(const std::vector<float> array_in, const unsigned &d);
	std::vector<float> anti_diff(const std::vector<float> array_in, const unsigned &d);
	void               Var_And_SD(const std::vector<float> array_in, float &expectation, float &sd);
	std::vector<float> normalize(const std::vector<float> array_in);
	std::vector<float> anti_normalize(const std::vector<float> array_in);
	std::vector<float> get_autocovariance(const std::vector<float> array_in);
	std::vector<float> get_autocorrelation(const std::vector<float> array_in);
	std::vector<float> AR_coeffs, MA_coeffs;
public:
	                   arima();
					   ~arima();
	std::vector<float> ARIMAmodel(const unsigned &p, const unsigned &d, const unsigned &q, std::vector<float> array_in);
	std::vector<float> AR_ols(const std::vector <float> array_in, const int &p);
	std::vector<float> AR_YuleWalkerEq(const std::vector<float> array_in, const int &p);
	std::vector<float> MA_ML_estimation(const std::vector<float> array_in, const int &q);
	// nstep forecast
	std::vector<float> AR_forecast(const std::vector<float> &array_in, const int &nsteps);
	std::vector<float> MA_forecast(const std::vector<float> &array_in, const int &nsteps);
	// single step forecast
	float              MA_forecast(const std::vector<float> &array_in);
	std::vector<float> ARIMAforecast(const std::vector<float> &array_in, const int &nsteps);
	unsigned           ARIMAget_p() {return ARIMA_p;};
	unsigned           ARIMAget_d() {return ARIMA_d;};
	unsigned           ARIMAget_q() {return ARIMA_q;};
	std::vector<float> ARIMAget_coeffs(const std::string &type);
};