#include "stdafx.h"
#include "global_trend.h"

GlobalTrend::GlobalTrend(stock_values& stock)
{
	stocks   = stock.values_time_major;
	stock_each = stock.values;
	ndays    = stock.ndays;
	ntimes   = stock.ntimes;
	nstocks  = stock.nvalues;
	global_value    = new float*[ndays];
	global_value[0] = new float [ndays*ntimes];
	for(int i = 0; i < ndays; ++i)global_value[i] = i*ntimes + global_value[0];
	each_sd         = new float*[ndays];
	each_sd[0]      = new float [ndays*nstocks];
	for(int i = 0; i < ndays; ++i)each_sd[i] = i*nstocks + each_sd[0];
	each_var        = new float*[ndays];
	each_var[0]     = new float [ndays*nstocks];
	for(int i = 0; i < ndays; ++i)each_var[i] = i*nstocks + each_var[0];
}

GlobalTrend::~GlobalTrend() {
	delete [] global_value[0];
	delete [] global_value;
	stocks = nullptr;
	//stock_each = nullptr;
}

float GlobalTrend::CDF_Fit(float *value, const int &length, float &fit_sigma) {
	//Gauss distribution
	//(1/sqrt(2.f*PI*pow(sigma,2)))*exp(-pow((x - mu),2)/(2.f*pow(sigma,2)))
	//CDF for Gauss
	//(1 + erf((x - mu)/sqrt(2.f*pow(sigma,2)))/2.f
	// prepare X and Y
	step_a = 0.02f;
	std::vector<float> X(value, value + length);
	std::sort (X.begin(), X.begin() + length);
	std::vector<float> Y(length);
	for(int i = 0; i < length; ++i) Y[i] = i/(length - 1.f);

	// inital guess of sigma and mu
	sigma = 1.f;
	mu = 0.f;
	float error_prev, error_now(10000.f);
	int niter(0);

	do {
		// prepare (x-mu)/sqrt(2*pow(sigma,2)) and residual
		float sqrt_2_sigma_sq_over_1 (1.f/sqrt(2.f*pow(sigma,2)));
		x_mu_sigma.resize(length);
		for(int i = 0; i < length; ++i)x_mu_sigma[i] = (X[i] - mu)*sqrt_2_sigma_sq_over_1;
		residual.resize(length);
		for(int i = 0; i < length; ++i)residual[i] = (1.f + erff(x_mu_sigma[i]))/2.f - Y[i];
		error_prev = error_now;
		error_now = 0.f;
		for(int i = 0; i < length; ++i)error_now += pow(residual[i],2);

		// calculate slope
		float slope_mu, slope_sigma;
		gd_slope(slope_mu, slope_sigma);

		// update sigma, mu
		mu    = mu    - step_a*slope_mu;
		sigma = sigma - step_a*slope_sigma;
		++niter;
		if(niter > 10000) {
			break;
			std::cout << "convergency error!" <<  std::endl;
		}
	} while (abs(error_prev - error_now)/abs(error_prev + error_now) > 0.000001f);
	fit_sigma = sigma;
	return mu;
}

bool GlobalTrend::gd_slope(float &slope_mu, float &slope_sigma) {
	int   nmax = residual.size();
	float tmp_sigma = 0.f, tmp_mu = 0.f;
	for(int i = 0; i < nmax; ++i) {
		tmp_sigma += 2.f*residual[i]*exp(-pow(x_mu_sigma[i],2))/sqrt(M_PI)*(-x_mu_sigma[i]/sigma);
		tmp_mu    += 2.f*residual[i]*exp(-pow(x_mu_sigma[i],2))/sqrt(M_PI)*(-1.f/sqrt(2.f*pow(sigma,2)));
	}
	slope_sigma = tmp_sigma/(float)nmax;
	slope_mu    = tmp_mu/(float)nmax;
	return true;
}

float GlobalTrend::erff(const float &x) {
	// constants
    float a1 =  0.254829592f;
    float a2 = -0.284496736f;
    float a3 =  1.421413741f;
    float a4 = -1.453152027f;
    float a5 =  1.061405429f;
    float p  =  0.3275911f;
    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    float x_abs = fabs(x);
    // A&S formula 7.1.26
    float t = 1.f/(1.f + p*x_abs);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x_abs*x_abs);
    return sign*y;
}

void GlobalTrend::GetGlobal() {
	float tmp;
	for(int i = 0; i < ndays; ++i) {
		std::cout << "processing day " << i << std::endl;
		for(int j = 0; j < ntimes; ++j) {
			if(j == 0) global_value[i][j] = 0.f;
			else global_value[i][j] = CDF_Fit(stocks[i][j], nstocks, tmp);
		}
	}
}

void GlobalTrend::get_each_sd() {
	float tmp;
	for(int i = 0; i < ndays; ++i) {
		for(int j = 0; j < nstocks; ++j) {
			tmp = 0.f;
			for(int k = 0; k < ntimes; ++k) {
				tmp += stocks[i][k][j];
			}
			each_var[i][j] = tmp/ntimes;
		}
	}
	for(int i = 0; i < ndays; ++i) {
		for(int j = 0; j < nstocks; ++j) {
			tmp = 0.f;
			for(int k = 0; k < ntimes; ++k) {
				tmp += pow((stocks[i][k][j] - each_var[i][j]),2);
			}
			if(each_var[i][j] > 0.f)each_sd[i][j] = sqrt(tmp/ntimes);
			else each_sd[i][j] = -sqrt(tmp/ntimes);
		}
	}
	std::ofstream file;
	file.open("all_sd.txt", std::ios_base::out);
	for(int i = 0; i < ndays; ++i) {
		for(int j = 0; j < nstocks; ++j) {
			file << each_sd[i][j] << ",";
		}
		file << std::endl;
	}

}