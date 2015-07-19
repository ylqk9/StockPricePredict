#include "stdafx.h"
#include "arima.h"
#include "regression.h"
using namespace std;

arima::arima() : autocov_order(15) {
	SD_Var_calculated = false;
}

arima::~arima() {
}

// TODO: copy constructor need to modify such SD_Var_calculated set to false

vector<float> arima::diff(const vector<float> array_in, const unsigned &d) {
	unsigned narray = array_in.size();
	vector<float> array_tmp;
	vector<float> array_out(array_in);
	for(int order = 1; order <= d; ++order) {
		array_tmp = array_out;
		for(int i = 0; i < (narray - order); ++i) {
			array_out[i] = array_tmp[i + 1]- array_tmp[i];
		}
	}
	// resize(x) will keep the first x elements
	array_out.resize(narray - d);
	return array_out;
}

vector<float> arima::anti_diff(const vector<float> array_in, const unsigned &d){
	return;
}

void arima::Var_And_SD(const std::vector<float> array_in, float &expectation, float &sd){
	unsigned narray = array_in.size();
	sd = 0.f;
	expectation = 0.f;
	float sd_sum(0.f), expectation_sum(0.f);
	for(int i = 0; i < narray; ++i)expectation_sum += array_in[i];
	expectation = expectation_sum/narray;
	for(int i = 0; i < narray; ++i)sd_sum += pow((array_in[i] - expectation_sum),2);
	sd = sqrt(sd_sum/narray);
	SD_Var_calculated = true;
	return;
}

vector<float> arima::normalize(const vector<float> array_in) {
	unsigned narray = array_in.size();
	vector<float> array_out(narray);
	Var_And_SD(array_in, expectation, sd);
	for(int i = 0; i < narray; ++i)array_out[i] = (array_in[i] - expectation)/sd;
	return array_out;
}

vector<float> arima::anti_normalize(const vector<float> array_in) {
}

vector<float> arima::get_autocovariance(const vector<float> array_in) {
	unsigned narray = array_in.size();
	if(!SD_Var_calculated)Var_And_SD(array_in, expectation, sd);
	float cov = pow(sd, 2);
	autocovariance.resize(autocov_order);
	for(int i = 1; i < autocov_order; ++i) {
		float sum_tmp(0.f);
		for(int j = 0; j < (narray - i); ++j) {
			sum_tmp = (array_in[j] - expectation)*(array_in[j+i] - expectation);
		}
		autocorrelation[i] = sum_tmp/(narray - i);
	}
	return autocovariance;
}

vector<float> arima::get_autocorrelation(vector<float> array_in) {
	Var_And_SD(array_in, expectation, sd);
	autocorrelation = get_autocovariance(array_in);
	for(int i = 0; i < autocov_order; ++i) {
		autocorrelation[i] = autocorrelation[i]/autocorrelation[0];
	}
	return autocorrelation;
}

vector<float> arima::AR_ols(const vector <float> array_in, const int &p){
	if(p == 0) {
		vector<float> a(1);
		a[0] = 0.f;
		return a;
	}
	ARIMA_p = p;
	unsigned narray = array_in.size();
	float **x = new float*[(narray - p)];
	x[0] = new float[(p + 1)*(narray - p)];
	for(int i = 1; i < (narray - p); ++i) x[i] = (p + 1)*i + x[0];
	float *y = new float[narray - p];
	for(int i = 0; i < narray - p; ++i) {
		for(int j = 0; j < (p + 1); ++j) {
			if(j == 0)x[i][j] = 1.f;
			else x[i][j] = array_in[i + j - 1];
		}
	}
	for(int i = 0; i < narray - p; ++i) y[i] = array_in[i + p];
	float *result = regression('r', x[0], y, (p + 1), (narray - p));
	AR_coeffs.resize(p + 1);
	for(int i = 0; i < p + 1; ++i)AR_coeffs[i] = *(result + i);
	free (result);
	return AR_coeffs;
}

vector<float> arima::AR_YuleWalkerEq(const vector<float> array_in, const int &p) {

}

vector<float> arima::MA_ML_estimation(const vector<float> array_in, const int &q) {
	if(q == 0) {
		vector<float> a(1);
		a[0] = 0.f;
		return a;
	}
	unsigned narray = array_in.size();
	MA_coeffs.resize(q);
	//initial guess. No reason.
	for(int i = 0; i < q; ++i)MA_coeffs[i] = 0.5f;
	float sd_MLE(1.f);
	float error_prev, error_now(10000.f);
	float step_a = 0.1f;

	vector<float> residual(narray, 0.f);
	vector<float> error(narray, 0.f);
	do {
		error_prev = error_now;
		error_now = 0.f;
		for(int i = 0; i < narray; ++i)error_now += pow(residual[i],2);
		error[0] = array_in[0];
		for(int i = 1; i < narray; ++i) {
			int step_q(q), step_i(i);
			do {
				if((step_q < 1) || (step_i < 1)) break;
				error[i] += MA_coeffs[step_q - 1]*error[step_i - 1];
				--step_q;
				--step_i;
			} while(1);
			error[i] = array_in[i] - error[i];
		}

		// calculate slope
	//	vector<float> slope(q, 0);
	//	for(int i = 0; i < q; ++i) {
	//		for(int j = 0; j < 
	//		slope[i] = 0.5f*error_now*error[;

	//	// update sigma, mu
	//	mu    = mu    - step_a*slope_mu;
	//	sigma = sigma - step_a*slope_sigma;
	//} 
	//fit_sigma = sigma;
	//return mu;
	} while (abs(error_prev - error_now)/abs(error_prev + error_now) > 0.000001f);
}


vector<float> arima::ARIMAmodel(const unsigned &p, const unsigned &d, const unsigned &q, vector<float> array_in) {
	ARIMA_p = p; ARIMA_d = d; ARIMA_q = q;
	vector<float> array_diff;
	if(d > 0) array_diff = diff(array_in, d);
	else array_diff = array_in;
	AR_coeffs = AR_ols(array_diff, p);
}

vector<float> arima::AR_forecast(const vector<float> &array_in, const int &nsteps) {
	unsigned p = ARIMAget_p();
	unsigned narray = array_in.size();
	vector<float> array_out((nsteps + p), 0.f);
	for(int i = 0; i < p; ++i)array_out[i] = array_in[narray - p + i];
	for(int i = 0; i < nsteps; ++i)array_out[i + p] = inner_product(AR_coeffs.begin(),AR_coeffs.end(), array_out.begin(), 0.f) + AR_coeffs[0];
	array_out.erase(array_out.begin(), array_out.begin() + p);
	return array_out;
}

vector<float> arima::MA_forecast(const vector<float> &array_in, const int &nsteps) {
	unsigned q = ARIMAget_q();
	unsigned narray = array_in.size();
	// calculate epsilon[t] = Y[t] - MA_coeffs*epsilon[t-]
	// MA_coeffis is backward stored for simiplicity
	vector<float> eps(narray, 0.f);
	eps[0] = array_in[0];
	for(int i = 1; i < narray; ++i) {
		if(i < q)eps[i] = array_in[i] - inner_product(MA_coeffs.begin(), MA_coeffs.begin() + i, eps.begin(), 0.f);
		else eps[i] = array_in[i] - inner_product(MA_coeffs.begin(), MA_coeffs.end(), eps.begin(), 0.f);
	}
	// Now we predict nsteps in the future
	vector<float> array_out(narray, 0.f);
	//for(int i = 0; i < nsteps; ++i) {
	//	array_out[i] = 
	//	eps[q - i - 1] = array_in[i + 1] - inner_product(MA_coeffs.begin(), MA_coeffs.begin() + i, eps.end() - i, 0.f);
}
