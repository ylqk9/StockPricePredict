#include "stdafx.h"
#include "dataload.h"
#include "regression.h"
#include <mkl.h>

void mean_normalization(const stock_features& features) {
	// feature(day_ind, time_ind, feature_ind);
	float *sum_feature = new float[features.nfeatures]();
	float *avg_feature = new float[features.nfeatures]();
	float *var_feature = new float[features.nfeatures]();
	float *sd_feature = new float[features.nfeatures]();
	float *scale_feature = new float[features.nfeatures]();
	float *shift_feature = new float[features.nfeatures]();
	float *result = new float[features.nfeatures]();
	//blas input prameters
	int   N(features.nfeatures);
	float a(1.f);
	int   incX(1), incY(1);
	//
	for(int i = 0; i < features.ndays; ++i) {
		for(int j = 0; j < features.ntimes; ++j) {
			for(int k = 0; k < features.nfeatures; ++k) {
				sum_feature[k] += features.features[i][j][k];
			}
		}
	}
	for(int k = 0; k < features.nfeatures; ++k) avg_feature[k] = sum_feature[k]/(features.ntimes)/(features.ndays);
	for(int i = 0; i < features.ndays; ++i) {
		for(int j = 0; j < features.ntimes; ++j) {
			for(int k = 0; k < features.nfeatures; ++k) {
				var_feature[k] += pow((features.features[i][j][k] - avg_feature[k]),2);
			}
		}
	}
	//(x-x_bar)/x_sd = x*(1/x_sd) + (-x_bar/x_sd) = x*scale + shift
	for(int k = 0; k < features.nfeatures; ++k) {
		var_feature[k] = var_feature[k]/(features.ntimes)/(features.ndays);
		sd_feature[k]  = sqrt(var_feature[k]);
		scale_feature[k] = 1.f/sd_feature[k];
		shift_feature[k] = - avg_feature[k]/sd_feature[k];
	}
	std::string file_name = "NormalizedFeature.csv";
	std::ofstream NormalizedData;
	NormalizedData.open(file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
	for(int i = 0; i < features.ndays; ++i) {
		for(int j = 0; j < features.ntimes; ++j) {
			vsMul (N, features.features[i][j], scale_feature, result);
			cblas_saxpy(N, a, shift_feature, incX, result, incY);
			cblas_scopy(N, result, incX, features.features[i][j], incY);
			for(int k = 0; k < (features.nfeatures - 1); ++k) NormalizedData << result[k] << ",";
			NormalizedData << result[features.nfeatures - 1] << std::endl;
		}
	}
	//for(int k = 0; k < features.nfeatures; ++k) std::cout << sd_feature[k] << std::endl;
	//std::cout << "-------------------------" << std::endl;
	//for(int k = 0; k < features.nfeatures; ++k) std::cout << avg_feature[k] << std::endl;
	NormalizedData.close();
	return;
}

float *regression(const char &mode, float* x, float *y, const int &nx, const int &ny){
	//x(ny, nx); y(ny)
	//
	// regression will perform (x'x)^(-1)*(x'y)
	float **xTx = new float*[nx];
	xTx[0] = new float[nx*nx]();
	for(int i = 1; i < nx; ++i)xTx[i] = i*nx + xTx[0];
	float *xTy = new float[ny]();
	for(int i = 0; i < nx; ++i) {
		for(int j = i; j < nx; ++j) {
			if(mode == 'c')	xTx[i][j] = cblas_sdot(ny, (x + i*ny), 1, (x + j*ny), 1);
			if(mode == 'r') xTx[i][j] = cblas_sdot(ny, (x + i), nx, (x + j), nx);
				//std::inner_product((x + i*ny), (x + ny + i*ny), (x + j*ny), 0.f); inner_product works differently in MSVS than STL
			xTx[j][i] = xTx[i][j];
		}
	}
	for(int i = 0; i < nx; ++i) {
		if(mode == 'c')xTy[i] = cblas_sdot(ny, (x + i*ny), 1, y, 1);
		if(mode == 'r')xTy[i] = cblas_sdot(ny, (x + i), nx, y, 1);
	}
		//std::inner_product((x + i*ny), (x + ny + i*ny), y, 0.f);
	int info;
	int *ipiv = new int[nx];
	info = LAPACKE_sgesv(LAPACK_COL_MAJOR, nx, 1, xTx[0], nx, ipiv, xTy, nx);
	return xTy;
}

//void multi_reg_single_day(const stock_features& features, const stock_values &value, const int& day_ind){
//	// value(stock_ind, day_ind, time_ind)
//	// feature(day_ind, time_ind, feature_ind);
//	//BETA = (xTy)/(xTx)
//	int   m(features.ntimes), n(features.nfeatures), nrhs(1);
//	int   lda(features.nfeatures), ldb(1), lwork(-1), info;
//	float *a = new float[features.ntimes*features.nfeatures];
//	for(int i = 0; i < features.ntimes*features.nfeatures; ++i) a[i] = *(features.features[day_ind][0]+i);
//	float *b = new float[value.ntimes];
//	for(int i = 0; i < features.ntimes; ++i) b[i] = value.values[0][day_ind][i];
//	float *reg_res = new float[m]();
//	for(int i = 0; i < m; ++i) {
//		std::cout << reg_res[i] << std::endl;
//	}
//
//	float wkopt;
//	float* work;
//
//	//for(int i = 0; i < m; ++i) {
//	//	for(int j = 0; j < n; ++j) {
//	//		std::cout << *(a + i*m + j) << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//for(int i = 0; i < m; ++i) std::cout << b[i] << std::endl;
//
//	
//	//info = LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, a, lda, b, ldb);
//
//
//
//	char mode('T');
//	sgels( &mode, &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info );
//	lwork = (int)wkopt;
//	work = (float*)malloc( lwork*sizeof(float) );
//	sgels( &mode, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info );
//
//
//
//	//lwork = (int)wkopt;
//	//work = (float*)malloc( lwork*sizeof(float) );
//	//sgels("No transpose", &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info);
//	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, 1.f, features.features[day_ind][0], lda, b, ldb, 0.f, reg_res, 1);
//}
