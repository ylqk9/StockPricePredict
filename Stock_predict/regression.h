#include "stdafx.h"
#include "dataload.h"

//mean normalization for stock
void mean_normalization(const stock_features& features);

//regression, mode = r (row-major), mode = c (col-major)
float             *regression(const char &mode, float* x, float *y, const int &nx, const int &ny);
std::vector<float> regression(const char &mode, const std::vector<std::vector<float> > &x, const std::vector<float> &y, const int &nx, const int &ny);