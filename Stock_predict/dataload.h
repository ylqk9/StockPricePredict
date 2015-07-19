#include "stdafx.h"

#ifndef _DATALOAD_H_
#define _DATALOAD_H_

class stock_values
{
public:
	float ***values;
	float ***values_time_major;
	int   nvalues;
	int   ndays;
	int   ntimes;
	stock_values(const int& value_number, const int& days_observ_number, const int& time_observ_number);
	~stock_values();
};

class stock_features
{
public:
	float ***features;
	int   nfeatures;
	int   ndays;
	int   ntimes;
	stock_features(const int& feature_number, const int& days_observ_number, const int& time_observ_number);
	~stock_features();
};

class StockData
{
public:
	stock_values   value;
	stock_features feature;
	StockData (const int& stock_number, const int& feature_number, const int& time_observ_number, const int& days_observ_number);
	~StockData ();
	void Load (const std::string& file_name, const int& file_ind);
	void FileList ();
	void BatchLoad ();
	void Dbg_array ();
private:
	int nfeatures, nstocks, ntimes, ndays;
	int load_success;
	float *value_store;
	float *feature_store;
	std::vector<std::string> files;
};

#endif