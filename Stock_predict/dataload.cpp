#include "stdafx.h"
#include "dataload.h"
using namespace std;

stock_values::stock_values(const int& value_number, const int& days_observ_number, const int& time_observ_number) :
nvalues(value_number), ndays(days_observ_number), ntimes(time_observ_number) {
	// value(stock_ind, day_ind, time_ind)
	values = new float**[nvalues];
	for (int i = 0; i < nvalues; ++i) values[i] = new float*[ndays];
	values[0][0]= new float[nvalues*ndays*ntimes];
	for (int i = 0; i < nvalues; ++i) {
		for (int j = 0; j < ndays; ++j) {
			values[i][j] = i*ndays*ntimes + j*ntimes + values[0][0];
		}
	}
	// values_time_major(day_ind, time_ind, stock_ind)
	values_time_major = new float**[ndays];
	for (int i = 0; i < ndays; ++i) values_time_major[i] = new float*[ntimes];
	values_time_major[0][0]= new float[nvalues*ndays*ntimes];
	for (int i = 0; i < ndays; ++i) {
		for (int j = 0; j < ntimes; ++j) {
			values_time_major[i][j] = i*ntimes*nvalues + j*nvalues + values_time_major[0][0];
		}
	}
}

stock_values::~stock_values() {
	delete [] values[0][0];
	for (int i = 0; i < nvalues; ++i) delete [] values[i];
	delete [] values;
	delete [] values_time_major[0][0];
	for (int i = 0; i < ndays; ++i) delete [] values_time_major[i];
	delete [] values_time_major;
}

stock_features::stock_features(const int& feature_number, const int& days_observ_number, const int& time_observ_number) :
nfeatures(feature_number), ndays(days_observ_number), ntimes(time_observ_number) {
	// feature(day_ind, time_ind, feature_ind);
	features    = new float**[ndays];
	for (int i = 0; i < ndays; ++i) features[i] = new float*[ntimes];
	features[0][0] = new float[ndays*ntimes*nfeatures];
	for (int i = 0; i < ndays; ++i) {
		for (int j = 0; j < ntimes; ++j) {
			features[i][j] = i*nfeatures*ntimes + j*nfeatures + features[0][0];
		}
	}
}

stock_features::~stock_features() {
	delete [] features[0][0];
	for (int i = 0; i < ndays; ++i) delete [] features[i];
	delete [] features;
}

StockData::StockData(const int& stock_number, const int& feature_number, const int& time_observ_number, const int& days_observ_number) : 
nstocks(stock_number), ntimes(time_observ_number), nfeatures(feature_number), ndays(days_observ_number),
value(stock_number, days_observ_number, time_observ_number), feature(feature_number, days_observ_number, time_observ_number){
}

StockData::~StockData() {
}

void StockData::FileList () {
	files.resize(0);
	string intsink;
	for (int i = 0; i < ndays; ++i) files.push_back(to_string((_Longlong)(i + 1)) + ".csv");
}

void StockData::Load (const string& file_name, const int& file_ind) {
	string databuf, cell;
	ifstream DataFile(file_name.c_str(), ifstream::in);
	DataFile.ignore(numeric_limits<streamsize>::max(), '\n' );
	for (int i = 0; i < ntimes; ++i) {
		getline(DataFile, databuf, '\n');
		int j = 0;
		stringstream linebuf(databuf);
		while (getline(linebuf, cell, ',')) {
			// value(stock_ind, day_ind, time_ind)
			// values_time_major(day_ind, time_ind, stock_ind)
			// feature(day_ind, time_ind, feature_ind);
			if (j < nstocks) {
				value.values[j][file_ind][i] = stof(cell);   //stof is in C++ 11. Here it could be a MS extension.
				value.values_time_major[file_ind][i][j] = stof(cell);
			} else feature.features[file_ind][i][j - nstocks] = stof(cell);
			++j;
		}
	}
	DataFile.close();
	return;
}

void StockData::BatchLoad () {
	FileList();
	//for (int i = 0; i < ndays; ++i) cout << files[i] << endl;
	for (int i = 0; i < ndays; ++i)	{
		Load (files[i], i);
	}
}