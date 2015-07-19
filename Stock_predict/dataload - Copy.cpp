#include "stdafx.h"
#include "dataload.h"
using namespace std;

StockData::StockData(const int& stock_number, const int& feature_number, const int& time_observ_number) : 
nstock(stock_number), ntime(time_observ_number), nfeature(feature_number) {
	value   = new float[nstock*ntime];
	feature = new float[nfeature*ntime];
	//for (int i = 1; i < nstock; ++i) value[i]   = i*nstock   + value[0];
	//for (int i = 1; i < nfeature; ++i) feature[i] = i*nfeature + feature[0];
}

StockData::~StockData() {
	delete [] value;
	delete [] feature;

}

void StockData::Load (const string& file_name) {
	string databuf, cell;
	ifstream DataFile(file_name.c_str(), ifstream::in);
	DataFile.ignore(numeric_limits<streamsize>::max(), '\n' );
	for (int i = 0; i < 1; ++i) {
		getline(DataFile, databuf, '\n');
		int j = 0;
		stringstream linebuf(databuf);
		while (getline(linebuf, cell, ',')) {
			if (j < nstock) value[i*nstock + j] = stof(cell); //stof is in C++ 11.
			else feature[i*nfeature + j - nstock] = stof(cell); //Here it could be a MS extension.
			++j;
		}
	}
	for (int i = 0; i < nstock; ++i) cout << value[i] <<endl;
	cout << "------------------------------------" << endl;
	for (int i = 0; i < nfeature; ++i) cout << feature[i] <<endl;
	return;
}