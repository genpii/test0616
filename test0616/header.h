#pragma once

#include "stdafx.h"

using namespace std;

/*functions*/
//int physio();
void Bsector(const vector<vector<float>>& env, float dangle);
void cairo(const vector<vector<float>>& env, float dangle);
void BSector2(const vector<vector<float>>& env, float dangle, float fs);
float BSector22(const vector<vector<double>>& env, float dangle, float fs, int frm);
//void Bsector3(const vector<vector<float>>& env, float dangle);

/*functions of DSP.cpp*/
void writev(vector<float> &v, string &str);
void writev(vector<float> &v, float &sc, string &str);
vector<float> NormCC(const vector<float> &x1, const vector<float> &x2, int lag, int offset);
vector<float> pwrspe(const vector<float> &v, const int &order);



/*class of fileopen.cpp*/
class file {
protected:
	ifstream fin;

public:	
	file(string filename);
	file();
	~file();

	void open(string filename);
	void start();
	void warp(int pos);
	void go(int pos);
};


class a10 : public file{

public:
	unsigned short len_record, frame, line, sample, ch;
	char* probe_name = (char*)malloc(8 + 1);
	float frq_probe, max_angle, offset_from_center, rad_of_cuv,
		frq_t, frq_r, frq_s, acq_start, acq_end, range;
	unsigned short probe_type, pole, wave, burst, line_start, line_end, max_beam;
	unsigned short focus_num, PRT;
	float focus_first, FR;
	double RF_size;
	vector<vector<vector<vector<short>>>> RF;
	vector<vector<vector<short>>> RF0;

	//function
	a10();
	a10(string filename);
	~a10();
	
	void loadheader();
	void printheader();
	void loadRF();
	void loadRF0(int frame);
	void freeRF();
	void freeRF0();
	void rmbias();
	short eledat(int frame, int line, int ch, int sample);

	//int plotRF0();
	int plotRF0(string dir);
	vector<vector<float>> calcenv(int frame, float max_angle, float frq_s);
};

/*class of physio.cpp*/
class physio : file{
public:
	string fn;
	vector<short> ECG;
	vector<short> PCG_min;
	vector<short> PCG_max;

	physio(string filename);
	~physio();
	int extract(int offset);

	void write();
};


/*class of GN.cpp*/
class GN{
	vector<vector<float>> jacob;
	vector<vector<float>> jacobt;
	vector<vector<float>> jj;
	vector<vector<float>> jjinv;
	vector<float> beta;
	vector<double> res;
	int m, n;
	Ipp32f *x, *y;

	float tmp;

public:
	GN(const vector<double> &x_ini, const vector<double> &y_ini);
	//~GN();
	void setj();
	void trans();
	void loadpoint();
	void mul(const vector<vector<float>> &a, const vector<vector<float>> &b);
	void mul(const vector<vector<float>> &a, const vector<float> &b);
	vector<double> solve(const double &r_ini, const double &c_ini);
};

class est_ss{
public:
	int grid_w, grid_h;
	int bpg; //beam per grid
	int line, ch, sample;

	const double c0 = 1540.0;
	double angle;
	double fs, f0;

	vector<int> b_ele = { 8, 9, 10, 11, 12, 13, 14, 19, 42, 47, 56, 79, 80 }; // broken element

	vector<double> elex;
	vector<double> theta; //beam oriented
	vector<double> phi; // lateral bound oriented
	vector<double> scatdep;

	vector<vector<vector<short>>> RF;

	
	vector<vector<double>> delay;
	vector<vector<double>> tau;
	vector<vector<double>> speed;
	vector<vector<vector<vector<double>>>> path;
	vector<vector<double>> path_near;
	
	//simulate
	vector<vector<double>> sp_g;
	double sp_n;

	//functions
	est_ss(int w, int h, int beam);
	void set_parameter(double dangle, double frq_s, double frq_t);
	void loadRF(const vector<vector<vector<short>>>& RF0);

	void sim_RF();
	void sim_delay_setspeed();
	void sim_delay_calcdelay();
	int calc_delay(double depth); // depth(mm);
	int calc_delay_ph(double depth); // depth(mm);
	int calc_path();

	int del_brokenelement();

	int write_mat();

	void SVD(); //singular value decomposition of path matrix
	
};

struct Point2D
{
	double x, y;
};