//main.cpp

//author: Gen Onodera
/* This C++ program is for the element RF data acquired by Hitachi Aloka alpha-10.
	<source>
	fileopen.cpp - open and save the binary data including signal and other information
	physio.cpp - open DICOM data in order to extract physio data
	Bsector_cairo.cpp - draw sector B-mode image using cairo graphic library
	Bsector_OpenCV.cpp - draw sector B-mode image using OpenCV
	DSP.cpp - functions of digital signal processing using Intel IPP
	<header>
	stdafx.h - include C++ STL
	header.h - header of hand-made classes and fuctions
	graphics.h - reprodution of libXG library using Windows API
	matplotlibcpp.h - plot graph data using matplotlib for Python

	*/

#include "stdafx.h"
#include "header.h"
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

int _tmain(int argc, _TCHAR* argv[])
{
	Py_SetPythonHome("C:/Python27/");
	//ダイナミックフォーカスを行うかどうか
	bool dynamic_f = true;



	//string RFdir = "D:/RFdata/study/20160622/";
	//string dirname = "D:/RFdata/study/20151026/";
	//physio phy(dirname + "2");
	//int nn = 9600;
	//vector<double> x(nn), y(nn);

	//for (int i = 0; i < nn; ++i)
	//	x[i] = i - nn / 2;
	//for (int i = 1; i <= 5; ++i){
	//	double depp = i * 20000;
	//	for (int j = 0; j < nn; ++j){
	//		y[j] = (depp + sqrt(pow(depp, 2) + pow(x[j], 2))) / 1540.0;
	//	}
	//	plt::plot(x, y, "-");
	//}
	////plt::plot(x, y, "--r");
	//plt::show();

#ifdef _DEBUG
	cout << "debugging now!\n";
#endif
	//CrossCorrNormExample();
	/* open data */
	cout << "Load started.\n";

	//よくエラーおきる(Releaseビルドで)
	//char *RFdirchar = "D:/RFdata/study/20160617/2.crf";
	//string *RFdir = new string(RFdirchar);

	//string RFdir = "D:/RFdata/study/20160617/2.crf";
	//string *RFdir = new string();
	//const string RFdir = "D:/RFdata/study/20160617/2.crf";
	//string RFdir("D:/RFdata/study/20160617/2.crf");
	//RFdir = "52101_1.crf"; //X220
	//a10 raw(*RFdir);
	//a10 raw("D:/RFdata/study/20160729/1.crf");
	a10 raw("./20160729/1.crf");


	raw.frq_s = 30.0; //30MHzにしておく
	raw.printheader();
	unsigned short frame = raw.frame;
	unsigned short line = raw.line;
	unsigned short sample = raw.sample;
	unsigned short ch = raw.ch;
	float max_angle = raw.max_angle;
	float frq_t = raw.frq_t;
	float frq_r = raw.frq_r;
	float frq_s = raw.frq_s;
	float FR = raw.FR;
	int physio_offset = 1000 / FR;
	float pitch;

	vector<int> b_ele; //故障した素子 <-後でクラスa10の方に組み込む予定
	string p_name = raw.probe_name;
	if (p_name == "52101"){ //52101だと端が使える
		b_ele = {};
		pitch = 200;
	}
	else if (p_name == "52105"){ //52105だと使えない
		b_ele = { 0, 1, 2, 3, 4, 5, 6, 7, 88, 89, 90, 91, 92, 93, 94, 95 };
		pitch = 240;
	}
	else{
		b_ele = {};
		pitch = 200;
	}
	int N_b_ele = b_ele.size();
	int finest_ele = 48; //基準素子
	vector<int> fine_ele(raw.ch, 0);
	for (int i = 0; i < raw.ch; ++i)
		fine_ele[i] = i;
	for (auto it = b_ele.begin(); it != b_ele.end(); ++it){
		auto res = remove(fine_ele.begin(), fine_ele.end(), *it);
		fine_ele.erase(res, fine_ele.end());
	}

	//phy.extract(physio_offset);
	//phy.write();

	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;

	//raw.loadRF();

	raw.loadRF0(0);
	//raw.plotRF0("line");
	//調べる箇所を限定する際の変数
	int cline = 60;//71;
	int carea = 2;

	raw.generate_AS(cline);
	
	vector<double> xx(4 * sample), yy(4 * sample);
	vector<float> yyy(4 * sample);
	for (int i = 0; i < 4 * sample; ++i){
		xx[i] = i;
		yy[i] = sqrt(pow(raw.ele0re[47][i], 2) + pow(raw.ele0im[47][i], 2));
		//yy[i] = log(sqrt(pow(raw.ele0re[47][i], 2) + pow(raw.ele0im[47][i], 2)));
		yyy[i] = static_cast<float>(yy[i]);
	}
	HanningMovingAverage(yyy, 21);
	vector<int> peak = PeakDetection(yyy);
	string str = "yyy.dat";
	PlotVector(yyy, str);
	str = "peak.dat";
	ofstream foutt(str, ios_base::out);
	for (int i = 0; i < peak.size(); ++i)
		foutt << peak[i] << " " << yyy[peak[i]] << "\n";
	foutt.close();
	//plt::plot(xx, yyy, "-");
	//plt::plot(peak, yyy, ".");
	//plt::grid(true);
	//plt::show();

	//相関計算に用いる定数変数の設定
	int tryN = 81; //繰り返し回数
	vector<float> cc(tryN, 0); //音速セット
	for (int i = 0; i < tryN; ++i)
		cc[i] = 1540.0 + (i - (tryN - 1) / 2) * 1.25; //[m/s]
	vector<float> xi(ch, 0); // x-coordinate of each element
	for (int i = 0; i < ch; ++i)
		//xi[i] = 0.24 * (39.5 - i) * 1e+3;
		xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	vector<float> theta(line, 0);
	for (int i = 0; i < line; ++i) //beam angle
		theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	int aarea = 4;
	vector<vector<int>> apf(line, vector<int>(aarea, 0)); //基準素子におけるラインごとの解析点
	
	int dep_interval = 1000;
	float dep, time_r, point_r_float, time_r2, point_r_float2;
	int point_r, point_r2;
	vector<vector<float>> qof(tryN, vector<float>(dep_interval, 0));
	vector<vector<float>> qof2(tryN, vector<float>(dep_interval, 0));
	ofstream fout("qof.dat", ios_base::out);
	ofstream fout2("qof2.dat", ios_base::out);
	for (int i = 0; i < tryN; ++i){
		for (int j = 0; j < dep_interval; ++j){
			qof[i][j] = 0.0;
			dep = static_cast<float>(j)* sample / frq_s / dep_interval; //us
			for (int k = 0; k < fine_ele.size(); ++k){
				time_r = (dep + sqrt(pow(dep, 2) + 4.0 * xi[fine_ele[k]] / cc[i] * (xi[fine_ele[k]] / cc[i] - dep * sin(theta[cline])))) / 2.0;
				
				point_r_float = time_r * (4.0 * frq_s);

				point_r = static_cast<int>(point_r_float) + 1;

				if (point_r - point_r_float < 0.5)
					--point_r;
				if (!(point_r < 0) && point_r < 4 * sample){
					qof[i][j] += (sqrt(pow(raw.ele0re[fine_ele[k]][point_r], 2) + pow(raw.ele0im[fine_ele[k]][point_r], 2)));
				}

			}
			if (dep >= 10.0){
				fout << cc[i] << " " << dep << " " << qof[i][j] << "\n";
			}
		}
		fout << "\n";
		fout2 << "\n";
	}
	cout << "finish!!!!!\n";
	fout.close();
	fout2.close();

	return 0;
}


