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
using namespace Eigen;
namespace plt = matplotlibcpp;

int _tmain(int argc, _TCHAR* argv[])
{
	Py_SetPythonHome("C:/Python27/");

	/* open data */
	cout << "Load started.\n";
	//ファイル名絡みでよくエラーおきるので注意(Releaseビルドで)
	//a10 raw("D:/RFdata/study/20160729/1.crf");
	//a10 raw("D:/RFdata/study/20160819/1.crf");
	a10 raw("D:/RFdata/study/20160823/1.crf");
	if (raw.probe_type == 3)
		raw.frq_s = 30.0; //30MHzにしておく
	raw.printheader();
	//a10 raw("./20160729/1.crf");
	//a10 raw("D:/RFdata/study/20160622/52101_1.crf");

	raw.loadRF0(0);
	//BLinear(raw);
	//簡単のためmain内の変数にする
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


	//raw.plotRF0(0);
	vector<int> b_ele; //故障した素子 <-後でクラスa10の方に組み込む予定
	string p_name = raw.probe_name;

	int N_b_ele = b_ele.size();
	int finest_ele = 48; //基準素子
	vector<int> fine_ele(raw.ch, 0);
	for (int i = 0; i < raw.ch; ++i)
		fine_ele[i] = i;
	for (auto it = b_ele.begin(); it != b_ele.end(); ++it){
		auto res = remove(fine_ele.begin(), fine_ele.end(), *it);
		fine_ele.erase(res, fine_ele.end());
	}
	//vector<float> xi(ch, 0); // x-coordinate of each element
	vector<float> xi = raw.getxi();

	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;

	//raw.loadRF();
	//raw.plotRF0("line");
	//調べる箇所を限定する際の変数
	int cline = 90;//71;
	int carea = 2;

	raw.generate_AS(cline);

	vector<double> xx(4 * sample), yy(4 * sample);
	vector<float> yyy(4 * sample);
	int ref = 70;
	for (int i = 0; i < 4 * sample; ++i){
		xx[i] = i;
		//yy[i] = raw.RF0[cline][47][i];
		//yy[i] = raw.ele0re[47][i];
		yy[i] = sqrt(pow(raw.ele0re[ref][i], 2) + pow(raw.ele0im[ref][i], 2));
		//yy[i] = log(sqrt(pow(raw.ele0re[47][i], 2) + pow(raw.ele0im[47][i], 2)));
		yyy[i] = static_cast<float>(yy[i]);
		//yyy[i] = raw.ele0im[47][i];
	}

	

	HanningMovingAverage(yyy, 601); //移動平均 引数は窓幅

	plt::plot(xx, yyy);
	plt::grid(true);
	plt::show();

	vector<int> peak = PeakDetection(yyy, 2048, 12000); //2048-12000以外のピークを除去

	vector<pair<float, int>> peakamp(peak.size()); //peak=(振幅,index)
	for (int i = 0; i < peak.size(); ++i){
		peakamp[i] = make_pair(yyy[peak[i]], peak[i]);
	}
	sort(peakamp.begin(), peakamp.end(), greater<pair<float, int>>()); //振幅順にソート

	//今だけ ゴミピーク消去
	/*peakamp.erase(peakamp.begin() + 4);
	peakamp.erase(peakamp.begin() + 2);
	peakamp.erase(peakamp.begin());*/
	//peakamp.erase(peakamp.begin() + 5);
	peakamp.erase(peakamp.begin() + 5);
	peakamp.erase(peakamp.begin() + 4);
	peakamp.erase(peakamp.begin() + 1);
	peakamp.erase(peakamp.begin());


	//

	int N_est = 3; //推定ピーク数
	vector<tuple<float, float, float>> val_est(N_est);
	ofstream fo("n.dat", ios_base::out);
	//ch = 48;
	//vector<int> depro = MakeDelayProfile(raw.ele0re, raw.ele0im, peakamp[1].second);


	cout << endl;
	for (int l = 0; l < N_est; ++l){
		vector<int> depro = MakeDelayProfile(raw.ele0re, raw.ele0im, peakamp[l].second);
		Matrix3f M;
		Vector3f v;
		float tmpmv;
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < 3; ++j){
				tmpmv = 0.0;
				for (int k = 0; k < ch; ++k)
					tmpmv += pow(xi[k], 4 - i - j);
				M(i, j) = tmpmv;
			}
			tmpmv = 0.0;
			for (int j = 0; j < ch; ++j)
				tmpmv += pow(xi[j], 2 - i) * pow(static_cast<float>((peakamp[l].second + depro[47]) / 2 + depro[j]) / (4 * frq_s), 2);
			v(i) = tmpmv;
		}
		Vector3f ans = M.inverse() * v;
		//FullPivLU<MatrixXf> LU(M);
		//cout << "LU Rank: " << LU.rank() << endl;
		//Vector3f ans = LU.solve(v);
		float SS = 1 / sqrt(ans(0));
		float DEP = sqrt(ans(2) / ans(0));
		float thetaR = asin(-ans(1) / 2.0 / sqrt(ans(2) * ans(0)));
		cout << "SV: " << SS << " m/s" <<endl;
		cout << "DEP: " << DEP / 1000.0 << " mm" << endl;
		cout << "thetaR: " << thetaR * 180.0 / M_PI << " deg" << endl;
		cout << "lat pos: " << DEP * sin(thetaR) << " um" << endl;
		cout << endl;

		val_est[l] = make_tuple(DEP, SS, thetaR);
	}

	fo << get<1>(val_est[1]) << " " << get<0>(val_est[1]) << "\n";


	fo.close();
	sort(val_est.begin(), val_est.end());

	int M = 3;
	int N = 3;
	int ll = 15000; //[um]
	MatrixXf D(M, N);
	D = MatrixXf::Zero(M, N);
	if (N == 3)
		D << get<0>(val_est[0]), 0, 0,
			2 * ll, get<0>(val_est[1]) - 2 * ll, 0,
			2 * ll, ll, get<0>(val_est[2]) - 3 * ll;
	else if(N == 4)
		D << ll, get<0>(val_est[0]) - ll, 0, 0,
			ll, ll, get<0>(val_est[1]) - 2 * ll, 0,
			ll, ll, ll, get<0>(val_est[2]) - 3 * ll;
	VectorXf T(M);
	T = VectorXf::Zero(M);
	for (int i = 0; i < M; ++i)
		T(i) = get<0>(val_est[i]) / get<1>(val_est[i]);
	
	cout << "D is:" << endl << D << endl;
	cout << "T is:" << endl << T << endl;


	JacobiSVD<MatrixXf> SVD(D, ComputeThinU | ComputeThinV);
	//cout << SVD.singularValues() << endl;
	VectorXf tau(N);
	tau = SVD.solve(T);
	//tau = D.inverse() * T;
	VectorXf SV(N);
	for (int i = 0; i < N; ++i)
		SV(i) = 1 / tau(i);
	
	cout << "SV is:" << endl << SV << endl;

	cout << "finished.\n" << "\n";
	return 0;
}


