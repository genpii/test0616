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
IppStatus CrossCorrNormExample(void);

int _tmain(int argc, _TCHAR* argv[])
{

	//ダイナミックフォーカスを行うかどうか
	bool dynamic_f = true;

	//string RFdir = "D:/RFdata/study/20160622/";
	//string dirname = "D:/RFdata/study/20151026/";
	//physio phy(dirname + "2");
	/*int nn = 5000;
	vector<double> x(nn), y(nn);
	for (int i = 0; i < nn; ++i){
	x[i] = i;
	y[i] = sin(2 * M_PI * i / 360.0);
	}
	plt::plot(x, y, "--r");
	plt::show();*/



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
	a10 raw("D:/RFdata/study/20160617/2.crf");

	vector<int> b_ele; //故障した素子 <-後でクラスa10の方に組み込む予定
	string p_name = raw.probe_name;
	if (p_name == "52101") //52101だと端が使える
		b_ele = { 8, 9, 10, 11, 12, 13, 14, 25, 41, 42, 47, 56, 79, 80, 88 };
	else if (p_name == "52105") //52105だと使えない
		b_ele = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 41, 42, 47, 56, 79, 80, 88, 89, 90, 91, 92, 93, 94, 95 };
	else
		b_ele = {};
	int N_b_ele = b_ele.size();
	int finest_ele = 48; //基準素子

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
	//phy.extract(physio_offset);
	//phy.write();

	/* load RF data */
	// RF[frame][line][ch][sample]
	cout << "initializing array...\n";
	short tmp = 0;

	//raw.loadRF();

	raw.loadRF0(0);
	

	//vector<vector<float>> env(line, vector<float>(sample, 0));
	//env = raw.calcenv(0, max_angle, frq_s);
	//BSector2(env, max_angle, frq_s);
	//raw.plotRF0("0627phantom");


	cout << "creating analytic signal...\n";

	vector<vector<vector<float>>> elere(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0)));
	vector<vector<vector<float>>> eleim(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0)));

	

	//調べる箇所を限定する際の変数
	int cline = 71;
	int carea = 2;


	//spec and buffer setting for FFT
	Ipp8u *specbuff, *initbuff, *workbuff;
	Ipp8u *specbufi, *initbufi, *workbufi;
	int size_specf, size_initf, size_workf;
	int size_speci, size_initi, size_worki;
	IppsFFTSpec_C_32fc *specf = 0;
	IppsFFTSpec_C_32fc *speci = 0;
	Ipp32fc *ipsrc = ippsMalloc_32fc((int)sample);
	Ipp32fc *ipdst = ippsMalloc_32fc((int)sample);
	Ipp32fc *ipsrc2 = ippsMalloc_32fc((int)(4 * sample));
	Ipp32fc *ipdst2 = ippsMalloc_32fc((int)(4 * sample));
	const int fftorder = (int)(log((double)sample) / log(2.0));
	const int ifftorder = (int)(log((double)(4 * sample)) / log(2.0));
	ippsFFTGetSize_C_32fc(fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_specf, &size_initf, &size_workf);
	ippsFFTGetSize_C_32fc(ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_speci, &size_initi, &size_worki);
	specbuff = ippsMalloc_8u(size_specf);
	specbufi = ippsMalloc_8u(size_speci);
	initbuff = ippsMalloc_8u(size_initf);
	initbufi = ippsMalloc_8u(size_initi);
	workbuff = ippsMalloc_8u(size_workf);
	workbufi = ippsMalloc_8u(size_worki);
	ippsFFTInit_C_32fc(&specf, fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuff, initbuff);
	ippsFFTInit_C_32fc(&speci, ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbufi, initbufi);
	//For non-dynamic focus
	Ipp32fc *ipdelay = ippsMalloc_32fc((int)(sample));
	float min_delay;
	


	//相関計算に用いる定数変数の設定
	int tryN = 11; //繰り返し回数
	vector<float> cc(tryN, 0); //音速セット
	for (int i = 0; i < tryN; ++i)
		cc[i] = 1540.0 + (i - (tryN - 1) / 2) * 10.0; //[m/s]
	vector<float> xi(ch, 0); // x-coordinate of each element
	for (int i = 0; i < ch; ++i)
		xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	vector<float> theta(line, 0);
	for (int i = 0; i < line; ++i) //beam angle
		theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	int aarea = 4;
	vector<vector<int>> apf(line, vector<int>(aarea, 0)); //基準素子におけるラインごとの解析点
	
	float dep; //解析深さ[um]
	int ap; //解析点
	float apdecimal; //小数あり
	bool for_ini; //開始素子とそれ以外を判別する
	int next_ele; //生きている隣接素子番号
	vector<int> fine_ele(ch, 0);
	for (int i = 0; i < ch; ++i)
		fine_ele[i] = i;
	for (auto it = b_ele.begin(); it != b_ele.end(); ++it){
		auto res = remove(fine_ele.begin(), fine_ele.end(), *it);
		fine_ele.erase(res, fine_ele.end());
	}



	//相関計算に用いる 相関はNCC(not ZNCC)
	IppStatus status;
	const int csrc1Len = 256;
	const int csrc2Len = 4 * sample;
	int lowLag = -256; // -lowLag * 2 + 1 が全遅延量
	const int cdstLen = -2 * lowLag + 1;
	int Lag;
	Ipp32fc *csrc1 = ippsMalloc_32fc(csrc1Len);
	Ipp64f csrc1Norm;
	Ipp32fc *csrc2 = ippsMalloc_32fc(csrc2Len);
	Ipp64f csrc2Normtmp;
	Ipp64f *csrc2Norm = ippsMalloc_64f(cdstLen);
	Ipp32f *csrcNorm = ippsMalloc_32f(cdstLen);
	Ipp32fc *csrcNormCplx = ippsMalloc_32fc(cdstLen);
	Ipp32f *cdstLenzeros = ippsMalloc_32f(cdstLen);
	ippsZero_32f(cdstLenzeros, cdstLen);
	Ipp32fc *cdst = ippsMalloc_32fc(cdstLen);
	Ipp32f *cdstMag = ippsMalloc_32f(cdstLen);
	Ipp32f *cdstPha = ippsMalloc_32f(cdstLen);
	Ipp32f cmax;
	int cmaxidx;
	IppEnum NormA = (IppEnum)(ippAlgAuto | ippsNormNone);
	int bufsize = 0;
	Ipp8u *pbuffer;
	status = ippsCrossCorrNormGetBufferSize(csrc1Len, csrc2Len, cdstLen, lowLag, ipp32fc, NormA, &bufsize);
	if (status != ippStsNoErr)
		return status;
	pbuffer = ippsMalloc_8u(bufsize);
	

	vector<vector<float>> ct(line, vector<float>(ch, 0));
	
	float ctmax = sqrt(pow(raw.focus_first * 1e+3, 2) + pow(xi[0], 2) - 2 * raw.focus_first * 1e+3 * xi[0] * sin(theta[line - 1]));
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch; ++j){
			ct[i][j] = ctmax - sqrt(pow(raw.focus_first * 1e+3, 2) + pow(xi[j], 2) - 2 * raw.focus_first * 1e+3 * xi[j] * sin(theta[i]));
		}
	}
	//ct[0][ch - 1] = 0.0;
	//ct[line - 1][0] = 0.0;


	vector<float> x(4 * sample), y(4 * sample);
	for (int i = 0; i < 4 * sample; ++i)
		x[i] = (float)i / (4 * frq_s);
	//cout << "[";
	//float anacount = 0.0;
	/*for (int j = 0; j < line; ++j){*/
	if (dynamic_f){
		for (int j = cline; j < cline + 1; ++j){

			for (int k = 0; k < ch; ++k){
				ippsZero_32fc(ipsrc, sample);
				ippsZero_32fc(ipdst, sample);
				ippsZero_32fc(ipsrc2, 4 * sample);
				ippsZero_32fc(ipdst2, 4 * sample);
				//set
				for (int l = 0; l < sample - 1; ++l)
					ipsrc[l].re = raw.RF0[j][k][l];
				//ipsrc[l].re = raw.RF[60][j][k][l] - raw.RF[59][j][k][l];

				//do FFT
				ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
				ippsZero_8u(workbuff, size_workf);
				//double positive part and delete negative part
				for (int l = 0; l < sample / 2; ++l){
					ipdst[l].re = ipdst[l].re * 2 / sample;
					ipdst[l].im = ipdst[l].im * 2 / sample;
					ipdst[l + sample / 2].re = 0.0;
					ipdst[l + sample / 2].im = 0.0;
				}
				for (int l = 0; l < 34; ++l){
					ipdst[l].re = 0.0;
					ipdst[l].im = 0.0;
				}


				//ipdst[0].re /= 2;
				//ipdst[0].im /= 2;

				for (int l = 0; l < sample; ++l){
					ipsrc2[l].re = ipdst[l].re;
					ipsrc2[l].im = ipdst[l].im;
				}


				//do IFFT
				ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
				ippsZero_8u(workbufi, size_worki);

				//save
				for (int l = 0; l < 4 * sample; ++l){
					elere[j][k][l] = 4 * sample * ipdst2[l].re;
					eleim[j][k][l] = 4 * sample * ipdst2[l].im;
				}
				ippsZero_32fc(ipdst2, 4 * sample);
			}
		}

		//free IPP array
		ippsFree(ipsrc);
		ippsFree(ipdst);
		ippsFree(ipdst2);
		//free RF
		raw.freeRF0();
		cout << "finished creating AS!\n";
	}
	else{
		for (int j = cline; j < cline + 1; ++j){

			for (int k = 0; k < ch; ++k){
				ippsZero_32fc(ipsrc, sample);
				ippsZero_32fc(ipdst, sample);
				





				//set
				for (int l = 0; l < sample - 1; ++l)
					ipsrc[l].re = raw.RF0[j][k][l];
				//ipsrc[l].re = raw.RF[60][j][k][l] - raw.RF[59][j][k][l];

				//do FFT
				ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
				ippsZero_8u(workbuff, size_workf);
				//double positive part and delete negative part
				for (int l = 0; l < sample / 2; ++l){
					ipdst[l].re = ipdst[l].re * 2 / sample;
					ipdst[l].im = ipdst[l].im * 2 / sample;
					ipdst[l + sample / 2].re = 0.0;
					ipdst[l + sample / 2].im = 0.0;
				}
				for (int l = 0; l < 34; ++l){
					ipdst[l].re = 0.0;
					ipdst[l].im = 0.0;
				}

				for (int l = 0; l < tryN; ++l){
					ippsZero_32fc(ipdelay, sample);
					ippsZero_32fc(ipsrc2, 4 * sample);
					ippsZero_32fc(ipdst2, 4 * sample);
					for (int m = 0; m < sample; ++m){
						ipdelay[m].re = cos(2 * M_PI * ct[j][k] / cc[l] * frq_s * m / sample);
						ipdelay[m].im = -sin(2 * M_PI * ct[j][k] / cc[l] * frq_s * m / sample);
					}
					ippsMul_32fc(ipdst, ipdelay, ipsrc2, sample);

					//do IFFT
					ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
					ippsZero_8u(workbufi, size_worki);
					for (int m = 0; m < 4 * sample; ++m)
						y[m] = ipdst2[m].re;
					if (l == 0 && (k == 28 || k == 48 || k == 68)){
						ostringstream ost;
						ost << "ch:" << k;
						plt::named_plot(ost.str(), x, y, "-");
						plt::legend();
						ost.str("");
						
					}

				}
				ippsZero_32fc(ipdst2, 4 * sample);
			}
		}
	}
	//plt::plot(x, elere[76][48], "-");
	//plt::grid(true);
	//plt::show();



	//plt::ylim(-400, 400);
	////plt::xlim(9900, 10200);
	////plt::grid(true);
	//plt::show();

	
	



	//plot analitic signal
	/*string call;
	vector<float> smpx(4 * sample);
	for (int i = 0; i < 4 * sample; ++i)
		smpx[i] = i;
	while (1){
		cout << "input line number\n";
		cin >> call;
		if (isdigit(call[0])){
			int callline;
			stringstream(call) >> callline;
			if (callline >= 0 && callline <= line - 1){
				cout << "input ch number\n";
				cin >> call;
				if (isdigit(call[0])){
					int callch;
					stringstream(call) >> callch;
					if (callch >= 0 && callch <= ch - 1){
						plt::plot(smpx, elere[callline][callch], "-");
						cin >> callch;
						plt::plot(smpx, elere[callline][callch], "-");
						plt::grid(true);
						plt::show();
					}
				}
				else if (call == "q")
					break;
			}
		}
		else if (call == "q")
			break;
	}*/



	int Lagf;
	vector<vector<double>> tof(aarea, vector<double>(fine_ele.size(), 0));
	//領域内の最大値を解析点と考える(ver.1.0)
	//int Lag;
	//if (dynamic_f){
		for (int i = 0; i < line; ++i){
			//解析点の決定
			for (int j = 0; j < aarea; ++j){
				//cout << i << " " << j << "\n";
				auto it = next(elere[i][finest_ele].begin(), j * (4 * sample) / aarea);
				auto maxit = max_element(it, next(it, (4 * sample) / aarea));
				if (j == 0)
					maxit = max_element(next(it, (4 * sample) / aarea / 2), next(it, (4 * sample) / aarea));
				else if (j == aarea - 1)
					maxit = max_element(it, next(it, (4 * sample) / aarea / 2));
				apf[i][j] = distance(elere[i][finest_ele].begin(), maxit);
			}
		}
	//}
	/*estimation distribution of sound speed(main part)*/

	for (int i = cline; i < cline + 1; ++i){
		
		for (int j = 0; j < fine_ele.size(); ++j){

			ippsZero_32fc(csrc2, csrc2Len);
			for (int k = 0; k < 4 * sample; ++k){
				csrc2[k].re = elere[i][fine_ele[j]][k];
				csrc2[k].im = eleim[i][fine_ele[j]][k];
			}
			
			for (int k = carea; k < carea + 1; ++k){
				ippsZero_32fc(csrc1, csrc1Len);
				Lagf = apf[i][k] - csrc1Len / 2 + 1;
				for (int l = 0; l < csrc1Len; ++l){
					csrc1[l].re = elere[i][finest_ele][Lagf + l];
					csrc1[l].im = eleim[i][finest_ele][Lagf + l];
				}
				csrc1Norm = 0.0;
				ippsNorm_L2_32fc64f(csrc1, csrc1Len, &csrc1Norm);

				ippsZero_32fc(cdst, cdstLen);
				ippsZero_32f(cdstMag, cdstLen);
				ippsZero_32f(cdstPha, cdstLen);
				ippsZero_64f(csrc2Norm, cdstLen);
				ippsZero_32f(csrcNorm, cdstLen);
				ippsZero_32fc(csrcNormCplx, cdstLen);
				status = ippsCrossCorrNorm_32fc(csrc1, csrc1Len, csrc2, csrc2Len, cdst, cdstLen, Lagf + lowLag, NormA, pbuffer);
				for (int m = 0; m < cdstLen; ++m){ //src2のノルム
					csrc2Normtmp = 0.0;
					ippsNorm_L2_32fc64f(csrc2 + Lagf + lowLag + m, csrc1Len, &csrc2Normtmp);
					csrc2Norm[m] = csrc2Normtmp;
				}
				ippsMulC_64f_I(csrc1Norm, csrc2Norm, cdstLen); //src1のノルム*src2のノルム
				ippsConvert_64f32f(csrc2Norm, csrcNorm, cdstLen); //double->float
				ippsRealToCplx_32f(csrcNorm, cdstLenzeros, csrcNormCplx, cdstLen); //実数列->複素列
				ippsDiv_32fc_I(csrcNormCplx, cdst, cdstLen); //正規化
				ippsMagnitude_32fc(cdst, cdstMag, cdstLen); //相関振幅
				ippsPhase_32fc(cdst, cdstPha, cdstLen); //相関位相
				ippsMaxIndx_32f(cdstMag, cdstLen, &cmax, &cmaxidx);

				
				tof[k][j] = static_cast<double>(cmaxidx + lowLag + apf[i][k]) / (4 * frq_s);
				cout << "line:" << i << " ch:" << fine_ele[j] << " tof:" << tof[k][j] << " cmax:" << cmax << "\n";

			}
		}

		//最小二乗法で音速・深さを推定




	}
	


	//ここまで

	//ビーム→テスト音速→チャンネル
	//ofstream cfout("corr.dat", ios_base::out);
	vector<vector<vector<float>>> corrre(ch, vector<vector<float>>(tryN, vector<float>(cdstLen, 0)));
	vector<vector<vector<float>>> corrim(ch, vector<vector<float>>(tryN, vector<float>(cdstLen, 0)));
	vector<vector<vector<float>>> corrmag(ch, vector<vector<float>>(tryN, vector<float>(cdstLen, 0)));
	vector<vector<vector<float>>> corrpha(ch, vector<vector<float>>(tryN, vector<float>(cdstLen, 0)));
	//for (int i = 0; i < line; ++i){ //ライン
	for (int i = cline; i < cline + 1; ++i){

		for (int j = 0; j < tryN; ++j){ //音速セット

			//for (int k = 0; k < aarea; ++k){ //エリア
			for (int k = carea; k < carea + 1; ++k){
				//src1セット
				ippsZero_32fc(csrc1, csrc1Len);
				Lagf = apf[i][k] - csrc1Len / 2 + 1;
				for (int l = 0; l < csrc1Len; ++l){
					csrc1[l].re = elere[i][finest_ele][Lagf + l];
					csrc1[l].im = eleim[i][finest_ele][Lagf + l];
				}
				csrc1Norm = 0.0;
				ippsNorm_L2_32fc64f(csrc1, csrc1Len, &csrc1Norm); //src1のノルム

				//深さを判定する
				dep = (pow(cc[5] * apf[i][k] / (4 * frq_s), 2) - pow(xi[finest_ele], 2)) / 2 / (cc[5] * apf[i][k] / (4 * frq_s) - xi[finest_ele] * sin(theta[i]));
				//for_ini = false;
				for (int l = 0; l < fine_ele.size(); ++l){
					if (l != finest_ele){
						//cout << "(" << i << " " << j << " " << k << " " << l << ")\n";
						ippsZero_32fc(csrc2, csrc2Len);
						ippsZero_32fc(cdst, cdstLen);
						ippsZero_32f(cdstMag, cdstLen);
						ippsZero_32f(cdstPha, cdstLen);
						ippsZero_64f(csrc2Norm, cdstLen);
						ippsZero_32f(csrcNorm, cdstLen);
						ippsZero_32fc(csrcNormCplx, cdstLen);
						//解析点の算出
						apdecimal = (4 * frq_s) / cc[j] * (dep + sqrt(pow(dep, 2) + pow(xi[fine_ele[l]], 2) - 2 * dep * xi[fine_ele[l]] * sin(theta[i])));
						ap = static_cast<int>(apdecimal);
						if (apdecimal - static_cast<int>(apdecimal) >= 0.5)
							++ap;
						Lag = ap - csrc1Len / 2 + 1;
						if (l == 30)
							cout << "ch:" << fine_ele[l] << "ss: " << cc[j] << "lag: " << dep << "\n";
						for (int m = 0; m < 4 * sample; ++m){ //csrc2を設定
							csrc2[m].re = elere[i][fine_ele[l]][m];
							csrc2[m].im = eleim[i][fine_ele[l]][m];
						}
						//相関
						status = ippsCrossCorrNorm_32fc(csrc1, csrc1Len, csrc2, csrc2Len, cdst, cdstLen, Lag + lowLag, NormA, pbuffer);

						for (int m = 0; m < cdstLen; ++m){ //src2のノルム
							csrc2Normtmp = 0.0;
							ippsNorm_L2_32fc64f(csrc2 + Lag + lowLag + m, csrc1Len, &csrc2Normtmp);
							csrc2Norm[m] = csrc2Normtmp;
						}
						ippsMulC_64f_I(csrc1Norm, csrc2Norm, cdstLen); //src1のノルム*src2のノルム
						ippsConvert_64f32f(csrc2Norm, csrcNorm, cdstLen); //double->float
						ippsRealToCplx_32f(csrcNorm, cdstLenzeros, csrcNormCplx, cdstLen); //実数列->複素列
						//ippsDiv_32fc_I(csrcNormCplx, cdst, cdstLen); //正規化
						ippsMagnitude_32fc(cdst, cdstMag, cdstLen); //相関振幅
						ippsPhase_32fc(cdst, cdstPha, cdstLen); //相関位相
						ippsMaxIndx_32f(cdstMag, cdstLen, &cmax, &cmaxidx);
						for (int m = 0; m < cdstLen; ++m){
							corrre[fine_ele[l]][j][m] = cdst[m].re;
							corrim[fine_ele[l]][j][m] = cdst[m].im;
							corrmag[fine_ele[l]][j][m] = cdstMag[m];
							corrpha[fine_ele[l]][j][m] = cdstPha[m];
						}


						//結果を出力
						//cout << "Max:" << cmax << ", Diff: " << cmaxidx + lowLag << "\n";

						//cfout << cc[j] << " " << fine_ele[l] << " " << cmaxidx + lowLag << "\n";
					}
				}
			}
			//cfout << "\n";
		}

	}
	//cfout.close();

	/*vector<float> dstx(cdstLen);
	for (int i = 0; i < cdstLen; ++i)
		dstx[i] = (i + lowLag) / (4 * frq_s);

	plt::named_plot("1490 m/s", dstx, corrmag[78][0], "-");
	plt::named_plot("1540 m/s", dstx, corrmag[78][5], "--");
	plt::named_plot("1590 m/s", dstx, corrmag[78][10], "-.");
	plt::legend();
	plt::grid(true);
	plt::show();

	cout << "finish!\n";
	cout << "!?\n";*/


	//frame = 2;
	//initialize focused RF array(vector)
	//vector<vector<vector<vector<float>>>> elere(frame,
	//	vector<vector<vector<float>>>(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0))));
	//vector<vector<vector<vector<float>>>> eleim(frame,
	//	vector<vector<vector<float>>>(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0))));

	////spec and buffer setting for FFT
	//Ipp8u *specbuff, *initbuff, *workbuff;
	//Ipp8u *specbufi, *initbufi, *workbufi;
	//int size_specf, size_initf, size_workf;
	//int size_speci, size_initi, size_worki;
	//IppsFFTSpec_C_32fc *specf = 0;
	//IppsFFTSpec_C_32fc *speci = 0;
	//Ipp32fc *ipsrc = ippsMalloc_32fc((int)sample);
	//Ipp32fc *ipdst = ippsMalloc_32fc((int)sample);
	//Ipp32fc *ipsrc2 = ippsMalloc_32fc((int)(4 * sample));
	//Ipp32fc *ipdst2 = ippsMalloc_32fc((int)(4 * sample));
	//const int fftorder = (int)(log((double)sample) / log(2.0));
	//const int ifftorder = (int)(log((double)(4 * sample)) / log(2.0));
	//ippsFFTGetSize_C_32fc(fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_specf, &size_initf, &size_workf);
	//ippsFFTGetSize_C_32fc(ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_speci, &size_initi, &size_worki);
	//specbuff = ippsMalloc_8u(size_specf);
	//specbufi = ippsMalloc_8u(size_speci);
	//initbuff = ippsMalloc_8u(size_initf);
	//initbufi = ippsMalloc_8u(size_initi);
	//workbuff = ippsMalloc_8u(size_workf);
	//workbufi = ippsMalloc_8u(size_worki);
	//ippsFFTInit_C_32fc(&specf, fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuff, initbuff);
	//ippsFFTInit_C_32fc(&speci, ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbufi, initbufi);
	//
	//for (int i = 0; i < frame; ++i){
	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < ch; ++k){
	//			ippsZero_32fc(ipsrc, sample);
	//			ippsZero_32fc(ipdst, sample);
	//			ippsZero_32fc(ipsrc2, 4 * sample);
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//			//set
	//			for (int l = 0; l < sample - 1; ++l)
	//				ipsrc[l].re = raw.RF[i][j][k][l];
	//				//ipsrc[l].re = raw.RF[60][j][k][l] - raw.RF[59][j][k][l];
	//
	//			//do FFT
	//			ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
	//			ippsZero_8u(workbuff, size_workf);
	//			//double positive part and delete negative part
	//			for (int l = 0; l < sample / 2; ++l){
	//				ipdst[l].re = ipdst[l].re * 2 / sample;
	//				ipdst[l].im = ipdst[l].im * 2 / sample;
	//				ipdst[l + sample / 2].re = 0.0;
	//				ipdst[l + sample / 2].im = 0.0;
	//			}
	//			for (int l = 0; l < 34; ++l){
	//				ipdst[l].re = 0.0;
	//				ipdst[l].im = 0.0;
	//			}
	//			
	//			
	//			//ipdst[0].re /= 2;
	//			//ipdst[0].im /= 2;

	//			for (int l = 0; l < sample; ++l){
	//				ipsrc2[l].re = ipdst[l].re;
	//				ipsrc2[l].im = ipdst[l].im;
	//			}


	//			//do IFFT
	//			ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
	//			ippsZero_8u(workbufi, size_worki);
	//			
	//			//save
	//			for (int l = 0; l < 4 * sample; ++l){
	//				elere[i][j][k][l] = 4 * sample * ipdst2[l].re;
	//				eleim[i][j][k][l] = 4 * sample * ipdst2[l].im;
	//			}
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//		}
	//	}
	//}

	////free IPP array
	//ippsFree(ipsrc);
	//ippsFree(ipdst);
	//ippsFree(ipdst2);
	////free RF
	////vector<vector<vector<vector<short>>>>().swap(RF);
	//raw.freeRF();

	///* interpolation */
	//cout << "interpolating...\n";

	//vector<vector<vector<float>>> RFre(frame, vector<vector<float>>(line, vector<float>(sample, 0)));
	//vector<vector<vector<float>>> RFim(frame, vector<vector<float>>(line, vector<float>(sample, 0)));

	////calculate delay
	//const float c0 = 1540.0;
	//int point, add;
	//float eledep; //in-bound(um)
	//float decimal; //decimal part in sampling point of round trip distance
	//vector<float> xi(ch, 0); // x-coordinate of each element
	//for (int i = 0; i < ch; ++i)
	//	xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	//vector<float> theta(line, 0);
	//for (int i = 0; i < line; ++i) //beam angle
	//	theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	//vector<float> cendep(sample, 0);
	//for (int i = 0; i < sample; ++i)
	//	cendep[i] = i * (c0 / (2 * frq_s)); //out-bound(um)

	////addition
	//for (int i = 0; i < frame; ++i){
	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < sample; ++k){
	//			add = 0;
	//			for (int l = 0; l < ch; ++l){
	//				eledep = sqrt(pow(xi[l], 2) + pow(cendep[k], 2) - 2 * xi[l] * cendep[k] * sin(theta[j]));
	//				point = static_cast<int>(((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)));
	//				decimal = ((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)) - point;
	//				if (point < 4 * sample - 1){
	//					RFre[i][j][k] += (elere[i][j][l][point] + (elere[i][j][l][point + 1] - elere[i][j][l][point]) * decimal);
	//					RFim[i][j][k] += (eleim[i][j][l][point] + (eleim[i][j][l][point + 1] - eleim[i][j][l][point]) * decimal);
	//					//RFre[i][j][k] += elere[i][j][l][point];
	//					//RFim[i][j][k] += eleim[i][j][l][point];
	//					++add;
	//				}
	//			}
	//			if (add != 0){
	//				RFre[i][j][k] /= add;
	//				RFim[i][j][k] /= add;
	//			}
	//			else{
	//				RFre[i][j][k] = 0.0;
	//				RFim[i][j][k] = 0.0;
	//			}
	//		}
	//	}
	//}

	////free eledata
	//vector<vector<vector<vector<float>>>>().swap(elere);
	//vector<vector<vector<vector<float>>>>().swap(eleim);

	//vector<vector<vector<float>>> env(frame, vector<vector<float>>(line, vector<float>(sample, 0)));
	//for (int i = 0; i < frame; ++i)
	//		for (int k = 0; k < sample; ++k)
	//			env[i][j][k] = sqrt(pow(RFre[i][j][k], 2) + pow(RFim[i][j][k], 2));

	//int bn;
	////cin >> bn;
	//ofstream fout2("./new/8.dat", ios_base::out);
	///*for (int i = 0; i < sample; ++i){
	//	fout2 << i << " " << raw.RF0[60][8][i] << "\n";
	//}
	//fout2.close();*/
	//
	//for (int i = 0; i < ch; ++i){
	//	/*ost << "./" << dir << "/" << i << ".dat";*/
	//	ostringstream ost;
	//	ost << "./new/" << i << ".dat";
	//	fout2.open(ost.str(), ios_base::out);
	//	ost.clear();
	//	ost.str("");
	//	for (int k = 0; k < sample; ++k){
	//		fout2 << k << " " << raw.RF0[0][i][k] << "\n";
	//	}
	//	fout2.close();
	//}

	//cout << "ok\n";
	//string savedir = "0622";
	//raw.plotRF0(savedir + "/" + RFname);
	//raw.loadRF0(2);
	//vector<vector<double>> env = raw.calcenv(1, max_angle, frq_s);
	//BSector22(raw.calcenv(1, max_angle, frq_s), max_angle, frq_s, 2);

	//est_ss est(5, 3, 3); // (grid_w, grid_h, grid/beam)

	//est.set_parameter(max_angle, frq_s, frq_t);
	//est.loadRF(raw.RF0);
	//est.sim_RFset();

	//est.calc_delay(40.0); //40.0mm
	//est.sim_delay_setspeed();
	//est.calc_path();
	//est.sim_delay_calcdelay();
	//est.del_brokenelement();
	//est.write_mat();

	//est.SVD();

	/*bool checkbr;
	int avframe = 12;
	int avline = 11;*/
	//vector<vector<short>> rf(ch, vector<short>(sample, 0));
	//int sel_beamnum = 81;
	//rf = raw.RF0[sel_beamnum];
	//raw.freeRF0();
	//
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < sample - 1; ++j){
	//		rf[i][j] = raw.RF[avframe + 1][avline][i][j] - raw.RF[avframe][avline][i][j];
	//	}
	//}

	//Ipp32f *hilx = ippsMalloc_32f(sample);
	//Ipp32fc *hily = ippsMalloc_32fc(sample);
	//IppsHilbertSpec_32f32fc *hilspec;
	//IppStatus hilst;
	//hilst = ippsHilbertInitAlloc_32f32fc(&hilspec, sample, ippAlgHintFast);
	//vector<vector<float>> anare(ch, vector<float>(sample, 0));
	//vector<vector<float>> anaim(ch, vector<float>(sample, 0));
	//for (int i = 0; i < ch; ++i){
	//	ippsZero_32f(hilx, sample);
	//	ippsZero_32fc(hily, sample);
	//	for (int j = 0; j < sample - 1; ++j){
	//		hilx[j] = rf[i][j];
	//	}
	//	hilst = ippsHilbert_32f32fc(hilx, hily, hilspec);
	//	for (int j = 0; j < sample; ++j){
	//		anare[i][j] = hily[j].re;
	//		anaim[i][j] = hily[j].im;
	//	}
	//}
	//vector<vector<short>>().swap(rf);
	//ippsHilbertFree_32f32fc(hilspec);
	//ippsFree(hilx);
	//ippsFree(hily);

	//
	//Ipp64f std1, std2;
	//Ipp32fc mean1, mean2;
	//IppStatus corrst;
	//IppEnum NormA = (IppEnum)(ippAlgAuto | ippsNormNone);
	//int bufsize = 0;
	//int lowlag = 1500;
	//Ipp8u *pbuffer;
	///*const int src1Len = 256;
	//const int src2Len = sample;
	//const int dstLen = 128;*/
	//const int src1Len = 500;
	//const int src2Len = 500;
	//const int dstLen = src1Len * 2;
	//Ipp32fc *pSrc1 = ippsMalloc_32fc(src1Len);
	//Ipp32fc *pSrc2 = ippsMalloc_32fc(src2Len);
	//Ipp32fc *pDst = ippsMalloc_32fc(dstLen);
	//Ipp32fc *ptmp = ippsMalloc_32fc(src1Len);
	//Ipp32f *ptmp2 = ippsMalloc_32f(src2Len);
	//Ipp32f *ptmp3 = ippsMalloc_32f(dstLen);
	//corrst = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, lowlag, ipp32fc, NormA, &bufsize);
	//pbuffer = ippsMalloc_8u(bufsize);
	//vector<float> norm(dstLen, 0);
	//vector<vector<float>> corrre(ch, vector<float>(dstLen, 0));
	//vector<vector<float>> corrim(ch, vector<float>(dstLen, 0));
	//vector<vector<float>> corramp(ch, vector<float>(dstLen, 0));
	//float normtmp;

	//for (int i = 0; i < src1Len; ++i){
	//	pSrc1[i].re = anare[ch - 1][i + lowlag];
	//	pSrc1[i].im = anaim[ch - 1][i + lowlag];
	//}
	////ippsMean_32fc(pSrc1, src1Len, &mean1, ippAlgHintFast);
	////ippsSubCRev_32fc(pSrc1, mean1, ptmp, src1Len);
	////ippsNorm_L2_32fc64f(ptmp, src1Len, &std1);
	//ippsNorm_L2_32fc64f(pSrc1, src1Len, &std1);

	//for (int i = 0; i < ch; ++i){
	//	ippsZero_32fc(pSrc2, sample);
	//	ippsZero_32fc(pDst, dstLen);
	//	ippsZero_32fc(ptmp, src1Len);
	//	ippsZero_32f(ptmp2, src2Len);
	//	for (int j = 0; j < sample; ++j){
	//		pSrc2[j].re = anare[i][j];
	//		pSrc2[j].im = anaim[i][j];
	//		if (j < src1Len){
	//			ptmp[j].re = anare[i][j];
	//			ptmp[j].im = anaim[i][j];
	//		}
	//	}
	//	corrst = ippsCrossCorrNorm_32fc(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowlag, NormA, pbuffer);
	//	ippsNorm_L2_32fc64f(ptmp, src1Len, &std2);
	//	ippsMagnitude_32fc(pSrc2, ptmp2, src2Len);
	//	ippsMagnitude_32fc(pDst, ptmp3, dstLen);
	//	
	//	for (int j = 0; j < dstLen; ++j){
	//		ippsNorm_L2_32fc64f(pSrc2 + j, src1Len, &std2);
	//		norm[j] = std2;
	//	}
	//	for (int j = 0; j < dstLen; ++j){
	//		corrre[i][j] = pDst[j].re / std1 / norm[j];
	//		corrim[i][j] = pDst[j].im / std1 / norm[j];
	//		corramp[i][j] = ptmp3[j] / std1 / norm[j];
	//	}
	//}
	//ofstream foutcorr("corr.dat", ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < dstLen; ++j){
	//		foutcorr << j << " " << corramp[i][j] - i * 2 << "\n";
	//	}
	//	foutcorr << "\n";
	//}
	//foutcorr.close();
	//vector<int> mv(ch, 0);
	//vector<float>::iterator itr;
	//for (int i = 0; i < ch; ++i){
	//	itr = max_element(corramp[i].begin(), corramp[i].end());
	//	mv[i] = distance(corramp[i].begin(), itr);
	//}
	//ofstream fout("correlation.dat", ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	fout << i << " " << mv[i] << "\n";
	//}

	///*frequency analytic*/
	//float han[60];
	//for (int i = 0; i < 60; ++i)
	//	han[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / 60);

	//Ipp8u *sbuf, *ibuf, *wbuf;
	//int s_s, s_i, s_w;
	//IppsFFTSpec_C_32fc *sp = 0;
	//Ipp32fc *src = ippsMalloc_32fc(1024);
	//ippsFFTGetSize_C_32fc(10, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &s_s, &s_i, &s_w);
	//sbuf = ippsMalloc_8u(s_s);
	//ibuf = ippsMalloc_8u(s_i);
	//wbuf = ippsMalloc_8u(s_w);
	//ippsFFTInit_C_32fc(&sp, 10, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, sbuf, ibuf);
	//for (int i = 0; i < 60; ++i)
	//	src[i].re = RF[60][17][90][i + 1270] * han[i];
	//for (int i = 0; i < 1024 - 60; ++i)
	//	src[i + 60].re = 0.0;
	//ippsFFTFwd_CToC_32fc(src, src, sp, wbuf);
	//float freqi, tmpfc;
	//ofstream spect("./spectrum.dat", ios_base::out);
	//for (int i = 0; i < 512; ++i){
	//	freqi = static_cast<float>(i) * 15 / 512;
	//	tmpfc = 10.0 * log10(pow(src[i].re, 2) + pow(src[i].im, 2));
	//	spect << freqi << " " << tmpfc << "\n";
	//}
	//spect.close();

	///*calculate differences between frames before and after*/
	//short tmps, th;
	//th = 500; // threshold
	//for (int i = 0; i < 10; ++i){
	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < ch; ++k){
	//			for (int l = 0; l < sample - 1; ++l){
	//				tmps = abs(RF[i + 55][j][k][l] - RF[i + 54][j][k][l]);
	//				if (tmps > th && i==5)
	//					cout << "over " << th << ":(" << i + 55 << " " << j << " " << k << " " << l << ")\n";
	//			}
	//		}
	//	}
	//}

	//raw.loadRF();
	//raw.loadRF0(60);
	/* do FFT and IFFT */
	/*string fname2 = "./raw.dat";
	ofstream fout2(fname2, ios_base::out);
	for (int i = 0; i < ch; ++i){
	for (int j = 0; j < sample - 1; ++j)
	fout2 << j << " " << raw.RF0[30][i][j] - i * 4096 << "\n";
	fout2 << "\n";
	}
	fout2.close();*/

	/*ostringstream fst;
	ofstream foutt(fst.str(), ios_base::out);
	for (int i = 0; i < 7; ++i){
	for (int j = 0; j < line; ++j){
	fst << "./diff0104/fr" << i + 10 << "/l" << j << ".dat";
	foutt.open(fst.str(), ios_base::out);
	fst.clear();
	fst.str("");
	for (int k = 0; k < ch; ++k){
	for (int l = 0; l < sample - 1; ++l){
	foutt << l << " " << raw.RF[i + 10][j][k][l] - raw.RF[i + 9][j][k][l] - k * 4096 << "\n";
	}
	foutt << "\n";
	}
	foutt.close();
	}
	}*/
	//cross correlation
	/*vector<int> b_ele = { 8, 9, 10, 11, 12, 13, 14, 19, 42, 47, 56, 79, 80 };
	IppStatus status;
	int lowlag = -500;
	IppEnum NormA = (IppEnum)(ippAlgAuto | ippsNormA);
	int bufSize = 0;
	Ipp8u *pBuffer;
	const int srcLen = 500;
	const int dstLen = 2 * srcLen;
	Ipp32f *pSrc1 = ippsMalloc_32f(srcLen);
	Ipp32f *pSrc2 = ippsMalloc_32f(srcLen);
	Ipp32f *pDst = ippsMalloc_32f(dstLen);
	status = ippsCrossCorrNormGetBufferSize(srcLen, srcLen, dstLen, lowlag, ipp32f, NormA, &bufSize);

	pBuffer = ippsMalloc_8u(bufSize);*/


	//vector<vector<vector<short>>> rfp(line, vector<vector<short>>(ch, vector<short>(sample - 1, 0)));
	//raw.loadRF0(0);
	//rfp = raw.RF0;

	////spec and buffer setting for FFT
	//Ipp8u *specbuff, *initbuff, *workbuff;
	//Ipp8u *specbufi, *initbufi, *workbufi;
	//int size_specf, size_initf, size_workf;
	//int size_speci, size_initi, size_worki;
	//IppsFFTSpec_C_32fc *specf = 0;
	//IppsFFTSpec_C_32fc *speci = 0;
	//Ipp32fc *ipsrc = ippsMalloc_32fc((int)sample);
	//Ipp32fc *ipdst = ippsMalloc_32fc((int)sample);
	//Ipp32fc *ipsrc2 = ippsMalloc_32fc((int)(4 * sample));
	//Ipp32fc *ipdst2 = ippsMalloc_32fc((int)(4 * sample));
	//const int fftorder = (int)(log((double)sample) / log(2.0));
	//const int ifftorder = (int)(log((double)(4 * sample)) / log(2.0));
	//ippsFFTGetSize_C_32fc(fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_specf, &size_initf, &size_workf);
	//ippsFFTGetSize_C_32fc(ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_speci, &size_initi, &size_worki);
	//specbuff = ippsMalloc_8u(size_specf);
	//specbufi = ippsMalloc_8u(size_speci);
	//initbuff = ippsMalloc_8u(size_initf);
	//initbufi = ippsMalloc_8u(size_initi);
	//workbuff = ippsMalloc_8u(size_workf);
	//workbufi = ippsMalloc_8u(size_worki);
	//ippsFFTInit_C_32fc(&specf, fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuff, initbuff);
	//ippsFFTInit_C_32fc(&speci, ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbufi, initbufi);

	//vector<vector<vector<float>>> elere(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0)));
	//vector<vector<vector<float>>> eleim(line, vector<vector<float>>(ch, vector<float>(4 * sample, 0)));

	////calculate delay
	//const float c0 = 1540.0;
	//int point, add;
	//float eledep; //in-bound(um)
	//float decimal; //decimal part in sampling point of round trip distance
	//vector<float> xi(ch, 0); // x-coordinate of each element
	//for (int i = 0; i < ch; ++i)
	//	xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	//vector<float> theta(line, 0);
	//for (int i = 0; i < line; ++i) //beam angle
	//	theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	//vector<float> cendep(sample, 0);
	//for (int i = 0; i < sample; ++i)
	//	cendep[i] = i * (c0 / (2 * frq_s)); //out-bound(um)

	//vector<vector<float>> RFre(line, vector<float>(sample, 0));
	//vector<vector<float>> RFim(line, vector<float>(sample, 0));

	//vector<vector<float>> env(line, vector<float>(sample, 0));

	//float maxa;
	//for (int i = 1; i < frame; ++i){
	//	raw.freeRF0();
	//	raw.loadRF0(i);

	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < ch; ++k){
	//			ippsZero_32fc(ipsrc, sample);
	//			ippsZero_32fc(ipdst, sample);
	//			ippsZero_32fc(ipsrc2, 4 * sample);
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//			//set
	//			for (int l = 0; l < sample - 1; ++l){
	//				ipsrc[l].re = rfp[j][k][l];
	//				ipsrc[l].im = 0.0;
	//			}

	//			//do FFT
	//			ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
	//			ippsZero_8u(workbuff, size_workf);
	//			//double positive part and delete negative part
	//			for (int l = 0; l < sample / 2; ++l){
	//				ipdst[l].re = ipdst[l].re * 2 / sample;
	//				ipdst[l].im = ipdst[l].im * 2 / sample;
	//				ipdst[l + sample / 2].re = 0.0;
	//				ipdst[l + sample / 2].im = 0.0;
	//			}
	//			for (int l = 0; l < 34; ++l){
	//				ipdst[l].re = 0.0;
	//				ipdst[l].im = 0.0;
	//			}

	//			for (int l = 0; l < sample; ++l){
	//				ipsrc2[l].re = ipdst[l].re;
	//				ipsrc2[l].im = ipdst[l].im;
	//			}


	//			//do IFFT
	//			ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
	//			ippsZero_8u(workbufi, size_worki);

	//			//save
	//			for (int l = 0; l < 4 * sample; ++l){
	//				elere[j][k][l] = 4 * sample * ipdst2[l].re;
	//				eleim[j][k][l] = 4 * sample * ipdst2[l].im;
	//			}
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//		}
	//	}


	//	//addition

	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < sample; ++k){
	//			add = 0;
	//			for (int l = 0; l < ch; ++l){
	//				eledep = sqrt(pow(xi[l], 2) + pow(cendep[k], 2) - 2 * xi[l] * cendep[k] * sin(theta[j]));
	//				point = static_cast<int>(((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)));
	//				decimal = ((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)) - point;
	//				if (point < 4 * sample - 1){
	//					RFre[j][k] += (elere[j][l][point] + (elere[j][l][point + 1] - elere[j][l][point]) * decimal);
	//					RFim[j][k] += (eleim[j][l][point] + (eleim[j][l][point + 1] - eleim[j][l][point]) * decimal);
	//					//RFre[i][j][k] += elere[i][j][l][point];
	//					//RFim[i][j][k] += eleim[i][j][l][point];
	//					++add;
	//				}
	//			}
	//			if (add != 0){
	//				RFre[j][k] /= add;
	//				RFim[j][k] /= add;
	//			}
	//			else{
	//				RFre[j][k] = 0.0;
	//				RFim[j][k] = 0.0;
	//			}
	//		}
	//	}

	//	for (int j = 0; j < line; ++j)
	//		for (int k = 0; k < sample; ++k)
	//			env[j][k] = sqrt(pow(RFre[j][k], 2) + pow(RFim[j][k], 2));

	//	maxa = BSector22(env, max_angle, frq_s, i);

	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < ch; ++k){
	//			ippsZero_32fc(ipsrc, sample);
	//			ippsZero_32fc(ipdst, sample);
	//			ippsZero_32fc(ipsrc2, 4 * sample);
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//			//set
	//			for (int l = 0; l < sample - 1; ++l){
	//				ipsrc[l].re = raw.RF0[j][k][l] - rfp[j][k][l];
	//				ipsrc[l].im = 0.0;
	//			}

	//			//do FFT
	//			ippsFFTFwd_CToC_32fc(ipsrc, ipdst, specf, workbuff);
	//			ippsZero_8u(workbuff, size_workf);
	//			//double positive part and delete negative part
	//			for (int l = 0; l < sample / 2; ++l){
	//				ipdst[l].re = ipdst[l].re * 2 / sample;
	//				ipdst[l].im = ipdst[l].im * 2 / sample;
	//				ipdst[l + sample / 2].re = 0.0;
	//				ipdst[l + sample / 2].im = 0.0;
	//			}
	//			for (int l = 0; l < 34; ++l){
	//				ipdst[l].re = 0.0;
	//				ipdst[l].im = 0.0;
	//			}

	//			for (int l = 0; l < sample; ++l){
	//				ipsrc2[l].re = ipdst[l].re;
	//				ipsrc2[l].im = ipdst[l].im;
	//			}


	//			//do IFFT
	//			ippsFFTInv_CToC_32fc(ipsrc2, ipdst2, speci, workbufi);
	//			ippsZero_8u(workbufi, size_worki);

	//			//save
	//			for (int l = 0; l < 4 * sample; ++l){
	//				elere[j][k][l] = 4 * sample * ipdst2[l].re;
	//				eleim[j][k][l] = 4 * sample * ipdst2[l].im;
	//			}
	//			ippsZero_32fc(ipdst2, 4 * sample);
	//		}
	//	}


	//	//addition

	//	for (int j = 0; j < line; ++j){
	//		for (int k = 0; k < sample; ++k){
	//			add = 0;
	//			for (int l = 0; l < ch; ++l){
	//				eledep = sqrt(pow(xi[l], 2) + pow(cendep[k], 2) - 2 * xi[l] * cendep[k] * sin(theta[j]));
	//				point = static_cast<int>(((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)));
	//				decimal = ((cendep[k] + eledep) / 2) / (c0 / (8 * frq_s)) - point;
	//				if (point < 4 * sample - 1){
	//					RFre[j][k] += (elere[j][l][point] + (elere[j][l][point + 1] - elere[j][l][point]) * decimal);
	//					RFim[j][k] += (eleim[j][l][point] + (eleim[j][l][point + 1] - eleim[j][l][point]) * decimal);
	//					//RFre[i][j][k] += elere[i][j][l][point];
	//					//RFim[i][j][k] += eleim[i][j][l][point];
	//					++add;
	//				}
	//			}
	//			if (add != 0){
	//				RFre[j][k] /= add;
	//				RFim[j][k] /= add;
	//			}
	//			else{
	//				RFre[j][k] = 0.0;
	//				RFim[j][k] = 0.0;
	//			}
	//		}
	//	}

	//	for (int j = 0; j < line; ++j)
	//		for (int k = 0; k < sample; ++k)
	//			env[j][k] = sqrt(pow(RFre[j][k], 2) + pow(RFim[j][k], 2));

	//	BSector2(env, max_angle, frq_s, i, maxa);

	//	rfp = raw.RF0;
	//}






	//ofstream fou("elediftime.dat", ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < sample - 1; ++j){
	//		fou << double(j) / double(frq_s) << " " << raw.RF[avframe + 1][avline][i][j] - raw.RF[avframe][avline][i][j] - i * 4096 << "\n";
	//	}
	//	fou << "\n";
	//}
	//fou.close();

	//ofstream fou2("elediftime95.dat", ios_base::out);
	//for (int i = 95; i < ch; ++i){
	//	for (int j = 0; j < sample - 1; ++j){
	//		fou2 << double(j) / double(frq_s) << " " << raw.RF[avframe + 1][avline][i][j] - raw.RF[avframe][avline][i][j] << "\n";
	//	}
	//	fou2 << "\n";
	//}
	//fou2.close();

	//for (int i = 0; i < srcLen; ++i){
	//	pSrc1[i] = raw.RF[avframe + 1][avline][ch - 1][1500 + i] - raw.RF[avframe][avline][ch - 1][1500 + i];
	//	//pSrc1[i] = 1.0;
	//}
	//
	//Ipp32f stdev1, stdev2;
	//ippsStdDev_32f(pSrc1, srcLen, &stdev1, ippAlgHintFast);

	//vector<vector<float>> cc(ch, vector<float>(dstLen, 0));
	//ofstream foutcorr("corr.dat", ios_base::out);



	//for (int i = 0; i < ch; ++i){
	//	checkbr = any_of(b_ele.begin(), b_ele.end(), [i](int x){return x == i; });
	//	ippsZero_32f(pSrc2, srcLen);
	//	ippsZero_32f(pDst, dstLen);
	//	stdev2 = 0.0;
	//	
	//	for (int j = 0; j < srcLen; ++j){
	//		pSrc2[j] = raw.RF[avframe + 1][avline][i][1500 + j] - raw.RF[avframe][avline][i][1500 + j];
	//		//pSrc2[j] = 1.0;
	//	}
	//	ippsStdDev_32f(pSrc2, srcLen, &stdev2, ippAlgHintAccurate);
	//	status = ippsCrossCorrNorm_32f(pSrc1, srcLen, pSrc2, srcLen, pDst, dstLen, lowlag, NormA, pBuffer);

	//	for (int j = 0; j < dstLen; ++j){
	//		if (stdev1 * stdev2 >= 1e-9)
	//			cc[i][j] = pDst[j] / (stdev1 * stdev2);
	//		else cc[i][j] = 0.0;

	//		//plot
	//		if (!checkbr && i >= 14)
	//			foutcorr << double(j - 500) / double(frq_s) << " " << cc[i][j] - i * 2 << "\n";
	//	}
	//	foutcorr << "\n";
	//}
	//foutcorr.close();


	//ippsFree(pBuffer);

	//vector<size_t> maxcc(ch, 0);
	//for (int i = 0; i < ch; ++i){
	//	vector<float>::iterator itrr = max_element(cc[i].begin(), cc[i].end());
	//	maxcc[i] = distance(cc[i].begin(), itrr);
	//}
	//
	////ippsAbs_32f_I(pSrc1, srcLen);
	//Ipp32f pMax;
	//int pIndx;
	//ippsMaxIndx_32f(pSrc1, srcLen, &pMax, &pIndx);
	//float maxoff = (pIndx + 1500) / frq_s;

	//ofstream maxc("maxc.dat", ios_base::out); //delay

	//
	//for (int i = 0; i < ch; ++i){
	//	checkbr = any_of(b_ele.begin(), b_ele.end(), [i](int x){return x == i; });
	//	if (!checkbr && i >= 14)
	//		maxc << i << " " << (static_cast<int>(maxcc[i]) - 500) / frq_s << "\n";
	//}
	//maxc.close();

	//ofstream fout3("elemplus.dat", ios_base::out); //delay (raw signal)
	//for (int i = 0; i < ch; ++i){
	//	checkbr = any_of(b_ele.begin(), b_ele.end(), [i](int x){return x == i; });
	//	if (!checkbr && i >= 14)
	//		fout3 << (static_cast<int>(maxcc[i]) - 500) / frq_s << " " << -2 * i + 1 << "\n";
	//		//fout3 << i << " " << (static_cast<int>(maxcc[i]) - 500) / frq_s + maxoff << "\n";
	//}
	//fout3.close();

	//string fname = "./ele.dat";
	//ofstream fout(fname, ios_base::out);
	//for (int i = 0; i < ch; ++i){
	//	for (int j = 0; j < dstLen; ++j){
	//		fout << j << " " << cc[i][j] - 2.0 * i << "\n";
	//	}
	//	fout << "\n";
	//}
	//fout.close();



	//Bsector(env[0], max_angle);
	//BSector2(env[0], max_angle, frq_s);
	//Bsector3(env[0], max_angle);
	//cairo(env[0], max_angle);
	
	return 0;
}

IppStatus CrossCorrNormExample(void) {
	IppStatus status;
	const int src1Len = 5, src2Len = 7, dstLen = 16;
	int lowLag = -5;
	Ipp32f pSrc1[src1Len] = { 1.f, 1.f, 3.f, 1.f, 1.f }, pSrc2[src2Len] = { 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f }, pDst[dstLen];
	IppEnum funCfgNormNo = (IppEnum)(ippAlgAuto | ippsNormNone);
	IppEnum funCfgNormA = (IppEnum)(ippAlgAuto | ippsNormA);
	IppEnum funCfgNormB = (IppEnum)(ippAlgAuto | ippsNormB);
	int bufSizeNo = 0, bufSizeA = 0, bufSizeB = 0, bufSizeMax = 0;
	Ipp8u *pBuffer;

	status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, -5, ipp32f, funCfgNormNo, &bufSizeNo);
	if (status != ippStsNoErr) return status;
	status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, -5, ipp32f, funCfgNormA, &bufSizeA);
	if (status != ippStsNoErr) return status;
	status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, dstLen, -5, ipp32f, funCfgNormB, &bufSizeB);
	if (status != ippStsNoErr) return status;

	bufSizeMax = IPP_MAX(bufSizeNo, IPP_MAX(bufSizeA, bufSizeB));// get max buffer size

	pBuffer = ippsMalloc_8u(bufSizeMax);

	status = ippsCrossCorrNorm_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, funCfgNormNo, pBuffer);
	//printf_32("pDst_NormNone", pDst, dstLen);
	cout << "NormNone -> " << fixed << setprecision(3);
	for (int i = 0; i < dstLen; ++i)
		cout << pDst[i] << " ";
	cout << "\n";

	status = ippsCrossCorrNorm_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, funCfgNormA, pBuffer);
	//printf_32f("pDst_NormA", pDst, dstLen);
	cout << "NormA    -> ";
	for (int i = 0; i < dstLen; ++i)
		cout << pDst[i] << " ";
	cout << "\n";

	status = ippsCrossCorrNorm_32f(pSrc1, src1Len, pSrc2, src2Len, pDst, dstLen, lowLag, funCfgNormB, pBuffer);
	//printf_32f("pDst_NormB", pDst, dstLen);
	cout << "NormB    -> ";
	for (int i = 0; i < dstLen; ++i)
		cout << pDst[i] << " ";
	cout << "\n";

	ippsFree(pBuffer);
	return status;
}

struct misra1a_functor
{
	misra1a_functor(int inputs, int values, double *x, double *y)
		: inputs_(inputs), values_(values), x(x), y(y) {}

	double *x;
	double *y;

	double p = 0.2e3;
	double l = 9.5e3;
	double sn = 0.15646446504;

	// 目的関数
	int operator()(const VectorXd& b, VectorXd& fvec) const
	{
		for (int i = 0; i < values_; ++i) {
			fvec[i] = (b[0] + sqrt(pow(p * x[i], 2) + 2.0 * (b[0] * sn - l)*p*x[i] + pow(b[0], 2) + pow(l, 2) - 2.0 * b[0] * l*sn)) / b[1] - y[i];
		}
		return 0;
	}
	// 微分,ヤコビアン
	int df(const VectorXd& b, MatrixXd& fjac)
	{
		for (int i = 0; i < values_; ++i) {
			fjac(i, 0) = 1 / b[1] + (b[0] + sn*(p*x[i] - l)) / sqrt(pow(p * x[i], 2) + 2.0 * (b[0] * sn - l)*p*x[i] + pow(b[0], 2) + pow(l, 2) - 2.0 * b[0] * l*sn);
			fjac(i, 1) = -(b[0] + sqrt(pow(p * x[i], 2) + 2.0 * (b[0] * sn - l)*p*x[i] + pow(b[0], 2) + pow(l, 2) - 2.0 * b[0] * l*sn)) / b[1] / b[1];
		}
		return 0;
	}

	const int inputs_;
	const int values_;
	int inputs() const { return inputs_; }
	int values() const { return values_; }
};