// fileopen.cpp
/*This program defines the class for opening and storing alpha-10 element data.
	declaration part is in "header.h".*/

#include "stdafx.h"
#include "header.h"

using namespace std;

file::file(string filename)
{
	cout << "open file=" << filename << "\n";
	fin.open(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
		
	}
}

file::file()
{

}

file::~file()
{
	cout << "delete object.\n";
}

void file::open(string filename)
{
	fin.open(filename, ios_base::in | ios_base::binary);
	if (!fin){
		cout << "couldn't load file.\n";
	}
}

void file::start()
{
	fin.clear();
	fin.seekg(0, ios_base::beg);
}

void file::warp(int pos)
{
	fin.seekg(pos, ios_base::beg);
}

void file::go(int pos)
{
	fin.seekg(pos, ios_base::cur);
}


/*a10*/
a10::a10(string filename) : file(filename)
{
	loadheader();
}

a10::a10()
{

}

a10::~a10()
{
}

void a10::loadheader()
{
	start();
	fin.seekg(32, ios_base::beg);

	fin.read((char*)&len_record, sizeof(unsigned short));
	fin.read((char*)&frame, sizeof(unsigned short));
	frame = frame - 1;
	fin.read((char*)&line, sizeof(unsigned short));

	fin.read((char*)&sample, sizeof(unsigned short));
	fin.read((char*)&ch, sizeof(unsigned short));

	fin.seekg(24, ios_base::cur);

	fin.read(probe_name, 8);
	probe_name[8] = '\0';
	fin.read((char*)&probe_type, sizeof(unsigned short));
	fin.read((char*)&frq_probe, sizeof(float));
	fin.read((char*)&pole, sizeof(unsigned short));
	fin.read((char*)&wave, sizeof(unsigned short));
	fin.read((char*)&max_angle, sizeof(float));
	fin.read((char*)&offset_from_center, sizeof(float));
	fin.read((char*)&rad_of_cuv, sizeof(float));
	fin.read((char*)&frq_t, sizeof(float));
	fin.read((char*)&frq_r, sizeof(float));
	fin.read((char*)&frq_s, sizeof(float));
	fin.read((char*)&burst, sizeof(unsigned short));
	fin.read((char*)&acq_start, sizeof(float));
	fin.read((char*)&acq_end, sizeof(float));
	fin.read((char*)&line_start, sizeof(unsigned short));
	fin.read((char*)&line_end, sizeof(unsigned short));
	fin.read((char*)&max_beam, sizeof(unsigned short));
	fin.read((char*)&range, sizeof(float));

	fin.seekg(24, ios_base::cur);

	fin.read((char*)&focus_num, sizeof(unsigned short));
	fin.read((char*)&focus_first, sizeof(float));
	fin.read((char*)&PRT, sizeof(unsigned short));
	fin.read((char*)&FR, sizeof(float));

	fin.seekg(352, ios_base::beg);

	fin.read((char*)&RF_size, sizeof(double));
}

void a10::printheader()
{
	cout << "-----RF Data Information-----\n" << "record length:" << len_record << "\n";
	cout << "number of frames:" << frame << "\n";
	cout << "number of lines:" << line << "\n";
	cout << "samples per line:" << sample << "\n";
	cout << "number of channel:" << ch << "\n";
	cout << "probe name:UST-" << probe_name << "\n";
	switch (probe_type)
	{
	case 1:
		cout << "probe type:linear\n";
		break;
	case 2:
		cout << "probe type:convex\n";
		break;
	case 3:
		cout << "probe type:sector\n";
		break;
	case 4:
		cout << "probe type:annular\n";
		break;
	default:
		cout << "probe type:unknown\n";
		break;
	}
	cout << "probe frequency[MHz]:" << frq_probe << "\n";
	cout << "transmit pole:" << pole << "\n";
	cout << "wave pattern:" << wave << "\n";
	cout << "max angle of probe:" << max_angle << "\n";
	cout << "offset from center[mm]:" << offset_from_center << "\n";
	cout << "radius of curvature[mm]:" << rad_of_cuv << "\n";
	cout << "transmit frequency[MHz]:" << frq_t << "\n";
	cout << "receiving frequency[MHz]:" << frq_r << "\n";
	cout << "sampling frequency[MHz]:" << frq_s << "\n";
	cout << "burst cycle:" << burst << "\n";
	cout << "top of ROI[mm]:" << acq_start << "\n";
	cout << "bottom of ROI[mm]:" << acq_end << "\n";
	cout << "beam number of acquire start:" << line_start << "\n";
	cout << "beam number of acquire end:" << line_end << "\n";
	cout << "max number of beams per frame:" << max_beam << "\n";
	cout << "display range[mm]:" << range << "\n";
	cout << "number of focusing:" << focus_num << "\n";
	cout << "first transmit focus[mm]:" << focus_first << "\n";
	cout << "PRT[us]:" << PRT << "\n";
	cout << "frame rate[Hz]:" << FR << "\n";
	cout << "RF data size:" << RF_size << endl;
}

void a10::loadRF()
{
	fin.seekg(360, ios_base::beg);
	fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur);

	RF = vector<vector<vector<vector<short>>>>(frame, vector<vector<vector<short>>>(line,
		vector<vector<short>>(ch, vector<short>(sample - 1, 0))));

	short tmp;
	cout << "loading RF...\n";
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j){
			for (int k = 0; k < ch - 16; ++k){ // back of 80 elements
				fin.seekg(8, ios_base::cur); //attribute 6byte channel number 2byte
				for (int l = 0; l < sample - 1; ++l){
					fin.read((char*)&tmp, sizeof(short));
					RF[i][j][k + 16][l] = tmp - 2048;
				}
			}
			for (int k = 0; k < 16; ++k){ // front of 16 elements
				fin.seekg(8, ios_base::cur);
				for (int l = 0; l < sample - 1; ++l){
					fin.read((char*)&tmp, sizeof(short));
					RF[i][j][k][l] = tmp - 2048;
				}
			}
		}
}

void a10::loadRF0(int frame)
{
	fin.seekg(360, ios_base::beg);
	fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur); //ignore first frame

	if (RF0.empty())
		RF0 = vector<vector<vector<short>>>(line, vector<vector<short>>(ch, vector<short>(sample - 1, 0)));

	for (int i = 0; i < frame; ++i)
		fin.seekg((line * ch * (sample + 3))* sizeof(short), ios_base::cur);

	short tmp;
	cout << "loading RF(frame: " << frame << ")...\n";
	for (int i = 0; i < line; ++i){
		for (int j = 0; j < ch - 16; ++j){ // back of 80 elements
			fin.seekg(8, ios_base::cur); //attribute 6byte channel number 2byte
			for (int k = 0; k < sample - 1; ++k){
				fin.read((char*)&tmp, sizeof(short));
				RF0[i][j + 16][k] = tmp - 2048;
			}
		}
		for (int j = 0; j < 16; ++j){ // front of 16 elements
			fin.seekg(8, ios_base::cur);
			for (int k = 0; k < sample - 1; ++k){
				fin.read((char*)&tmp, sizeof(short));
				RF0[i][j][k] = tmp - 2048;
			}
		}
	}

}

void a10::freeRF()
{
	if (!RF.empty()){
		vector<vector<vector<vector<short>>>>().swap(RF);
		cout << "free RF data!\n";
	}
	else cout << "RF data is empty.\n";
}

void a10::freeRF0()
{
	if (!RF0.empty()){
		vector<vector<vector<short>>>().swap(RF0);
		cout << "free RF data!\n";
	}
	else cout << "RF data is empty.\n";
}

void a10::rmbias()
{
	cout << "removing bias...\n";
	int bias = 0;
	for (int i = 0; i < frame; ++i)
		for (int j = 0; j < line; ++j)
			for (int k = 0; k < ch; ++k){
				bias = accumulate(RF[i][j][k].begin(), RF[i][j][k].end(), 0);
				bias = bias / (sample - 1);
				for (int l = 0; l < sample - 1; ++l)
					RF[i][j][k][l] -= bias;
			}
}

int a10::plotRF0(string dir)
{
	if (RF0.empty()){
		cout << "RF0 is empty.\n";
		return 1;
	}

	ostringstream ost;
	ofstream fout;
	int line = RF0.size();
	int ch = RF0[0].size();
	int sample = RF0[0][0].size();

	for (int i = 60; i < 61; ++i){
		/*ost << "./" << dir << "/" << i << ".dat";*/
		ost << "./" << dir << ".dat";
		fout.open(ost.str(), ios_base::out);
		ost.clear();
		ost.str("");
		for (int j = 0; j < ch; ++j){
			for (int k = 0; k < sample; ++k){
				fout << k << " " << RF0[i][j][k] - 4096 * j << "\n";
			}
			fout << "\n";
		}
		fout.close();
	}

	return 0;
}

vector<vector<double>> a10::calcenv(int frame, float max_angle, float frq_s)
{
	if (RF0.empty()){
		cout << "RF0 is empty.\n";
		vector<vector<double>> dummy;
		return dummy;
	}

	int line = RF0.size();
	int ch = RF0[0].size();
	int sample = RF0[0][0].size() + 1; //fit power of two for FFT

	//spec and buffer setting for FFT
	Ipp8u *specbuff, *initbuff, *workbuff;
	Ipp8u *specbufi, *initbufi, *workbufi;
	int size_specf, size_initf, size_workf;
	int size_speci, size_initi, size_worki;
	IppsFFTSpec_C_64fc *specf = 0;
	IppsFFTSpec_C_64fc *speci = 0;
	Ipp64fc *ipsrc = ippsMalloc_64fc((int)sample);
	Ipp64fc *ipdst = ippsMalloc_64fc((int)sample);
	Ipp64fc *ipsrc2 = ippsMalloc_64fc((int)(4 * sample));
	Ipp64fc *ipdst2 = ippsMalloc_64fc((int)(4 * sample));
	const int fftorder = (int)(log((double)sample) / log(2.0));
	const int ifftorder = (int)(log((double)(4 * sample)) / log(2.0));
	ippsFFTGetSize_C_64fc(fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_specf, &size_initf, &size_workf);
	ippsFFTGetSize_C_64fc(ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &size_speci, &size_initi, &size_worki);
	specbuff = ippsMalloc_8u(size_specf);
	specbufi = ippsMalloc_8u(size_speci);
	initbuff = ippsMalloc_8u(size_initf);
	initbufi = ippsMalloc_8u(size_initi);
	workbuff = ippsMalloc_8u(size_workf);
	workbufi = ippsMalloc_8u(size_worki);
	ippsFFTInit_C_64fc(&specf, fftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbuff, initbuff);
	ippsFFTInit_C_64fc(&speci, ifftorder, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, specbufi, initbufi);

	//calculate delay
	const double c0 = 1540.0;
	int point, add;
	double eledep; //in-bound(um)
	double decimal; //decimal part in sampling point of round trip distance
	vector<double> xi(ch, 0); // x-coordinate of each element
	for (int i = 0; i < ch; ++i)
		xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
		//xi[i] = 0.2 * (i - 47.5) * 1e+3;
	vector<double> theta(line, 0);
	for (int i = 0; i < line; ++i) //beam angle
		theta[i] = max_angle * ((line - 1) / 2 - i) * (M_PI / 180.0);
	vector<double> cendep(sample, 0);
	for (int i = 0; i < sample; ++i)
		cendep[i] = i * (c0 / (2 * frq_s)); //out-bound(um)

	//vector<vector<vector<double>>> elere(line, vector<vector<double>>(ch, vector<double>(4 * sample, 0)));
	//vector<vector<vector<double>>> eleim(line, vector<vector<double>>(ch, vector<double>(4 * sample, 0)));

	vector<vector<double>> RFre(line, vector<double>(sample, 0));
	vector<vector<double>> RFim(line, vector<double>(sample, 0));

	vector<int> addcnt(sample, 0);

	for (int j = 0; j < line; ++j){
		for (int k = 0; k < ch; ++k){
			ippsZero_64fc(ipsrc, sample);
			ippsZero_64fc(ipdst, sample);
			ippsZero_64fc(ipsrc2, 4 * sample);
			ippsZero_64fc(ipdst2, 4 * sample);
			//set
			for (int l = 0; l < sample - 1; ++l){
				ipsrc[l].re = RF0[j][k][l];
				ipsrc[l].im = 0.0;
			}

			//do FFT
			ippsFFTFwd_CToC_64fc(ipsrc, ipdst, specf, workbuff);
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

			for (int l = 0; l < sample; ++l){
				ipsrc2[l].re = ipdst[l].re;
				ipsrc2[l].im = ipdst[l].im;
			}


			//do IFFT
			ippsFFTInv_CToC_64fc(ipsrc2, ipdst2, speci, workbufi);
			ippsZero_8u(workbufi, size_worki);

			//save
			for (int l = 0; l < 4 * sample; ++l){
				ipdst2[l].re *= 4 * sample;
				ipdst2[l].im *= 4 * sample;
			}
			//ippsZero_64fc(ipdst2, 4 * sample);

			for (int l = 0; l < sample; ++l){
				eledep = sqrt(pow(xi[k], 2) + pow(cendep[l], 2) - 2 * xi[k] * cendep[l] * sin(theta[j]));
				point = static_cast<int>(((cendep[l] + eledep) / 2) / (c0 / (8 * frq_s)));
				decimal = ((cendep[l] + eledep) / 2) / (c0 / (8 * frq_s)) - point;
				if (point < 4 * sample - 1){
					RFre[j][l] += ipdst2[point].re + (ipdst2[point + 1].re - ipdst2[point].re) * decimal;
					RFim[j][l] += ipdst2[point].im + (ipdst2[point + 1].im - ipdst2[point].im) * decimal;
					++addcnt[l];
				}
			}
		}
		for (int k = 0; k < sample; ++k){
			if (addcnt[k] != 0){
				RFre[j][k] /= addcnt[k];
				RFim[j][k] /= addcnt[k];
			}
			else{
				RFre[j][k] = 0.0;
				RFim[j][k] = 0.0;
			}
			addcnt[k] = 0.0;
		}
	}



	//vector<vector<vector<double>>>().swap(elere);
	//vector<vector<vector<double>>>().swap(eleim);

	vector<vector<double>> env(line, vector<double>(sample, 0));

	for (int j = 0; j < line; ++j)
		for (int k = 0; k < sample; ++k)
			env[j][k] = sqrt(pow(RFre[j][k], 2) + pow(RFim[j][k], 2));

	vector<vector<double>>().swap(RFre);
	vector<vector<double>>().swap(RFim);

	return env;
}





physio::physio(string filename) : file(filename)
{
	fn = filename;
}

physio::~physio()
{
}

int physio::extract(int offset)
{
	vector<unsigned char> check;
	char tmp;
	while (!fin.eof()){
		fin.read((char*)&tmp, sizeof(unsigned char));
		check.push_back(tmp);
	}
	fin.close();

	/*searching key*/
	list<unsigned char> key = { 0xFF, 0x53, 0x21, 0x10, 0x4F, 0x42, 0x00, 0x00 };
	vector<unsigned char>::iterator it = find_end(check.begin(), check.end(), key.begin(), key.end());
	int dist = distance(check.begin(), it);
	if (it == check.end()){
		cout << "not found key\n";
		return 1;
	}
	else{
		cout << "key is at check[" << dist << "]\n";
	}
	int size = (check.size() - dist - 12 - 1) / 32;
	vector<unsigned char>().swap(check);

	/*extract physio data*/
	open(fn);
	fin.seekg(dist + 12, ios_base::beg);

	short tmp2 = 0;

	ECG.reserve(size);
	PCG_min.reserve(size);
	PCG_max.reserve(size);
	for (int i = 0; i < size; ++i){
		fin.seekg(8, ios_base::cur);
		fin.read((char*)&tmp2, sizeof(short));
		ECG.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_min.push_back(tmp2 & 0x03FF);
		fin.read((char*)&tmp2, sizeof(short));
		PCG_max.push_back(tmp2 & 0x03FF);
		fin.seekg(18, ios_base::cur);
	}

	ECG.erase(ECG.begin(), ECG.begin() + offset);
	PCG_min.erase(PCG_min.begin(), PCG_min.begin() + offset);
	PCG_max.erase(PCG_max.begin(), PCG_max.begin() + offset);

	return 0;
}

void physio::write()
{
	ofstream pecg("./ECG.dat", ios_base::out);
	ofstream ppcgmin("./PCG_min.dat", ios_base::out);
	ofstream ppcgmax("./PCG_max.dat", ios_base::out);

	float ratio = 1000 / 998;
	for (int i = 0; i < ECG.size(); ++i){
		pecg << i * ratio << " " << ECG[i] << "\n";
		ppcgmin << i * ratio << " " << PCG_min[i] << "\n";
		ppcgmax << i * ratio << " " << PCG_max[i] << "\n";
	}
	pecg.close();
	ppcgmin.close();
	ppcgmax.close();

	cout << "wrote physio data!\n";
}