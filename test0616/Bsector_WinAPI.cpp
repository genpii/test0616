#include "stdafx.h"
#include "graphics.h"


using namespace std;

void BSector2(const vector<vector<float>>& env, float dangle, float fs)
{
	float gain = 50.0;
	COLORREF gscolor;
	
	// calculate max value for normalization
	int line = static_cast<int>(distance(env.begin(), env.end()));
	int sample = static_cast<int>(distance(env[0].begin(), env[0].end()));
	vector<float> maxvec(line, 0);
	vector<int>::iterator it;
	int tmp;
	for (int i = 0; i < line; ++i){
		maxvec[i] = *max_element(env[i].begin(), env[i].end());
	}
	float max = *max_element(maxvec.begin(), maxvec.end());
	//cout << "frame:" << frm << " " << max << "\n";
	//max = 1.99135e007;
	float deltad = 1540.0 / (2.0 * (fs * 1e6));

	int wwB = 700;
	unsigned long wid = pageb(0, 0, 1000, wwB, "B-mode");
	start();

	int Bx_ini = 500;
	int By_ini = 30;
	float ini_d = 0.0;
	float b_angld = (2.0 * M_PI) * dangle / 360.; //20151117
	int dadd = static_cast<int>(4. * 240. * b_angld / (2.0 * M_PI));
	int ini_dn = static_cast<int>(1e-3 * ini_d / deltad);
	//nstp = (ncount + ini_dn) / By_ini;
	int nstp = (sample + ini_dn) / 650;
	++nstp;

	int n_angl, ic;
	float b_angl;
	float addmin = -120;
	for (int i = 0; i < line; ++i){
		n_angl = 2 * (sample + ini_dn) * sin(b_angld) / 2;
		for (int j = 0; j < sample; j += nstp){
			ic = static_cast<int>(256*(20.*log10(env[i][j] / max) + gain) / gain);
			if (ic < 0) ic = 0;
			if (ic > 255) ic = 255;
			gscolor = gscol256(ic, ic, ic);
			for (int k = 0; k < n_angl; ++k){
				b_angl = 0.125 * (2.0 * M_PI) * (addmin + dadd * (line - 1 - i) - 0) / 120. + k * b_angld / (float)n_angl;
				gsline(Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)),
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)),
					Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)) + 1,
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)), gscolor);
				gspt(Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)) + 1,
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)), gscolor);
			}
		}
	}
	char str[128];
	int frm = 0;
	sprintf(str, "./%d.jpg", frm);
	funcSaveRect(str, hdc, 0, 0, 1000, 700);
	end();
	InvalidateRect(hwnd, NULL, 1);


}

float BSector22(const vector<vector<double>>& env, float dangle, float fs, int frm)
{
	float gain = 60.0;
	COLORREF gscolor;

	// calculate max value for normalization
	int line = static_cast<int>(distance(env.begin(), env.end()));
	int sample = static_cast<int>(distance(env[0].begin(), env[0].end()));
	vector<double> maxvec(line, 0);
	vector<int>::iterator it;
	int tmp;
	for (int i = 0; i < line; ++i){
		maxvec[i] = *max_element(env[i].begin(), env[i].end());
	}
	double max = *max_element(maxvec.begin(), maxvec.end());
	cout << "frame:" << frm << " " << max << "\n";
	//max = 1.99135e007;
	double deltad = 1540.0 / (2.0 * (fs * 1e6));

	int wwB = 700;
	unsigned long wid = pageb(0, 0, 1000, wwB, "B-mode");
	start();

	int Bx_ini = 500;
	int By_ini = 30;
	double ini_d = 0.0;
	double b_angld = (2.0 * M_PI) * dangle / 360.; //20151117
	int dadd = static_cast<int>(4. * 240. * b_angld / (2.0 * M_PI));
	int ini_dn = static_cast<int>(1e-3 * ini_d / deltad);
	//nstp = (ncount + ini_dn) / By_ini;
	int nstp = (sample + ini_dn) / 650;
	++nstp;

	int n_angl, ic;
	double b_angl;
	double addmin = - 10 - line * dangle;
	for (int i = 0; i < line; ++i){
		n_angl = 2 * (sample + ini_dn) * sin(b_angld) / 2;
		for (int j = 0; j < sample; j += nstp){
			ic = static_cast<int>(256 * (20.*log10(env[i][j] / max) + gain) / gain);
			if (ic < 0) ic = 0;
			if (ic > 255) ic = 255;
			gscolor = gscol256(ic, ic, ic);
			for (int k = 0; k < n_angl; ++k){
				b_angl = 0.125 * (2.0 * M_PI) * (addmin + dadd * (line - 1 - i) - 0) / 120. + k * b_angld / (float)n_angl;
				gsline(Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)),
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)),
					Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)) + 1,
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)), gscolor);
				gspt(Bx_ini + (int)(((j + ini_dn) / nstp)*sin(b_angl)) + 1,
					By_ini + (int)(((j + ini_dn) / nstp)*cos(b_angl)), gscolor);
			}
		}
	}
	char str[128];
	sprintf(str, "phantom0202.jpg");
	//funcSaveRect(str, hdc, 0, 0, 1000, 700);
	end();
	InvalidateRect(hwnd, NULL, 1);
	scanf("%s", str);
	return max;
}