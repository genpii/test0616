#include "stdafx.h"
#include <opencv2/opencv.hpp>
#include <opencv2/opencv_lib.hpp>
#include "header.h"

using namespace std;
using namespace cv;

void BLinear(a10 &raw){ //参照でa10にアクセス

	int line = raw.line;
	int sample = raw.sample;
	int ch = raw.ch;
	int fs = raw.frq_s;
	int linescale = 2;
	float pitch = 0.2;
	float c = 1540.0;
	float itvl = c / fs / 2.0 / 1000; //[mm/sample]
	float lppx = pitch / linescale; //[mm/px]

	vector<float> xi(ch, 0); // x-coordinate of each element
	for (int i = 0; i < ch; ++i)
		xi[i] = 0.2 * (47.5 - i) * 1e+3; //um
	vector<float> cendep(sample, 0);
	for (int i = 0; i < sample; ++i)
		cendep[i] = i * (c / (2 * fs));
	int point, add;
	float eledep; //in-bound(um)
	float decimal; //decimal part in sampling point of round trip distance

	int row = sample * itvl / lppx + 1;
	int col = line * linescale;
	float gain = 30.0;
	int color;

	vector<float> re(sample, 0);
	vector<float> im(sample, 0);
	vector<float> env(sample, 0);
	vector<int> addcnt(sample, 0);

	vector<int> chbegin(line, 0);
	vector<int> chend(line, ch);
	for (int i = 0; i < line; ++i){
		if (i < 42)
			chbegin[i] = 42 - i;

		if (i > 138)
			chend[i] = ch - (i - 138);
	}

	Mat dst;
	dst = Mat::zeros(row, col, CV_8UC1);
	for (int i = 0; i < line; ++i){
		raw.generate_AS(i);

		for (int j = chbegin[i]; j < chend[i]; ++j){
			for (int k = 0; k < sample; ++k){
				eledep = sqrt(pow(xi[j], 2) + pow(cendep[k], 2));
				point = static_cast<int>(((cendep[k] + eledep) / 2) / (c / (8 * fs)));
				decimal = ((cendep[k] + eledep) / 2) / (c / (8 * fs)) - point;
				if (point < 4 * sample - 1){
					re[k] += raw.ele0re[j][point] + (raw.ele0re[j][point + 1] - raw.ele0re[j][point]) * decimal;
					im[k] += raw.ele0im[j][point] + (raw.ele0im[j][point + 1] - raw.ele0im[j][point]) * decimal;
					++addcnt[k];
				}
			}
		}
		for (int k = 0; k < sample; ++k){
			if (addcnt[k] != 0){
				re[k] /= addcnt[k];
				im[k] /= addcnt[k];
			}
			else{
				re[k] = 0.0;
				im[k] = 0.0;
			}
			addcnt[k] = 0;

			env[k] = sqrt(pow(re[k], 2) + pow(im[k], 2));
		}

		float max = *max_element(env.begin(), env.end());

		for (int j = 0; j < row; ++j){
			int index = static_cast<int>((float)j / row * sample);
			//cout << index << endl;
			color = static_cast<int>(256 * (20 * log10(env[index] / max) + gain) / gain);
			if (color > 255) color = 255;
			if (color < 0) color = 0;
			for (int k = 0; k < linescale; ++k){
				dst.at<unsigned char>(j, col - (i * linescale + k) - 1) = color;
			}
		}
		
	}
	imwrite("B.png", dst);
	imshow("B-mode", dst);
	waitKey(0);

}
