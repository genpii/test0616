#include "stdafx.h"
#include "header.h"

/* Intel IPP */
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

void SimpleMovingAverage(vector<float> &src, int N){
	int length = src.size();
	int half = (N - 1) / 2;
	float tmp;
	vector<float> dst(length, 0);

	for (int i = 0; i < half; ++i){ //ëOï˚ã´äE
		tmp = 0.0;
		for (int j = -i; j <= half + i; ++j)
			tmp += src[i + j];
		dst[i] = tmp / static_cast<float>(half + i + 1);
	}
	for (int i = half; i < length - half; ++i){ //îÒã´äE
		tmp = 0.0;
		for (int j = -half; j <= half; ++j)
			tmp += src[i + j];
		dst[i] = tmp / static_cast<float>(N);
	}
	for (int i = length - half; i < length; ++i){ //å„ï˚ã´äE
		tmp = 0.0;
		for (int j = -half; j < length - i; ++j)
			tmp += src[i + j];
		dst[i] = tmp / static_cast<float>(length - i + half);
	}

	src = dst;
}

void HanningMovingAverage(vector<float> &src, int N){
	int length = src.size();
	int half = (N - 1) / 2;
	float tmp;
	vector<float> dst(length, 0);
	vector<float> han(N, 0);
	for (int i = 0; i < N; ++i)
		han[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / (N - 1));
	
	for (int i = 0; i < half; ++i){ //ëOï˚ã´äE
		tmp = 0.0;
		for (int j = -i; j <= half; ++j)
			tmp += src[i + j] * han[j + half];
		dst[i] = tmp / static_cast<float>(half + i + 1);
	}
	for (int i = half; i < length - half; ++i){ //îÒã´äE
		tmp = 0.0;
		for (int j = -half; j <= half; ++j)
			tmp += src[i + j] * han[j + half];
		dst[i] = tmp / static_cast<float>(N);
	}
	for (int i = length - half; i < length; ++i){ //å„ï˚ã´äE
		tmp = 0.0;
		for (int j = -half; j < length - i; ++j)
			tmp += src[i + j] * han[j + half];
		dst[i] = tmp / static_cast<float>(length - i + half);
	}

	src = dst;
}

vector<int> PeakDetection(const vector<float> &src){
	vector<int> dst;
	for (int i = 1; i < src.size() - 1; ++i){
		if (src[i] > src[i - 1] && src[i] > src[i + 1])
			dst.push_back(i);
	}
	return dst;
}

vector<int> PeakDetection(const vector<float> &src, int N){
	vector<int> dst;
	int cut = N;
	for (int i = 1; i < src.size() - 1; ++i){
		if (src[i] > src[i - 1] && src[i] > src[i + 1])
			dst.push_back(i);
	}

	auto result = remove_if(dst.begin(), dst.end(), [cut](int x){return cut >= x; });
	dst.erase(result, dst.end());

	return dst;
}

void PlotVector(const vector<float> &src, string &str){
	ofstream fout(str, ios_base::out);
	for (int i = 0; i < src.size(); ++i){
		fout << i << " " << src[i] << "\n";
	}
	fout.close();
}

void PlotVector(const vector<float> &src, float scale, string &str){
	ofstream fout(str, ios_base::out);
	for (int i = 0; i < src.size(); ++i){
		fout << static_cast<float>(i) * scale << " " << src[i] << "\n";
	}
	fout.close();
}


vector<int> MakeDelayProfile(const vector<vector<float>> &re, const vector<vector<float>> &im, int refcent){
	int ch = re.size();
	int sample = re[0].size();
	int ref = 47; //47th element is reference
	vector<int> dst(ch, 0);

	//ëää÷éZèoópÇÃê›íË
	IppStatus status;
	const int csrc1Len = 256; //ëää÷ëãïù
	const int csrc2Len = sample; //ÇﬂÇÒÇ«Ç≠Ç≥Ç¢Ç©ÇÁêMçÜëSïîì¸ÇÍÇøÇ·Ç§
	int lowLag = -256; // -lowLag * 2 + 1 Ç™ëSíxâÑó 
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
	pbuffer = ippsMalloc_8u(bufsize);

	//ÉäÉtÉ@ÉåÉìÉXÇÃê›íË
	ippsZero_32fc(csrc1, csrc1Len);
	int csrc1start = refcent - csrc1Len / 2 + 1;
	for (int i = 0; i < csrc1Len; ++i){
		csrc1[i].re = re[ref][csrc1start + i];
		csrc1[i].im = im[ref][csrc1start + i];
	}
	csrc1Norm = 0.0;
	ippsNorm_L2_32fc64f(csrc1, csrc1Len, &csrc1Norm);

	//ëää÷éZèo
	for (int i = 0; i < ch; ++i){
		ippsZero_32fc(csrc2, csrc2Len);
		ippsZero_32fc(cdst, cdstLen);
		ippsZero_32f(cdstMag, cdstLen);
		ippsZero_32f(cdstPha, cdstLen);
		ippsZero_64f(csrc2Norm, cdstLen);
		ippsZero_32f(csrcNorm, cdstLen);
		ippsZero_32fc(csrcNormCplx, cdstLen);

		for (int j = 0; j < csrc2Len; ++j){
			csrc2[j].re = re[i][j];
			csrc2[j].im = im[i][j];
		}
		status = ippsCrossCorrNorm_32fc(csrc1, csrc1Len, csrc2, csrc2Len, cdst, cdstLen, csrc1start + lowLag, NormA, pbuffer);

		for (int j = 0; j < cdstLen; ++j){
			csrc2Normtmp = 0.0;
			ippsNorm_L2_32fc64f(csrc2 + csrc1start + lowLag + j, csrc1Len, &csrc2Normtmp);
			csrc2Norm[j] = csrc2Normtmp;
		}
		ippsMulC_64f_I(csrc1Norm, csrc2Norm, cdstLen);
		ippsConvert_64f32f(csrc2Norm, csrcNorm, cdstLen);
		ippsRealToCplx_32f(csrcNorm, cdstLenzeros, csrcNormCplx, cdstLen);
		ippsDiv_32fc_I(csrcNormCplx, cdst, cdstLen);
		ippsMagnitude_32fc(cdst, cdstMag, cdstLen);
		ippsPhase_32fc(cdst, cdstPha, cdstLen);
		ippsMaxIndx_32f(cdstMag, cdstLen, &cmax, &cmaxidx);

		dst[i] = cmaxidx + lowLag;
	}

	return dst;

}

/* Eigen */