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



/* Eigen */