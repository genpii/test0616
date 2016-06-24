// stdafx.h : 標準のシステム インクルード ファイルのインクルード ファイル、または
// 参照回数が多く、かつあまり変更されない、プロジェクト専用のインクルード ファイル
// を記述します。
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>


// TODO: プログラムに必要な追加ヘッダーをここで参照してください。
#define _USE_MATH_DEFINES
#include <iostream> //コンソール入出力
#include <string> //string文字列
#include <fstream> //ファイル入力
#include <vector> //配列ベクトル用
#include <valarray> //多次元配列
#include <numeric> //演算
#include <complex> //複素数
#include <ipp.h> //IPP 信号処理
#include <ippm.h> //IPP 行列
#include <list> //list
#include <algorithm> //アルゴリズム
#include <direct.h> //階層走査 mkdir用
#include <array> //C++ STL配列
#include <sstream> //stringstream
#include <map>
//#include <mkl.h> //IntelMKL
#include "Eigen\core"
#include "unsupported\Eigen\NonLinearOptimization"