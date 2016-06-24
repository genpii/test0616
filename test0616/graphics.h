#pragma once
#include <windows.h>
#include <math.h>

HINSTANCE hInstance = GetModuleHandle(0);
HWND hwnd;
HDC hdc;
PAINTSTRUCT ps;
HPEN hpen;
HBRUSH brush;
int x_coordinate=0, y_coordinate=0;

int get_x(), get_y();

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {

	switch (msg) {
	default:
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	case WM_MOUSEMOVE:
		x_coordinate = get_x();
		y_coordinate = get_y();
		break;
	}
	return DefWindowProc(hwnd, msg, wp, lp);
}

int pageb(int x1, int y1, int nx, int ny, char *ch) {

	WNDCLASS winc;

	winc.style = CS_HREDRAW | CS_VREDRAW;
	winc.lpfnWndProc = WndProc;
	winc.cbClsExtra = winc.cbWndExtra = 0;
	winc.hInstance = hInstance;
	winc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	winc.hCursor = LoadCursor(NULL, IDC_ARROW);
	winc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	winc.lpszMenuName = NULL;
	winc.lpszClassName = TEXT("window");

	if (!RegisterClass(&winc)) return 0;

	hwnd = CreateWindow(TEXT("window"), ch, WS_OVERLAPPEDWINDOW|WS_VISIBLE, x1, y1, nx, ny, NULL, NULL, hInstance, NULL);

	if (hwnd == NULL) return 0;
	
	return 1;
}

void gspt(int x1, int y1, COLORREF color) {

	SetPixel(hdc, x1, y1, color);
}

void gsline(int x1, int y1, int x2, int y2, COLORREF color) {

	hpen = CreatePen(PS_SOLID, 0, color);
	SelectObject(hdc, hpen);
	MoveToEx(hdc, x1, y1, NULL);
	LineTo(hdc, x2, y2);
	DeleteObject(hpen);
}

void gsline_d(int x1, int y1, int x2, int y2, COLORREF color) {

	hpen = CreatePen(PS_DASH, 0, color);
	SelectObject(hdc, hpen);
	MoveToEx(hdc, x1, y1, NULL);
	LineTo(hdc, x2, y2);
	DeleteObject(hpen);
}

void gsrect(int x1, int y1, int x2, int y2, COLORREF color_in, COLORREF color_out) {

	hpen = CreatePen(PS_SOLID, 0, color_out);
	SelectObject(hdc, hpen);
	brush = CreateSolidBrush(color_in);
	SelectObject(hdc, brush);
	Rectangle(hdc, x1, y1, x2, y2);
	DeleteObject(hpen);
	DeleteObject(brush);
}

void gsarc(int x0, int y0, float r, float theta1, float theta2, int thickness, COLORREF color) {

	hpen = CreatePen(PS_SOLID, thickness, color);
	SelectObject(hdc, hpen);
	float temp = theta1;
	theta1 = theta2;
	theta2 = temp + 360;
	int x1 = (int)round(x0 - r);
	int y1 = (int)round(y0 - r);
	int x2 = (int)round(x0 + r);
	int y2 = (int)round(y0 + r);
	int x3 = (int)round(x0 + r*cos(theta1*M_PI / 180));
	int y3 = (int)round(y0 + r*sin(theta1*M_PI / 180));
	int x4 = (int)round(x0 + r*cos(theta2*M_PI / 180));
	int y4 = (int)round(y0 + r*sin(theta2*M_PI / 180));
	Arc(hdc, x1, y1, x2, y2, x3, y3, x4, y4);
	DeleteObject(hpen);
}

COLORREF gscol256(int ir, int ig, int ib) {

	COLORREF color;

	color = RGB(ir, ig, ib);
	return color;
}

void ptext(int x1, int y1, int len, char *ch) {

	SetTextColor(hdc, RGB(0, 0, 0));
	TextOut(hdc, x1, y1, ch, len);
}

int get_x() {

	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(hwnd, &pt);
	return pt.x;
}

int get_y() {

	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(hwnd, &pt);
	return pt.y;
}

void cursorpos(int x, int y) {

	POINT pt;
	pt.x = x;
	pt.y = y;
	ClientToScreen(hwnd, &pt);
	SetCursorPos(pt.x, pt.y);
}

void eraseg() {

	EndPaint(hwnd, &ps);
	InvalidateRect(hwnd, NULL, TRUE);
	hdc = BeginPaint(hwnd, &ps);
}

int messageloop(void) {

	WNDCLASS winc;
	MSG msg;

	//winc.lpfnWndProc = WndProc;
	while (GetMessage(&msg, NULL, 0, 0)>0) {
		DispatchMessage(&msg);
	}
	return msg.wParam;
}

void start(void) {

	hdc = BeginPaint(hwnd, &ps);
}

void end(void) {

	EndPaint(hwnd, &ps);
}

extern BOOL funcSaveRect(LPCTSTR lpFname, HDC hDC, LONG cx, LONG cy, LONG sx, LONG sy)
{
	// windowのフレーム分だけ補正
	sx = sx - 16;
	sy = sy - 39;

	HANDLE hFile = CreateFile(lpFname, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

	if (hFile != INVALID_HANDLE_VALUE) {
		LONG    lHeadSize = (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFO));
		LONG    lWidthSize = (sx * sizeof(DWORD));
		LONG    lImageSize = (lWidthSize * sy);
		DWORD   dwSize;

		// BITMAPFILEHEADERの初期化
		BITMAPFILEHEADER bmpHead = { 0 };
		bmpHead.bfType = 0x4D42;       // 識別子(BM)
		bmpHead.bfSize = lHeadSize + lImageSize;
		bmpHead.bfReserved1 = 0;
		bmpHead.bfReserved2 = 0;
		bmpHead.bfOffBits = lHeadSize;

		// BITMAPINFOの初期化
		BITMAPINFO bmpInfo = { 0 };
		bmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmpInfo.bmiHeader.biWidth = sx;
		bmpInfo.bmiHeader.biHeight = sy;
		bmpInfo.bmiHeader.biPlanes = 1;
		bmpInfo.bmiHeader.biBitCount = 32;
		bmpInfo.bmiHeader.biCompression = BI_RGB;
		bmpInfo.bmiHeader.biSizeImage = 0;
		bmpInfo.bmiHeader.biXPelsPerMeter = 0;
		bmpInfo.bmiHeader.biYPelsPerMeter = 0;
		bmpInfo.bmiHeader.biClrUsed = 0;
		bmpInfo.bmiHeader.biClrImportant = 0;

		// DIBセクションの作成
		LPDWORD     lpPixel;    // ピクセル配列
		HBITMAP     hBitmap;    // ビットマップ
		HDC         hSaveDC;    // 保存スクリーン
		hBitmap = CreateDIBSection(NULL, &bmpInfo, DIB_RGB_COLORS, (LPVOID*)&lpPixel, NULL, 0);
		hSaveDC = CreateCompatibleDC(hDC);
		SelectObject(hSaveDC, hBitmap);

		// 保存領域のコピー
		BitBlt(hSaveDC, 0, 0, sx, sy, hDC, cx, cy, SRCCOPY);

		// ファイルに書き込む
		WriteFile(hFile, &bmpHead, sizeof(BITMAPFILEHEADER), &dwSize, NULL);
		WriteFile(hFile, &bmpInfo, sizeof(BITMAPINFO), &dwSize, NULL);
		WriteFile(hFile, lpPixel, lImageSize, &dwSize, NULL);

		// DIBセクションの破棄
		DeleteDC(hSaveDC);
		DeleteObject(hBitmap);
		CloseHandle(hFile);
		return TRUE;
	}
	return FALSE;
}