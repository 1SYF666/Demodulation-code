#pragma once
#include "distinguish.h"
#include "../Algorithm_Demodulation.h"
#include "../../demoduletion_factroy.h"
#include"../../Interface/computing_interface.h"
#include <tchar.h>
#include <iostream>
#include <vector>
#include <Windows.h>
#include <iostream>
#include <stdlib.h>
#include "math.h"
#include "omp.h"
#include "ipp.h"
#pragma comment(lib,"ippsmt.lib")
#pragma comment(lib,"ippcoremt.lib")
#pragma comment(lib,"ippvmmt.lib")
#include <QMutex>
#include <QDateTime>
#include <QFile>
#include <QByteArray>

using namespace std;

class Demo_CDMA : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_CDMA(int index);

	~Demo_CDMA();

	void max_index_array(float* data, int len, float* output);

	void ippfft(float* data_real, float* data_imag, int FFT_LENGTH1, int flag1, float* Ioutput, float* Qoutput, float* output);

	void ippifft(float* data_real, float* data_imag, int FFT_LENGTH1, float* ifft_real, float* ifft_imag, float* Module);

	void Capture(float* inputI, float* inputQ);

	void Tracking(float* inputI, float* inputQ, int symbol_num, Complex* demoResult);

	void Tracking_new(float* inputI, float* inputQ, int symbol_num, Complex* PNsyndata);

	void Initial_Tracking();

	void InitCostasPLL();

	void CDMADemo(float* input_I, float* input_Q, int signLen, DemodulationResult* signal_16APSK_Out);

	void InitialDemodulation(const DemodulationInitParamater& signalinit);

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	float mean_function(float* data, int start, int end);

	void AGC(Complex* BuffIn, float target_power);

	float meanArray(float* data, int len);

	// add on 20240622
	void AGCCDMA(Complex* BuffIn, float target_power, int len);

	void InitBlockCFilter();

	void Initial_Tracking_new();

	// add on 20240803
	void estCDMASpreadCode(float* input_I, float* input_Q, float fs,int pnEstLen,int* pncode);
	int maxIndex(float* data, int start, int end);
	int minIndex(float* data, int start, int end);
	float max_function(float* data, int start, int end);
public:
	CostasPLL CostasPll;
	CFILTERPARAM CFilter;
private:
	Complex* cFilterBuff;
	int pnLength;				// PN码长度
	float fs;					// 采样率
	int* pnCode;				// PN码序列
	int m_index;
	SAlgDemInit m_sAlgDemInit;  // 算法初始化结构

	float fAGCPastVc = 0.5;		//AGC中的VC值缓存

	// add on 20240429
	int flag = 0;               // 缓存数据标志位
	int slice = 0;             // 扩频解调总片数
	int m_SamleSize = 0;     // 其它信号解调一次点数为8192
	int CDMA_Samplesize = 0; // 扩频解调一次总点数
	Complex* Databuf;


	int chip_sps;				// 码片采样点数
	int flag_acqusiation;   // 捕获标志位 1--进行捕获、0--不进行捕获
	int* delayDot_buf;			// 捕获延迟点数及载频寄存器 0--载频；1--延迟点
	int capture_fc;			// 捕获载频值
	int pncodelength_sps;       // PN码采样后长度
	int fc_acqusiation;         // 最大频差
	int step_acqusiation;		// 扫频步长
	int step_fc_acqusiation;    // 扫频值
	int *fc_acq_buff;           // 扫频寄存器s
	float* PNcode_est_sps;      // PN码采样寄存器
	//float* PNfftdataI;			// PN码fft结果
	//float* PNfftdataQ;
	//float* PNfftdata;
	//float* CapturedataI;       // 捕获信号数据
	//float* CapturedataQ;
	//float* SignalFFTI;			// 捕获信号数据fft	
	//float* SignalFFTQ;
	//float* Signaldata;
	//float* DotMultI;			// 捕获自相关频谱
	//float* DotMultQ;
	//float* SignalifftI;
	//float* SignalifftQ;
	//float* Signaldata_temp;



	Algorithm_Demodulation De_CDMA;
	ComputingInterface Computer;    // add on 20240807

	struct TrackInit
	{
		float ThresholdPositive;
		float ThresholdNegative;
		float TrackingLoopCoef1;
		float TrackingLoopCoef2;
		float TrackingLoopCoef_new1;
		float TrackingLoopCoef_new2;
	};

	TrackInit trackcoef;

	int sum_temp1;				// 扩频信号第一次积分点数
	int sum_temp2;             // 扩频信号第二次积分点数

	float fPLLNCO;								//锁相环中的本地NCO
	float fPLLPastFreqPart;						//锁相环中的频率跟踪曲线
	FreqLLBuffer FLLBuffer;

	// add on 2024/05/06
	int* quotient;						// 商 
	int* remainder;						// 余数 
	int flag_trcking = 0;
	float* data_track_temp_I;			// 跟踪数据缓存器
	float* data_track_temp_Q;			// 跟踪数据缓存器


	// add on 20240803

	int flag1 = 1;                     // 扩频序列估计标志位



};




