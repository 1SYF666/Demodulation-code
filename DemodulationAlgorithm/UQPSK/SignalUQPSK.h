#pragma once
#include "distinguish.h"
#include "../Algorithm_Demodulation.h"
#include "../../demoduletion_factroy.h"
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



class Demo_UQPSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;


public:

	Demo_UQPSK(int index);

	~Demo_UQPSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_UQPSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(const DemodulationInitParamater& info);

	void UQPSKInit();

	void InitSRC();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[],int bufflen, char buffout[]);

	void ippfft(float* data_real, float* data_imag, int fftLength, float* output);

	void fft_shift(float* fft_s, int FFT_size);

	int getFFTmax(float* a, int start, int end);

	// add on 20240726
	void PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart);

	void estuqpsk(float* s_2_real_temp, int len_temp, float* rb1, float* rb2, int fs, int fftLength);

	float mean(float* data, int start, int end);

	void Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout);

public:
	SAlgDemInit m_sAlgDemInit;									//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;						// 对应大码速那部
	SymbolSyncFactor symbolsyncsactorInit2;						// 对应小码速那部
	CFILTERPARAM CFilter;
	CBLOCKFILTERPARAM CBlockFilter;
	SymbolSyncBuffer SySyncBufferI;
	SymbolSyncBuffer SySyncBufferQ;
	CostasPLL CostasPll;
	CostasPLL CostasPll2;
	CostasPLL CostasPll3;
	SRCPARAM SRCParam;

private:

	int fs;                                                     //采样率
	int nSilceCount = 0;										//当前处理的信号所在第几片
	int m_SamleSize;
	int m_Rb;
	int m_RbQ;
	float fAGCPastVc =0.5;										//AGC中的VC值缓存
	float fc[2] = { 0 };
	double fOrthDownConversionPhase = 0;						//正交下变频中的相位值
	int nDownSampClock;											//降采样时的时钟控制  
	float fDownConversionPhase;									//下变频中的相位值
	float fPLLNCO;												//锁相环中的本地NCO
	float fPLLPastFreqPart;										//锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	Complex* cFilterBuff;										//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;									//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;										//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;
	float rbEst1=0;
	float rbEst2=0;
	int flag1 = 0;
	Algorithm_Demodulation De_UQPSK;
	Algorithm_Demodulation De_UQPSK2;							// add on 20240727
	int m_index;

	// add on 20240621
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int UQPSK_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize


};