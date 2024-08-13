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

class Demo_PM_BP : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_PM_BP(int index);

	~Demo_PM_BP();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(const DemodulationInitParamater& info);

	void PM_BPInit();

	void InitSymbolSync();

	void InitCostasPLL();
	void PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize,float* fPLLNCO, float* fPLLPastFreqPart,float* PLL_Freq_Part_2, float* PLL_Phase_Part_2);

	void InitCostasPLL_bpsk();
	void PLLCostas_bpsk(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart);

	void Judgment(Complex buffin[], int bufflen, char buffout[]);

	void downphase(float BuffIn[], Complex BuffOut[], int nBuffSize, float* fDownConversionFc, double* fOrthDownConversionPhase);

	float mean_function(float* data, int start, int end);

	// add on 20240730
	void est_pmbp_rb(Complex* data, int len, int fftlen, float fs, float* est_rb);
	void ippfft(float* data_real, float* data_imag, int fftLength, float* output);
	int findmax1(float* a, int p, float FFT_size, int start);
	void Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout);
	void dgz_sort_pmbpqp(float arr[], int index[], int n);
	void downphaseFrequence(Complex* data_input, float phase, Complex* downout, int dataLength);

public:
	SAlgDemInit m_sAlgDemInit;							//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
	CostasPLL CostasPll2;
	CostasPLL CostasPll3;

	CostasPLL CostasPll_bpsk;
	CostasPLL CostasPll2_bpsk;
	CostasPLL CostasPll3_bpsk;

private:
	int nSilceCount = 0;								//当前处理的信号所在第几片
	int m_SamleSize;
	int m_Rb;
	float fc2 = 0;
	float fAGCPastVc = 0.5;								//AGC中的VC值缓存

	double fOrthDownConversionPhase = 0;				//正交下变频中的相位值
	int nDownSampClock;									//降采样时的时钟控制  
	float fDownConversionPhase;							//下变频中的相位值
	float fPLLNCO;										//PM锁相环中的本地NCO
	float fPLLNCO_1;									//BPSK锁相环中的本地NCO
	float fPLLPastFreqPart;								//PM锁相环中的频率跟踪曲线
	float fPLLPastFreqPart_1;							//bpsk锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	Complex* cFilterBuff;								//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;							//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;								//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;

	float Judgth1 = 0;
	int m_index;
	Algorithm_Demodulation PM_BP;

	// add on 20240730
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int PMBP_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize

	float fcsub;                    // 副载波
	int fs;
	int band;
	int flag1 = 0;					// 码速估计标志位		




};