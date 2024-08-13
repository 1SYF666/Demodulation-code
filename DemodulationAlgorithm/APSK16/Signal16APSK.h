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


class Demo_16APSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;
	
	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;
public:

	Demo_16APSK(int index);

	~Demo_16APSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void APSKInit();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[], float power[], int bufflen, int buffout[]);

	float mean_function(float* data, int start, int end);

	void diff_code_16apsk(int* y, char* c_symbol,int len);

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	float GetThreshold(Complex Buffin[], int bufflen, float Buffout[]);

	float max_function(float* data, int start, int end);

public:
	SAlgDemInit m_sAlgDemInit;					//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 1;						//当前处理的信号所在第几片
	int m_SamleSize;
	int m_Rb;

	float fAGCPastVc =0.5;						//AGC中的VC值缓存
	float mean_sort_theta_1 = 0;

	double fOrthDownConversionPhase = 0;		//正交下变频中的相位值
	int nDownSampClock;							//降采样时的时钟控制  
	float fDownConversionPhase;					//下变频中的相位值

	float fPLLNCO;								//锁相环中的本地NCO
	float fPLLPastFreqPart;						//锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	Complex* cFilterBuff;						//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;					//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;						//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;

	float Judgth1 = 0;
	int difftemp = 0;
	Algorithm_Demodulation De_APSK;

	int m_index;

	// add on 20240514
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int APSK_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize
};
