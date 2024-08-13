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

class Demo_64QAM :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:

	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_64QAM(int index);

	~Demo_64QAM();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	void InitialDemodulation(int Rb);

	void QAM64_Init();

	void InitSymbolSync();

	void InitCostasPLL(float* fBLcoef);

	void InitBlockCFilter();

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	void Judgment_64qam(Complex buffin[], int bufflen, int buffout[]);
	
	float mean_function(float* data, int start, int end);
		
	void diff_code_64qam(int* y, char* c_symbol, int len);

	void AGCQAM64(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);

public:
	SAlgDemInit m_sAlgDemInit;						//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 0;							//第几次信号解调标志位
	int m_SamleSize;
	int m_Rb;

	float fAGCPastVc = 0.3;							//AGC中的VC值缓存

	double fOrthDownConversionPhase = 0;			//正交下变频中的相位值
	int nDownSampClock;								//降采样时的时钟控制  
	float fDownConversionPhase;						//下变频中的相位值

	float fPLLNCO;									//锁相环中的本地NCO
	float fPLLPastFreqPart;							//锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	Complex* cFilterBuff;							//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;						//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;							//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;
	int difftemp = 0;
	float Judgth1 = 0;
	Algorithm_Demodulation De_QAM;
	int m_index;

	// add on 20240621
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int QAM64_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize
};
