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

class Demo_32QAM :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:

	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_32QAM(int index);

	~Demo_32QAM();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	float mean_function(float* data, int start, int end);

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	void Judgment32qam(Complex buffin[], int bufflen, int buffout[]);

	void diff_code_32qam(int* y, int* c_symbol, int len);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void QAM32Init();

	void InitSymbolSync();

	void InitCostasPLL();

	void AGCQAM32(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);



public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
    SAlgDemInit m_sAlgDemInit;              //算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

	Algorithm_Demodulation De_QAM;
	int m_SamleSize;
	Complex* cFilterBuff;
	float fAGCPastVc = 0.5;
	float fPLLNCO;							//锁相环中的本地NCO
	float fPLLPastFreqPart;					//锁相环中的频率跟踪曲线
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;

    int m_Rb;								//码速用于设置PLL参数
	int nDownSampClock;						//降采样时的时钟控制  
	float fDownConversionPhase;             //下变频中的相位值
	int nPLLBuffSize;

	int m_index;

	// add on 20240621
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int QAM32_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize

	// add on 20240624
	int difftemp = 0;

};




