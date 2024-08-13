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

class Demo_8PSK : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:

	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_8PSK(int index);

	~Demo_8PSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void PSK_8_Init();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[], int bufflen, char buffout[]);

public:
	SAlgDemInit m_sAlgDemInit;							//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
	CostasPLL CostasPll2;                               
	CostasPLL CostasPll3;                               
private:
	int nSilceCount = 0;								//当前处理的信号所在第几片
	int m_SamleSize;
	int m_Rb;
	float fAGCPastVc = 0.5;								//AGC中的VC值缓存
	double fOrthDownConversionPhase = 0;				//正交下变频中的相位值
	int nDownSampClock;									//降采样时的时钟控制  
	float fDownConversionPhase;							//下变频中的相位值
	float fPLLNCO;										//锁相环中的本地NCO
	float fPLLPastFreqPart;								//锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	Complex* cFilterBuff;								//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;							//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;								//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	int m_index;
	int difftemp = 0;
	Algorithm_Demodulation De_8PSK;

	// add on 20240621
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int PSK8_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize

};