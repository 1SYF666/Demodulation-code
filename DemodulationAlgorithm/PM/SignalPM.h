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



class Demo_PM : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_PM(int index);

	~Demo_PM();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_AM_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void AMInit();

	void InitCostasPLL();

	float mean_function(float* data, int start, int end);

	void PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart, float* PLL_Freq_Part, float* PLL_Phase_Part);

public:
	SAlgDemInit m_sAlgDemInit;						//算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 0;							//当前处理的信号所在第几片
	int m_SamleSize;
	int m_Rb;

	float fAGCPastVc =0.5;							//AGC中的VC值缓存

	double fOrthDownConversionPhase = 0;			//正交下变频中的相位值
	int nDownSampClock;								//降采样时的时钟控制  
	float fDownConversionPhase;						//下变频中的相位值
	float fPLLNCO;									//锁相环中的本地NCO
	float fPLLPastFreqPart;							//锁相环中的频率跟踪曲线
	int nPLLBuffSize;
	int nDelaySign ;								//标记延迟是否为第一片
	float* cDelayBuff;								//延迟需要的缓存长度

	Complex* cFilterBuff;							//复数滤波器缓存区(正交下变频)
	Complex* cBlockFilterBuff;						//复数滤波器缓存区(SRC)
	Complex* cFllterBuff;							//复数滤波器缓存区(SRC)
	FreqLLBuffer FLLBuffer;
	float m_JudgI;
	float m_JudgQ;									//判决缓存值
	int m_index;
	Algorithm_Demodulation Demo;
};