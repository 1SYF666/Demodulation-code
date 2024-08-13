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

class Demo_16QAM :public QObject, public DemoduletionFactroy
{
public:

    Q_OBJECT
        // DemoduletionFactroy interface
public:
    bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;
    void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
    Demo_16QAM(int index);
	
	~Demo_16QAM();
	
    void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);
	
	float mean_function(float* data, int start, int end);

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	void Judgment16qam(Complex buffin[], int bufflen, int buffout[]);

	void diff_code_16qam(int* y, int* c_symbol, int len);

	void InitBlockCFilter();

	void QAM16Init();

	void InitialDemodulation(int Rb);

	void InitSymbolSync();

	void InitCostasPLL();

public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
    SAlgDemInit m_sAlgDemInit;                  //算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

    Algorithm_Demodulation De_QAM;
	int m_SamleSize;
	Complex* cFilterBuff;  
	float fAGCPastVc = 0.5;        
    float fPLLNCO;                              //锁相环中的本地NCO
    float fPLLPastFreqPart;                     //锁相环中的频率跟踪曲线
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	int m_Rb;
    int nDownSampClock;                         //降采样时的时钟控制
    float fDownConversionPhase;                 //下变频中的相位值
	int nPLLBuffSize;
    int m_index;

	//add on 20240624
	int difftemp = 0;

};




