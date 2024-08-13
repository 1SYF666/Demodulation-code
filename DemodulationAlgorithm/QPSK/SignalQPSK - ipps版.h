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
#include "DemodulationModule/Interface/computing_interface.h"

using namespace std;
class Demo_QPSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;
public:
	Demo_QPSK(int index);

	~Demo_QPSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_QPSK_Out);

	void Judgmentqpsk(Complex buffin[], int bufflen, char buffout[]);

	void InitBlockCFilter();

	void QPSKInit();

	void InitialDemodulation(int Rb);

	void InitSymbolSync();

	void InitCostasPLL();

	float mean_function(float* data, int start, int end);

	void QPSKAGC(Ipp32fc* BuffIn, float target_power, int len);
	void BlockFilter_QPSK(Complex BuffIn[], int nBuffLen, Complex BuffOut[]);
	void InitBlockCFilter_QPSK();
public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
	SAlgDemInit m_sAlgDemInit;                     
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

	Algorithm_Demodulation De_QPSK;
	ComputingInterface Ippcom;
	int m_SamleSize;
	Complex* cFilterBuff;
	float fAGCPastVc = 0.5;
	float fPLLNCO;                                  
	float fPLLPastFreqPart;                        
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	int m_Rb;
	int nDownSampClock;                             
	float fDownConversionPhase;                     
	int nPLLBuffSize;
	int m_index;

	// add on 20240625
	int difftemp = 0;
};

