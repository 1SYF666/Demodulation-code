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

class Demo_4FSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_4FSK(int index);
	
	~Demo_4FSK();
	
	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);
	
	float mean_function(float* data, int start, int end);

	float max_function(float* data, int start, int end);

	int max_index_function(float* data, int start, int end);

	void Judgment_4fsk(Complex buffin[], int bufflen, char buffout[], float* threshold);

	void InitBlockCFilter();

	void FSK4Init();

    void InitialDemodulation(const DemodulationInitParamater& signalinit);

	void InitSymbolSync();

	void InitCostasPLL(float* fBLcoef);

	void Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout);

	void AGCFSK4(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);

public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
	SAlgDemInit m_sAlgDemInit;			//�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

	Algorithm_Demodulation De_FSK;
	int m_SamleSize;
	Complex* cFilterBuff;  
	float fAGCPastVc = 0.5;        
	float fPLLNCO;                    // ���໷�еı���NCO
	float fPLLPastFreqPartyemp;		  	
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	float m_Rb;
	float m_fs;
    float m_fc[8];
	int nDownSampClock;              //������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;      //�±�Ƶ�е���λֵ
	int nPLLBuffSize;
	int m_index;

	// add on 20240514
	int flag ;						// �������ݱ�־λ
	int slice;							// slice ����Ƭ��		
	int FSK_Samplesize;                // ����������
	Complex* Databuff0;					// ����������slice*m_SamleSize

	//add on 20240720
	int nSilceCount = 0;							//�ڼ����źŽ����־λ





};




