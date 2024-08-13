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

class Demo_8QAM :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:

	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_8QAM(int index);
	
	~Demo_8QAM();
	
	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);
	
	void Judgment8qam(Complex buffin[], int bufflen, char buffout[]);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void QAM8Init();

	void InitSymbolSync();

	void InitCostasPLL();

	void AGCQAM8(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);

	float mean_function(float* data, int start, int end);

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	void Judgment(Complex buffin[], int bufflen, char buffout[]); //add on 20240624

	void Phase_correction2(Complex Buffin[], Complex Buffout[], int bufflen);

public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
	SAlgDemInit m_sAlgDemInit;					//�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

	Algorithm_Demodulation De_QAM;
	int m_SamleSize;
	Complex* cFilterBuff;  
	float fAGCPastVc = 0.5;        
	float fPLLNCO;								//���໷�еı���NCO
	float fPLLPastFreqPart;						//���໷�е�Ƶ�ʸ�������
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	int m_Rb;
	int nDownSampClock;							//������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;					//�±�Ƶ�е���λֵ
	int nPLLBuffSize;

	// add on 20240623
	int m_index;
	int flag_phase;
	int flag_size ;
	//float* mean_sort_theta_1;
	float mean_sort_theta_1 = 0;
	float* mean_sort_theta_one;
	float *m_temp ;
	// add on 20240621
	int flag;						// �������ݱ�־λ
	int slice;							// slice ����Ƭ��		
	int QAM8_Samplesize;                // ����������
	Complex* Databuff0;					// ����������slice*m_SamleSize

	// add on 20240624
	int difftemp = 0;

};




