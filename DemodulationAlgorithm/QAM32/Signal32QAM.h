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
    SAlgDemInit m_sAlgDemInit;              //�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:

	Algorithm_Demodulation De_QAM;
	int m_SamleSize;
	Complex* cFilterBuff;
	float fAGCPastVc = 0.5;
	float fPLLNCO;							//���໷�еı���NCO
	float fPLLPastFreqPart;					//���໷�е�Ƶ�ʸ�������
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;

    int m_Rb;								//������������PLL����
	int nDownSampClock;						//������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;             //�±�Ƶ�е���λֵ
	int nPLLBuffSize;

	int m_index;

	// add on 20240621
	int flag;						// �������ݱ�־λ
	int slice;							// slice ����Ƭ��		
	int QAM32_Samplesize;                // ����������
	Complex* Databuff0;					// ����������slice*m_SamleSize

	// add on 20240624
	int difftemp = 0;

};




