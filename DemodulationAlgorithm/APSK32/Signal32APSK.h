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

class Demo_32APSK: public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;
	
	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_32APSK(int index);

	~Demo_32APSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_32APSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void APSKInit();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[], float power[], int bufflen, int buffout[]);

	float mean_function(float* data, int start, int end);

	void diff_code_32apsk(int* y, char* c_symbol,int len);

	void Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen);

	void GetThreshold(Complex Buffin[], int bufflen, float Buffout[]);

	float max_function(float* data, int start, int end);

	void diff_fuction(float* input, float* output, int len);
public:
	SAlgDemInit m_sAlgDemInit;					//�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 0;						//��ǰ������ź����ڵڼ�Ƭ
	int m_SamleSize;
	int m_Rb;
	float fAGCPastVc =0.5;						//AGC�е�VCֵ����
	double fOrthDownConversionPhase = 0;		//�����±�Ƶ�е���λֵ
	int nDownSampClock;							//������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;					//�±�Ƶ�е���λֵ
	int Diffnum;								//��ֽ��뻺��ֵ
	float fPLLNCO;								//���໷�еı���NCO
	float fPLLPastFreqPart;						//���໷�е�Ƶ�ʸ�������
	int nPLLBuffSize;
	Complex* cFilterBuff;						//�����˲���������(�����±�Ƶ)
	Complex* cBlockFilterBuff;					//�����˲���������(SRC)
	Complex* cFllterBuff;						//�����˲���������(SRC)
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	float Judgth2 = 0;
	Algorithm_Demodulation De_APSK;
	int difftemp = 0;

	int m_index;

	// add on 20240514
	int flag;						// �������ݱ�־λ
	int slice;							// slice ����Ƭ��		
	int APSK_Samplesize;                // ����������
	Complex* Databuff0;					// ����������slice*m_SamleSize
};