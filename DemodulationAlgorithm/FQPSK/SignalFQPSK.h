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

class Demo_FQPSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:

	Demo_FQPSK(int index);

	~Demo_FQPSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_FQPSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void FQPSKInit();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[], int bufflen, char buffout[]);

	void DelayDelete(Complex buffin[], int bufflen, Complex buffout[], float DelayBuff[], int& nDelaySign);

	float mean_function(float* data, int start, int end);

	// modify on 20240601
	void AGCFQPSK(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);

public:
	SAlgDemInit m_sAlgDemInit;							//�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 0;								//��ǰ������ź����ڵڼ�Ƭ
	int m_SamleSize;
	int m_Rb;
	float fAGCPastVc = 0.5;								//AGC�е�VCֵ����
	double fOrthDownConversionPhase = 0;				//�����±�Ƶ�е���λֵ
	int nDownSampClock;									//������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;							//�±�Ƶ�е���λֵ
	float fPLLNCO;										//���໷�еı���NCO
	float fPLLPastFreqPart;								//���໷�е�Ƶ�ʸ�������
	int nPLLBuffSize;

	int nDelaySign;										//����ӳ��Ƿ�Ϊ��һƬ
	float* cDelayBuff;									//�ӳ���Ҫ�Ļ��泤��

	Complex* cFilterBuff;								//�����˲���������(�����±�Ƶ)
	Complex* cBlockFilterBuff;							//�����˲���������(SRC)
	Complex* cFllterBuff;								//�����˲���������(SRC)
	FreqLLBuffer FLLBuffer;
	int m_Judg;
	Algorithm_Demodulation De_FQPSK;

	int m_index;


	// add on 20240531
	int flag;								// �������ݱ�־λ
	int slice;								// slice ����Ƭ��		
	int FQPSK_Samplesize;						// ����������
	Complex* Databuff0;						// ����������slice*m_SamleSize

	// add on 20240604
	int demodulation_flag;

};