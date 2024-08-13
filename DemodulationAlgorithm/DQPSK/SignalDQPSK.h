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



class Demo_DQPSK :public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;
public:

	Demo_DQPSK(int index);

	~Demo_DQPSK();

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_DQPSK_Out);

	void InitBlockCFilter();

	void InitialDemodulation(int Rb);

	void DQPSKInit();

	void InitSymbolSync();

	void InitCostasPLL();

	void Judgment(Complex buffin[],int bufflen, char buffout[]);

	void AGCDQPSK(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);


public:
	SAlgDemInit m_sAlgDemInit;						//�㷨��ʼ���ṹ
	SymbolSyncFactor symbolsyncsactorInit;
	CFILTERPARAM CFilter;
	SymbolSyncBuffer SySyncBuffer;
	CostasPLL CostasPll;
private:
	int nSilceCount = 0;							//��ǰ������ź����ڵڼ�Ƭ
	int m_SamleSize;
	int m_Rb;

	float fAGCPastVc =0.5;							//AGC�е�VCֵ����

	double fOrthDownConversionPhase = 0;			//�����±�Ƶ�е���λֵ
	int nDownSampClock;								//������ʱ��ʱ�ӿ���  
	float fDownConversionPhase;						//�±�Ƶ�е���λֵ
	float fPLLNCO;									//���໷�еı���NCO
	float fPLLPastFreqPart;        					//���໷�е�Ƶ�ʸ�������
	int nPLLBuffSize;
	Complex* cFilterBuff;							//�����˲���������(�����±�Ƶ)
	Complex* cBlockFilterBuff;						//�����˲���������(SRC)
	Complex* cFllterBuff;							//�����˲���������(SRC)
	FreqLLBuffer FLLBuffer;
	float m_JudgI;
	float m_JudgQ;									//�о�����ֵ
	Algorithm_Demodulation De_DQPSK;
	int m_index;

	// add on 20240607
	int flag;								// �������ݱ�־λ
	int slice;								// slice ����Ƭ��		
	int FQPSK_Samplesize;						// ����������
	Complex* Databuff0;						// ����������slice*m_SamleSize
}; 