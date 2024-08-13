#pragma once
#include "distinguish.h"
#include "../Algorithm_Demodulation.h"
#include "../../demoduletion_factroy.h"
#include"../../Interface/computing_interface.h"
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

class Demo_CDMA : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:
	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_CDMA(int index);

	~Demo_CDMA();

	void max_index_array(float* data, int len, float* output);

	void ippfft(float* data_real, float* data_imag, int FFT_LENGTH1, int flag1, float* Ioutput, float* Qoutput, float* output);

	void ippifft(float* data_real, float* data_imag, int FFT_LENGTH1, float* ifft_real, float* ifft_imag, float* Module);

	void Capture(float* inputI, float* inputQ);

	void Tracking(float* inputI, float* inputQ, int symbol_num, Complex* demoResult);

	void Tracking_new(float* inputI, float* inputQ, int symbol_num, Complex* PNsyndata);

	void Initial_Tracking();

	void InitCostasPLL();

	void CDMADemo(float* input_I, float* input_Q, int signLen, DemodulationResult* signal_16APSK_Out);

	void InitialDemodulation(const DemodulationInitParamater& signalinit);

	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

	float mean_function(float* data, int start, int end);

	void AGC(Complex* BuffIn, float target_power);

	float meanArray(float* data, int len);

	// add on 20240622
	void AGCCDMA(Complex* BuffIn, float target_power, int len);

	void InitBlockCFilter();

	void Initial_Tracking_new();

	// add on 20240803
	void estCDMASpreadCode(float* input_I, float* input_Q, float fs,int pnEstLen,int* pncode);
	int maxIndex(float* data, int start, int end);
	int minIndex(float* data, int start, int end);
	float max_function(float* data, int start, int end);
public:
	CostasPLL CostasPll;
	CFILTERPARAM CFilter;
private:
	Complex* cFilterBuff;
	int pnLength;				// PN�볤��
	float fs;					// ������
	int* pnCode;				// PN������
	int m_index;
	SAlgDemInit m_sAlgDemInit;  // �㷨��ʼ���ṹ

	float fAGCPastVc = 0.5;		//AGC�е�VCֵ����

	// add on 20240429
	int flag = 0;               // �������ݱ�־λ
	int slice = 0;             // ��Ƶ�����Ƭ��
	int m_SamleSize = 0;     // �����źŽ��һ�ε���Ϊ8192
	int CDMA_Samplesize = 0; // ��Ƶ���һ���ܵ���
	Complex* Databuf;


	int chip_sps;				// ��Ƭ��������
	int flag_acqusiation;   // �����־λ 1--���в���0--�����в���
	int* delayDot_buf;			// �����ӳٵ�������Ƶ�Ĵ��� 0--��Ƶ��1--�ӳٵ�
	int capture_fc;			// ������Ƶֵ
	int pncodelength_sps;       // PN������󳤶�
	int fc_acqusiation;         // ���Ƶ��
	int step_acqusiation;		// ɨƵ����
	int step_fc_acqusiation;    // ɨƵֵ
	int *fc_acq_buff;           // ɨƵ�Ĵ���s
	float* PNcode_est_sps;      // PN������Ĵ���
	//float* PNfftdataI;			// PN��fft���
	//float* PNfftdataQ;
	//float* PNfftdata;
	//float* CapturedataI;       // �����ź�����
	//float* CapturedataQ;
	//float* SignalFFTI;			// �����ź�����fft	
	//float* SignalFFTQ;
	//float* Signaldata;
	//float* DotMultI;			// ���������Ƶ��
	//float* DotMultQ;
	//float* SignalifftI;
	//float* SignalifftQ;
	//float* Signaldata_temp;



	Algorithm_Demodulation De_CDMA;
	ComputingInterface Computer;    // add on 20240807

	struct TrackInit
	{
		float ThresholdPositive;
		float ThresholdNegative;
		float TrackingLoopCoef1;
		float TrackingLoopCoef2;
		float TrackingLoopCoef_new1;
		float TrackingLoopCoef_new2;
	};

	TrackInit trackcoef;

	int sum_temp1;				// ��Ƶ�źŵ�һ�λ��ֵ���
	int sum_temp2;             // ��Ƶ�źŵڶ��λ��ֵ���

	float fPLLNCO;								//���໷�еı���NCO
	float fPLLPastFreqPart;						//���໷�е�Ƶ�ʸ�������
	FreqLLBuffer FLLBuffer;

	// add on 2024/05/06
	int* quotient;						// �� 
	int* remainder;						// ���� 
	int flag_trcking = 0;
	float* data_track_temp_I;			// �������ݻ�����
	float* data_track_temp_Q;			// �������ݻ�����


	// add on 20240803

	int flag1 = 1;                     // ��Ƶ���й��Ʊ�־λ



};




