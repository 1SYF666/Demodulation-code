#pragma once
#pragma once
#include "distinguish.h"
#include <Windows.h>
#include <iostream>
#include <stdlib.h>
#include "math.h"
#include "DemodulationModule/libHead4Demodulation.h"
#include "DemodulationModule/Interface/computing_interface.h"
#include "omp.h"
#include "ipp.h"
#pragma comment(lib,"ippsmt.lib")
#pragma comment(lib,"ippcoremt.lib")
#pragma comment(lib,"ippvmmt.lib")



class  Algorithm_Demodulation
{
public:
	SAlgDemInit m_sAlgDemInit;
	CostasPLL CostasPll = { 0,0 };
	CostasPLL CostasPll2 = { 0,0 };
	CostasPLL CostasPll3 = { 0,0 };
	FreqLL Fll;
	CFILTERPARAM CFllFilter;
	CFILTERPARAM CFilter;
	CBLOCKFILTERPARAM CBlockFilter;
	FILTERPARAM SFilter;
	SymbolSyncFactor symbolsyncsactorInit;
	SRCPARAM SRCParam; //2024/04/09 add

public:
	float meanArray(float* data, int len);

	void AGC(Complex* BuffIn, float target_power, int len = 8192);

	bool Fft(Complex m_SlipFFTBuff[], int FFT_size, float* out_FFT_data);

	void exechangNum(Complex* Inputexechang, int pos, Complex* Outputexechang, int FFT_size);

	int CountOne(unsigned char Buffin[], int nBuffSize);

	void down_orth_conversion(float Buffin[], Complex BuffOut[], double& fOrthDownConversionPhase);

	void down_conversion(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float fDownConversionFc, float* fDownConversionPhase);

	void Quick_Sort(float BufferIn[], int nQuickSortStart, int nQuickSortEnd);

	void CarrierSync(Complex BuffIn[], Complex BuffOutPLL[], Complex BuffOutFLL[], int nBuffSize, int NSignalType, float* fPLLNCO, float* fPLLPastFreqPart, FreqLLBuffer* FLLBuffer);

	void CarrierSync_32_or_64qam(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart); // nSignalType=1 32qam; nSignalType=2 64qam

	void CarrierSync_FSK(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart, float* fPLLPastFreqParttemp); //nSignalType=2--4FSK

	void PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart);

	void FLL(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, FreqLLBuffer* FLLBuffer);

	void SymbolSync(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int* nBuffSizeOut, SymbolSyncBuffer* SySyncBuffer);

	int  FindPlace(float BuffIn[], int nBuffSize, float FindNum);

	void Differential_Decoding(BYTE BuffIn[], BYTE Buffout[], int nSamleSize, int& BCache_Value);
	void normolize(Complex BuffIn[], int nBuffInLen, Complex BuffOut[], int nSignalType, NormolizeBuffer* m_NormolizeBuffer);
	//SRC
	void BlockFilter(Complex BuffIn[], int nBuffLen, Complex BuffOut[], Complex* cFilterBuff);
	void SingleBlockFilter(float BuffIn[], int nBuffLen, float BuffOut[], float* fFilterBuff);
	void FLLFilter(Complex BuffIn[], int nBuffLen, Complex BuffOut[], Complex* cFilterBuff);

	// 2024/3/5 add
	float mean_function(float* data, int start, int end);

	//2024/4/9 add for MSK
	void normolize(Complex BuffIn[], int nBuffInLen, Complex BuffOut[]);

	//2024/4/9 add for UQPSK
	void SRC(Complex cSrcBuffIn[], int nBuffInLen, Complex cSrcBuffOut[], int* nBuffOutLen, int* nDownSampClock, Complex* cFilterBuff);
	bool UpSampleFilter(Complex BuffIn[], int nBuffLen, Complex BuffOut[], int nUpSampRate, Complex* cFilterBuff);
	void DownSampleFilter(Complex BuffIn[], int nBuffInLen, Complex BuffOut[], int nDownSampRate, int nUpSampRate, int* nBuffOutLen, int* nDownSampClock, Complex* cFilterBuff);



	//--------------------------------需要继承---------------------------------------------------
public:
	double fOrthDownConversionPhase;  //0    //正交下变频中的相位值

	float fDownConversionPhase;             //下变频中的相位值
	float fDownConversionFc;                //下变频中下的频偏值(频偏估计的频偏值)
	float fPLLNCO;                   //0    //锁相环中的本地NCO
	float fPLLPastFreqPart;          //0    //锁相环中的频率跟踪曲线
	Complex* cFilterBuff;                   //复数滤波器缓存区(正交下变频)
	int   m_index = 0;                      //当前次数
	float* fFilterBuff;                     //实数滤波器缓存区
	BYTE BCache_Value;						//差分解码缓存值
	FreqLLBuffer FLLBuffer;
	SymbolSyncBuffer SySyncBuffer;
	//归一寄存
	float NormolizeMax;				//0		//缓存用于下一片归一
	float NormolizeMEAN;			//0		//缓存用于下一片归一
	int NormolizeFlag;				//0		//标识是否使用寄存值

	// add on 20240720
	ComputingInterface* interface = nullptr;


};
