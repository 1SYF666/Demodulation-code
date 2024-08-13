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

class Demo_8FSK : public QObject, public DemoduletionFactroy
{
	Q_OBJECT
		// DemoduletionFactroy interface
public:

	bool DemoduletionInit(const DemodulationInitParamater& demoParameter) override;

	void demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult) override;

public:
	Demo_8FSK(int index);
	
	~Demo_8FSK();
	
	void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);
	
	float mean_function(float* data, int start, int end);

	float max_function(float* data, int start, int end);

	int max_index_function(float* data, int start, int end);

	void Judgment_8fsk(Complex buffin[], int bufflen, char buffout[], float* threshold);

	void InitBlockCFilter();

	void FSK8Init();

	void InitialDemodulation(const DemodulationInitParamater& signalinit);

	void InitSymbolSync();

	void InitCostasPLL();

	void Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout);

	void AGCFSK8(Complex* BuffIn, float target_power, int len);

	float meanArray(float* data, int len);

public:

	SymbolSyncBuffer SySyncBuffer;
	CFILTERPARAM CFilter;
	SAlgDemInit m_sAlgDemInit;  //算法初始化结构
	SymbolSyncFactor symbolsyncsactorInit;
	CostasPLL CostasPll;

private:
	Algorithm_Demodulation De_FSK;
	int m_SamleSize;
	Complex* cFilterBuff;  
	float fAGCPastVc = 0.5;        
	float fPLLNCO;							//锁相环中的本地NCO
	float fPLLPastFreqPartemp;
	FreqLLBuffer FLLBuffer;
	float Judgth1 = 0;
	int m_Rb;
	int m_fs;
	float m_fc[8];

	int nDownSampClock;						//降采样时的时钟控制  
	float fDownConversionPhase;             //下变频中的相位值
	int nPLLBuffSize;
	int m_index;

	// add on 20240514
	int flag;						// 缓存数据标志位
	int slice;							// slice 缓存片数		
	int FSK_Samplesize;                // 缓存总数据
	Complex* Databuff0;					// 缓存总数据slice*m_SamleSize

};


template <typename T>
vector<size_t>sort_indexes_(vector<T>& v)
{
	vector<size_t>idx(v.size());
	iota(idx.begin(), idx.end(), 0); //iota算法用于填充一个区间,相当于初始化
	sort
	(
		idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2)
		{
			return v[i1] < v[i2]; //从小到大排序
		}
	);
	return idx;
}
