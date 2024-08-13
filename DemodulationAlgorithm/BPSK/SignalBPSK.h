
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


class Demo_BPSK :public QObject, public DemoduletionFactroy {
    Q_OBJECT
    // DemoduletionFactroy interface
public:
    bool DemoduletionInit(const DemodulationInitParamater &demoParameter) override;
    void demoduletion(Complex *dataInputRealSlice, long long baseTime, DemodulationResult *demodulationResult) override;
public:

    Demo_BPSK(int index);

    ~Demo_BPSK();

    void Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out);

    void InitialDemodulation(int Rb);

    void BPSKInit();

    void InitSymbolSync();

    void InitCostasPLL();

	void InitBlockCFilter();

    void Judgment(Complex buffin[], int bufflen, char buffout[]);

	float mean_function(float* data, int start, int end);
    
    void downphaseFrequence(Complex* data_input, float phase, Complex* downout, int dataLength);

public:
    SAlgDemInit m_sAlgDemInit;                      //算法初始化结构
    SymbolSyncFactor symbolsyncsactorInit;
    CFILTERPARAM CFilter;
    SymbolSyncBuffer SySyncBuffer;                  //符号同步初始化变量声明
    CostasPLL CostasPll;
private:
    int nSilceCount = 0;                            //当前处理的信号所在第几片
    int m_SamleSize;
    int m_Rb;
    float fAGCPastVc = 0.5;                         //AGC中的VC值缓存
    double fOrthDownConversionPhase = 0;            //正交下变频中的相位值
    int nDownSampClock;                             //降采样时的时钟控制
    float fDownConversionPhase;                     //下变频中的相位值
    float fPLLNCO;                                  //锁相环中的本地NCO
    float fPLLPastFreqPart;                         //锁相环中的频率跟踪曲线
    int nPLLBuffSize;
    Complex* cFilterBuff;                           //复数滤波器缓存区(正交下变频)
    Complex* cBlockFilterBuff;                      //复数滤波器缓存区(SRC)
    Complex* cFllterBuff;                           //复数滤波器缓存区(SRC)
    FreqLLBuffer FLLBuffer;
    float Judgth1 = 0;
    Algorithm_Demodulation De_BPSK;
    int m_index;

    // add on 20240721
    Complex* DataBuff ;
    Complex* DataPLLBuff ;
    Complex* DataFLLBuff ;
    Complex* DataSymbolSyncBuff;

};
