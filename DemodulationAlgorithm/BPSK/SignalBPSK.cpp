#include "SignalBPSK.h"

bool Demo_BPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb); return 1;
}

void Demo_BPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

Demo_BPSK::Demo_BPSK(int index) : m_index(index) {

}

Demo_BPSK::~Demo_BPSK() {

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataSymbolSyncBuff);

}

void Demo_BPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_BPSK_Out) {


	//Complex* DataBuff = new Complex[m_SamleSize];
	//memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_BPSK.BlockFilter(data_input_slice, m_SamleSize, DataBuff, cFilterBuff);

	De_BPSK.AGC(DataBuff, fAGCPastVc);


	//Complex* DataPLLBuff = new Complex[m_SamleSize];
	//memset(DataPLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
	//Complex* DataFLLBuff = new Complex[m_SamleSize];
	//memset(DataFLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
	De_BPSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, m_SamleSize, 1, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);


	int nSymbolSyncSize = 0;
	//Complex* DataSymbolSyncBuff = new Complex[m_SamleSize];
	//memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);
	De_BPSK.SymbolSync(DataPLLBuff, DataSymbolSyncBuff, m_SamleSize, &nSymbolSyncSize, &SySyncBuffer);


	char* c_symbol = new char[nSymbolSyncSize];
	Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbol);

	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].IData;
	//}
	////调用星座图绘制接口
	//signal_BPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	//memcpy(signal_BPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex)) * nSymbolSyncSize);
	//signal_BPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;


	int starlen = nSymbolSyncSize;
	Complex* StarBuff = new Complex[starlen];
	memset(StarBuff, 0x00, sizeof(Complex) * starlen);

	for (int i = 0; i < starlen; i++)
	{
		StarBuff[i].IData = DataSymbolSyncBuff[i].IData;
		StarBuff[i].QData = DataSymbolSyncBuff[i].QData;
	}
	// 移动相位pi/4
	Complex* StarBuffout = new Complex[starlen];
	memset(StarBuffout, 0x00, sizeof(Complex) * starlen);
	float phase = Pi / 4;
	downphaseFrequence(StarBuff, phase, StarBuffout, starlen);

	signal_BPSK_Out->burstData->softDistinguishData = new Complex[starlen];
	memcpy(signal_BPSK_Out->burstData->softDistinguishData, StarBuffout, (sizeof(Complex)) * starlen);
	signal_BPSK_Out->burstData->softDistinguishDataLen = starlen;





	signal_BPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_BPSK_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);

	signal_BPSK_Out->burstData->nBW = m_Rb;

	//DELETE_ARR(DataBuff);
	//DELETE_ARR(DataPLLBuff);
	//DELETE_ARR(DataFLLBuff);
	//DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(c_symbol);
}

void Demo_BPSK::InitialDemodulation(int Rb)
{

	m_Rb = Rb;
	BPSKInit();
	InitBlockCFilter();
	InitCostasPLL();
	InitSymbolSync();
	// add on 20240721
	DataBuff = new Complex[m_SamleSize];
	memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);

	DataPLLBuff = new Complex[m_SamleSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
	DataFLLBuff = new Complex[m_SamleSize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * m_SamleSize);

	DataSymbolSyncBuff = new Complex[m_SamleSize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);


}

void Demo_BPSK::BPSKInit()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_BPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                 //降采样时的时钟控制
	fDownConversionPhase = 0;                           //下变频中的相位值
	fPLLNCO = 0;                                        //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                               //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;
	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;                    //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;                    //符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;                //符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;                    //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;                     //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;                  //符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;                  //符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;             //符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;             //符号同步time_error缓存值
}


void Demo_BPSK::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 9.137764697780231e-04;            //更改BL，0.05
	symbolsyncsactorInit.fSymbolSyncFactor2 = 4.176198395515377e-07;        //更改BL，0.05
	De_BPSK.m_sAlgDemInit.nSampPerSymb = 4;                                    //每个码元中的采样点数
	De_BPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_BPSK::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;     //环路阻尼系数
	int Ko = 1;              //压控振荡器增益
	int Kd = 1;              //鉴相器增益
	float fBLcoef = 0.002;
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb * 4;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_BPSK.CostasPll = CostasPll;
}

void Demo_BPSK::InitBlockCFilter()
{

	// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

	CFilter.nFilterTaps = 9;

	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

	// 使用数组初始化的值
	float tempCoef[] =
	{
		-0.042463605937143	,1.462582099439346e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
		0.500262571713977,0.212318029685715,1.462582099439346e-17,-0.042463605937143
	};//QPSK

	// 将初始化值复制到动态分配的数组中
	for (int i = 0; i < CFilter.nFilterTaps; ++i)
	{
		CFilter.fFilterCoef[i] = tempCoef[i];
	}

	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);

	CFilter.bCoefEvenSym = false;
	De_BPSK.CFilter = CFilter;
}



void Demo_BPSK::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData > 0)
		{
			buffout[i] = 1;
		}
		else
		{
			buffout[i] = -1;
		}
	}
}

float Demo_BPSK::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_BPSK::downphaseFrequence(Complex* data_input, float phase, Complex* downout, int dataLength)
{
	//两路
	for (int i = 0; i < dataLength; i++)
	{
		downout[i].IData = data_input[i].IData * cos(phase) - data_input[i].QData * sin(phase);
		downout[i].QData = data_input[i].QData * cos(phase) + data_input[i].IData * sin(phase);
	}
}