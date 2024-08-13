#include "SignalDQPSK.h"

Demo_DQPSK::Demo_DQPSK(int index) :m_index(index)
{

}

Demo_DQPSK::~Demo_DQPSK()
{

}

bool Demo_DQPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb);
	return 1;
}

void Demo_DQPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	memcpy(Databuff0 + flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;

	Demodulation(Databuff0, baseTime, demodulationResult);

	// Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_DQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_DQPSK_Out)
{
	int datalength = FQPSK_Samplesize;

	//FILE*fp;
	//char file_path1[500];

	//sprintf(file_path1, "D:/data/inputI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FQPSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/inputQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FQPSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp);

	Complex* DataBuff = new Complex[datalength];
	memset(DataBuff, 0x00, sizeof(Complex) * datalength);

	De_DQPSK.BlockFilter(data_input_slice, datalength, DataBuff, cFilterBuff);

	//De_DQPSK.AGC(DataBuff, fAGCPastVc);
	AGCDQPSK(DataBuff, fAGCPastVc, datalength);

	Complex* DataPLLBuff = new Complex[datalength];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * datalength);
	Complex* DataFLLBuff = new Complex[datalength];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * datalength);

	De_DQPSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, datalength, 13, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	int nSymbolSyncSize = 0;
	Complex* DataSymbolSyncBuff = new Complex[datalength];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * datalength);

	De_DQPSK.SymbolSync(DataPLLBuff, DataSymbolSyncBuff, datalength, &nSymbolSyncSize, &SySyncBuffer);


	//sprintf(file_path1, "D:/data/DQPSK_output_symbolI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/DQPSK_output_symbolQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].QData);
	//}
	//fclose(fp);

	signal_DQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_DQPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	signal_DQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	char* DataResult = new char[nSymbolSyncSize];
	memset(DataResult, 0x00, sizeof(char) * nSymbolSyncSize);
	Judgment(DataSymbolSyncBuff, nSymbolSyncSize, DataResult);

	//sprintf(file_path1, "D:/data/DQPSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);
	signal_DQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_DQPSK_Out->burstData->nDemodulationByteI, DataResult, (sizeof(char)) * nSymbolSyncSize);

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataResult);
	DELETE_ARR(DataSymbolSyncBuff);

}

void Demo_DQPSK::InitBlockCFilter()
{
	if (m_Rb <= 23e3)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
		CFilter.nFilterTaps = 17;
		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// 使用数组初始化的值
		float tempCoef[] = { -0.0300349916229418, -0.0283172616659510, 1.03450095994087e-17, 0.0606798464270379, 0.150174958114709,
			0.254855354993559, 0.353841408874743, 0.424758924989266, 0.450524874344127, 0.424758924989266, 0.353841408874743,
			0.254855354993559, 0.150174958114709, 0.0606798464270379, 1.03450095994087e-17, -0.0283172616659510, -0.03003499162294183 };

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i) {
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{
		CFilter.nFilterTaps = 13;
		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// 使用数组初始化的值
		//float tempCoef[] =
		//{ 
		//	1.29978452894083e-17,0.00643082701099247,-4.38566444872506e-18,-0.0101055853029882,5.68512058168064e-18,
		//	0.0181900535453787,-8.12160083097234e-18,-0.0424434582725503,1.46188814957502e-17,0.212217291362751,
		//	0.500025212632458,0.636651874088254,0.500025212632458,0.212217291362751,1.46188814957502e-17,-0.0424434582725503,
		//	-8.12160083097234e-18,0.0181900535453787,5.68512058168064e-18,-0.0101055853029882,-4.38566444872506e-18,0.00643082701099247,
		//	1.29978452894083e-17
		//};

		float tempCoef[] =
		{
			0.0181926639868083,-8.12276635603421e-18,-0.0424495493025527,1.46209794408616e-17,0.212247746512764,
			0.500096970889115,0.636743239538291,0.500096970889115,0.212247746512764,
			1.46209794408616e-17,-0.0424495493025527,-8.12276635603421e-18,0.0181926639868083
		};

		//CFilter.nFilterTaps = 9;
		//CFilter.fFilterCoef = new float[CFilter.nFilterTaps];
		//float tempCoef[] =
		//{
		//	-0.0547932924058631, 0.0283751393117299, 0.247788232659709, 0.498647739990532, 
		//	0.610154029297649,0.498647739990532, 0.247788232659709, 0.0283751393117299,
		//	-0.0547932924058631
		//};


		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i) 
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
	CFilter.bCoefEvenSym = false;
	De_DQPSK.CFilter = CFilter;
}

void Demo_DQPSK::InitialDemodulation(int Rb)
{
	m_Rb = Rb;
	InitBlockCFilter();
	DQPSKInit();
	InitCostasPLL();
	InitSymbolSync();
}

void Demo_DQPSK::DQPSKInit() {
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_DQPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                     //降采样时的时钟控制  
	fDownConversionPhase = 0;                               //下变频中的相位值
	fPLLNCO = 0;                                            //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                   //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;

	// add on 20240607
	flag = 0;
	slice = 1;
	FQPSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[FQPSK_Samplesize];

	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;                        //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;                        //符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;                    //符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;                        //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;                         //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;                      //符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;                      //符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;                 //符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;                 //符号同步time_error缓存值

	m_JudgI = 0;
	m_JudgQ = 0;
}

void Demo_DQPSK::InitSymbolSync() {
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823 / 10;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 100;
	De_DQPSK.m_sAlgDemInit.nSampPerSymb = 4;                //每个码元中的采样点数
	De_DQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_DQPSK::InitCostasPLL() {
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;                                    //环路阻尼系数
	int Ko = 1;                                             //压控振荡器增益
	int Kd = 1;                                             //鉴相器增益
	float fBLcoef = 0.02;
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_DQPSK.CostasPll = CostasPll;
}

void Demo_DQPSK::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	Complex* Judgm = new Complex[bufflen];

	for (int i = 0; i < bufflen; i++)
	{
		Judgm[i].IData = buffin[i].IData * m_JudgI + buffin[i].QData * m_JudgQ;
		Judgm[i].QData = buffin[i].QData * m_JudgI - buffin[i].IData * m_JudgQ;
		m_JudgI = buffin[i].IData;
		m_JudgQ = buffin[i].QData;
		if (Judgm[i].IData > 0 && Judgm[i].QData > 0) buffout[i] = 0;
		else if (Judgm[i].IData < 0 && Judgm[i].QData > 0) buffout[i] = 1;
		else if (Judgm[i].IData < 0 && Judgm[i].QData < 0) buffout[i] = 3;
		else if (Judgm[i].IData > 0 && Judgm[i].QData < 0) buffout[i] = 2;
	}

	delete(Judgm);
}


void Demo_DQPSK::AGCDQPSK(Complex* BuffIn, float target_power, int len)
{

	float* data_pow2 = new float[len];
	memset(data_pow2, 0, sizeof(float) * len);
	for (int i = 0; i < len; i++)
	{
		//data_pow2[i] = BuffIn[i].IData * BuffIn[i].IData                           //之前
		data_pow2[i] = (BuffIn[i].IData * BuffIn[i].IData + BuffIn[i].QData * BuffIn[i].QData) / 2;    //4.26修改
	}
	float input_power = meanArray(data_pow2, len);
	float gain_factor = sqrt(target_power / input_power);
	for (int i = 0; i < len; i++)
	{
		BuffIn[i].IData = BuffIn[i].IData * gain_factor;
		BuffIn[i].QData = BuffIn[i].QData * gain_factor;
	}
	DELETE_ARR(data_pow2);
}

float Demo_DQPSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}