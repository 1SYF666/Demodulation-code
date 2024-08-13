#include"SignalQPSK.h"

Demo_QPSK::Demo_QPSK(int index) : m_index(index)
{
}

Demo_QPSK::~Demo_QPSK()
{
	// add on 20240719
	Ippcom.LhtFirFreeFloatLp();
}

bool Demo_QPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb); return 1;
}

void Demo_QPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_QPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_QPSK_Out)
{

	//FILE*fp;
	//char file_path1[500];

	Complex* DataBuff = new Complex[m_SamleSize];
	memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);

	//De_QPSK.BlockFilter(data_input_slice, m_SamleSize, DataBuff, cFilterBuff);
	BlockFilter_QPSK(data_input_slice, m_SamleSize, DataBuff);

	//De_QPSK.AGC(DataBuff, fAGCPastVc);

	//20240719
	//Ipp32fc* DataBufftem = new Ipp32fc[m_SamleSize];
	//memset(DataBufftem, 0x00, sizeof(Ipp32fc) * m_SamleSize);
	//memcpy(DataBufftem, DataBuff, sizeof(Ipp32fc)*m_SamleSize);

	QPSKAGC((Ipp32fc*)DataBuff, fAGCPastVc, m_SamleSize);

	Complex* DataPLLBuff = new Complex[m_SamleSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
	Complex* DataFLLBuff = new Complex[m_SamleSize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_QPSK.CarrierSync((Complex*)DataBuff, DataPLLBuff, DataFLLBuff, m_SamleSize, 2, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	// add on 20240420 ·ÀÖ¹ÔØ²¨Í¬²½Êý¾Ý¹ý´óÔì³É·ûºÅÍ¬²½Êý¾ÝÊä³öÎª0 
	//float* Signal_sort = new float[m_SamleSize];
	//memset(Signal_sort, 0x00, sizeof(float) * m_SamleSize);

	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	Signal_sort[i] = fabs(DataPLLBuff[i].IData);
	//}

	//De_QPSK.Quick_Sort(Signal_sort, 0, m_SamleSize - 1);

	//float mean_power = mean_function(Signal_sort, m_SamleSize * 8 / 11, m_SamleSize * 9 / 11);

	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	DataPLLBuff[i].IData = DataPLLBuff[i].IData / mean_power;
	//	DataPLLBuff[i].QData = DataPLLBuff[i].QData / mean_power;
	//}

	int nSymbolSyncSize = 0;										//突发在当前片中的长度
	Complex* DataSymbolSyncBuff = new Complex[m_SamleSize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);
	De_QPSK.SymbolSync(DataPLLBuff, DataSymbolSyncBuff, m_SamleSize, &nSymbolSyncSize, &SySyncBuffer);

	char* c_symbol = new char[nSymbolSyncSize];
	Judgmentqpsk(DataSymbolSyncBuff, nSymbolSyncSize, c_symbol);

	//FILE* fp;
	//char file_path1[500];

	//sprintf(file_path1, "D:/data/QPSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%d\n", c_symbol[i]);
	//}
	//fclose(fp);


	// IÂ·Êä³öºÍQÂ·Êä³ö add on 20240531
	char* symbol_I = new char[nSymbolSyncSize];
	char* symbol_Q = new char[nSymbolSyncSize];
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		symbol_I[i] = DataSymbolSyncBuff[i].IData > 0 ? 1 : -1;
		symbol_Q[i] = DataSymbolSyncBuff[i].QData > 0 ? 1 : -1;
	}

	//signal_QPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	//memcpy(signal_QPSK_Out->burstData->nDemodulationByteI, symbol_I, (sizeof(char))*nSymbolSyncSize);
	//signal_QPSK_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSize];
	//memcpy(signal_QPSK_Out->burstData->nDemodulationByteQ, symbol_Q, (sizeof(char))*nSymbolSyncSize);

	/*FILE* fp;
	char file_path1[500];
	sprintf(file_path1, "D:/data/QPSK_output_symbol_I.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSyncSize; i++) {
		fprintf(fp, "%d\n", symbol_I[i]);
	}
	fclose(fp);

	sprintf(file_path1, "D:/data/QPSK_output_symbol_Q.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSyncSize; i++) {
		fprintf(fp, "%d\n", symbol_Q[i]);
	}
	fclose(fp);*/


	signal_QPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_QPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex)) * nSymbolSyncSize);
	signal_QPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	signal_QPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_QPSK_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(c_symbol);
	DELETE_ARR(symbol_I);
	DELETE_ARR(symbol_Q);
}

// modify on 20240625
void Demo_QPSK::Judgmentqpsk(Complex buffin[], int bufflen, char buffout[])
{
	float th1 = 0.0;
	int* y = new int[bufflen];
	memset(y, 0x00, sizeof(int) * bufflen);

	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData < th1 && buffin[i].QData < th1)
		{
			y[i] = 0;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData < th1)
		{
			y[i] = 1;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= th1)
		{
			y[i] = 2;
		}
		else if (buffin[i].IData < th1 && buffin[i].QData >= th1)
		{
			y[i] = 3;
		}
	}

	for (int i = 0; i < bufflen; i++)
	{
		if (i == 0)
		{
			buffout[i] = ((y[i] - difftemp) + 4) % 4;
		}
		else
		{
			if (y[i] - y[i - 1] < 0)
			{
				buffout[i] = y[i] - y[i - 1] + 4;
			}
			else
			{
				buffout[i] = y[i] - y[i - 1];
			}
		}
	}

	difftemp = y[bufflen - 1];

	DELETE_ARR(y);

}


void Demo_QPSK::InitBlockCFilter()
{

	// ¼ÙÉèCFilterÊÇÒ»¸öÒÑ¶¨ÒåµÄ½á¹¹Ìå»òÀà£¬²¢ÇÒnFilterTapsºÍfFilterCoefÊÇÆä³ÉÔ±

	CFilter.nFilterTaps = 9;

	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

	// Ê¹ÓÃÊý×é³õÊ¼»¯µÄÖµ
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
	De_QPSK.CFilter = CFilter;
}

void Demo_QPSK::InitialDemodulation(int Rb)
{
	m_Rb = Rb;

	//InitBlockCFilter();

	QPSKInit();

	InitCostasPLL();

	InitSymbolSync();
	InitBlockCFilter_QPSK();
}

void Demo_QPSK::QPSKInit()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_QPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                                 //降采样时的时钟控制
	fDownConversionPhase = 0;                                           //下变频中的相位值
	fPLLNCO = 0;                                                        //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                               //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;
	//******************·ûºÅÍ¬²½³õÊ¼»¯***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;                                    //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;                                    //符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;                                //符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;                                    //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;                                     //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;                                  //符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;                                  //符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;                             //符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;                             //符号同步time_error缓存值
}

void Demo_QPSK::InitSymbolSync()
{

	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.003332997710206 / 10;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 5.556114814792136e-06 / 100;    // QPSK

	De_QPSK.m_sAlgDemInit.nSampPerSymb = 4;                             //每个码元中的采样点数
	m_sAlgDemInit.nSampPerSymb = 4;                                     //每个码元中的采样点数

	De_QPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_QPSK::InitCostasPLL()
{
	float BL = 0;
	float Wn = 0;
	float T_nco = 0;
	float sigma = 0.707;                                                //环路阻尼系数
	int Ko = 1;                                                         //压控振荡器增益
	int Kd = 1;                                                         //鉴相器增益
	float fBLcoef = 0.00350;                                            //QPSK
	//float fBLcoef = 0.0350;                                            //QPSK
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_QPSK.CostasPll = CostasPll;
}


float Demo_QPSK::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}
void Demo_QPSK::QPSKAGC(Ipp32fc* BuffIn, float target_power, int len)
{
	float* data_pow2 = new float[len];
	memset(data_pow2, 0, sizeof(float) * len);
	Ippcom.LhtCalculateRootOfSumOfSquaresFloat2(BuffIn, len, data_pow2);  // 求复数模值
	float input_power;
	Ippcom.LhtMeanFloat(data_pow2, len, &input_power);                    //  求均值
	float gain_factor = sqrt(target_power / input_power);
	Ippcom.LhtMulCFloat(BuffIn, gain_factor, len, nullptr);

	DELETE_ARR(data_pow2);
}
void Demo_QPSK::InitBlockCFilter_QPSK()
{
	CFilter.nFilterTaps = 9;
	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];
	float tempCoef[] =
	{
		-0.042463605937143	,1.462582099439346e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
		0.500262571713977,0.212318029685715,1.462582099439346e-17,-0.042463605937143
	};//QPSK
	for (int i = 0; i < CFilter.nFilterTaps; ++i)
	{
		CFilter.fFilterCoef[i] = tempCoef[i];
	}
	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
	CFilter.bCoefEvenSym = false;
	De_QPSK.CFilter = CFilter;
	Ippcom.LhtFirInitFloatLp(CFilter.fFilterCoef, CFilter.nFilterTaps);
}
void Demo_QPSK::BlockFilter_QPSK(Complex BuffIn[], int nBuffLen, Complex BuffOut[])
{
	Ipp32f* REAL = new Ipp32f[nBuffLen];
	memset(REAL, 0, sizeof(Ipp32f) * nBuffLen);
	Ipp32f* IMAG = new Ipp32f[nBuffLen];
	memset(IMAG, 0, sizeof(Ipp32f) * nBuffLen);
	float* REALout = new float[nBuffLen];
	memset(REALout, 0, sizeof(float) * nBuffLen);
	float* IMAGout = new float[nBuffLen];
	memset(IMAGout, 0, sizeof(float) * nBuffLen);
	Ippcom.LhtCopyRealAndimag((Ipp32fc*)BuffIn, REAL, IMAG, nBuffLen);									// 取结构体中的实部和虚部
	Ippcom.LhtFirFloatLp(REAL, REALout, nBuffLen, false);												// 实部滤波				
	Ippcom.LhtFirFloatLp(IMAG, IMAGout, nBuffLen, false);												// 虚部滤波
	Ippcom.LhtCopyRealAndimag2Complex((Ipp32fc*)BuffOut, (Ipp32f*)REALout, (Ipp32f*)IMAGout, nBuffLen);  //将实部虚部存入结构体
	DELETE_ARR(REAL);
	DELETE_ARR(REALout);
	DELETE_ARR(IMAG);
	DELETE_ARR(IMAGout);
}
