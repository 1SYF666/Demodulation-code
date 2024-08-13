#include "SignalOQPSK.h"

Demo_OQPSK::Demo_OQPSK(int index):m_index(index)
{

}

Demo_OQPSK::~Demo_OQPSK() 
{

}

bool Demo_OQPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb);
    return 1;
}

void Demo_OQPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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


void Demo_OQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_OQPSK_Out)
{

	//FILE*fp;
	//char file_path1[500];


	//sprintf(file_path1, "D:/data/inputI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < OQPSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/inputQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < OQPSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp);


    Complex* DataBuff = new Complex[OQPSK_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * OQPSK_Samplesize);

    De_OQPSK.BlockFilter(data_input_slice, OQPSK_Samplesize, DataBuff, cFilterBuff);
  
    //De_OQPSK.AGC(DataBuff, fAGCPastVc);
	AGCOQPSK(DataBuff, fAGCPastVc, OQPSK_Samplesize);

	Complex* DataPLLBuff = new Complex[OQPSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * OQPSK_Samplesize);
	Complex* DataFLLBuff = new Complex[OQPSK_Samplesize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * OQPSK_Samplesize);

    De_OQPSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, OQPSK_Samplesize, 4, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

    int DelayOutLen = OQPSK_Samplesize - nDelaySign;
    Complex* DataDelayOut = new Complex[DelayOutLen];
    memset(DataDelayOut, 0x00, sizeof(Complex) * DelayOutLen);
    DelayDelete(DataPLLBuff, OQPSK_Samplesize, DataDelayOut, cDelayBuff, nDelaySign);
  
	int nSymbolSyncSize = 0;                                //突发在当前片中的长度
	Complex* DataSymbolSyncBuff = new Complex[DelayOutLen];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * DelayOutLen);

    De_OQPSK.SymbolSync(DataDelayOut, DataSymbolSyncBuff, DelayOutLen, &nSymbolSyncSize, &SySyncBuffer);

	signal_OQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_OQPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	signal_OQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
  
	//sprintf(file_path1, "D:/data/SymbolSyncI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	//}
	//fclose(fp);


	//sprintf(file_path1, "D:/data/OQPSK_output_symbolI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%d\n", DataSymbolSyncBuff[i].IData > 0 ? 1 : -1);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/OQPSK_output_symbolQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%d\n", DataSymbolSyncBuff[i].QData > 0 ? 1 : -1);
	//}
	//fclose(fp);

	char* symbol_I = new char[nSymbolSyncSize];
	char* symbol_Q = new char[nSymbolSyncSize];
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		symbol_I[i] = DataSymbolSyncBuff[i].IData > 0 ? 1 : -1;
		symbol_Q[i] = DataSymbolSyncBuff[i].QData > 0 ? 1 : -1;
	}
	signal_OQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_OQPSK_Out->burstData->nDemodulationByteI, symbol_I, (sizeof(char))*nSymbolSyncSize);
	signal_OQPSK_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSize];
	memcpy(signal_OQPSK_Out->burstData->nDemodulationByteQ, symbol_Q, (sizeof(char))*nSymbolSyncSize);

    int* DataResult = new int[nSymbolSyncSize];
    memset(DataResult, 0x00, sizeof(int) * nSymbolSyncSize);
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, DataResult);
	//sprintf(file_path1, "D:/data/OQPSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);
	//signal_OQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	//memcpy(signal_OQPSK_Out->burstData->nDemodulationByteI, DataResult, (sizeof(char)) * nSymbolSyncSize);



	//sprintf(file_path1, "D:/data/OQPSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);



	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataResult);
    DELETE_ARR(DataDelayOut);
	DELETE_ARR(symbol_I);
	DELETE_ARR(symbol_Q);
	// 
	DELETE_ARR(DataSymbolSyncBuff);
}

void Demo_OQPSK::InitBlockCFilter()
{
	if (m_Rb <= 23e3) {
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
	else {
		CFilter.nFilterTaps = 23;
		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// 使用数组初始化的值
		float tempCoef[] = { 1.29978452894083e-17,0.00643082701099247,-4.38566444872506e-18,-0.0101055853029882,5.68512058168064e-18,0.0181900535453787,-8.12160083097234e-18,-0.0424434582725503,1.46188814957502e-17,0.212217291362751,0.500025212632458,0.636651874088254,0.500025212632458,0.212217291362751,1.46188814957502e-17,-0.0424434582725503,-8.12160083097234e-18,0.0181900535453787,5.68512058168064e-18,-0.0101055853029882,-4.38566444872506e-18,0.00643082701099247,1.29978452894083e-17 };

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i) {
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    De_OQPSK.CFilter = CFilter;
}

void Demo_OQPSK::InitialDemodulation(int Rb) 
{
    m_Rb = Rb;
	InitBlockCFilter();
    OQPSKInit();
    InitCostasPLL();
    InitSymbolSync();
}

void Demo_OQPSK::OQPSKInit() {
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
    De_OQPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

    nDownSampClock = 1;                                             //降采样时的时钟控制  
    fDownConversionPhase = 0;                                       //下变频中的相位值
    fPLLNCO = 0;                                                    //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                           //锁相环中的频率跟踪曲线
    nPLLBuffSize = 0;

	//add on 20240607
	flag = 0;
	slice = 3;
	OQPSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[OQPSK_Samplesize];


    //******************符号同步初始化***********
    SySyncBuffer.CSymbolSyncBuff = new Complex[4];
    SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
    memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
    memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
    SySyncBuffer.fSymbolSyncW = 0.5;                                //符号同步环路滤波器输出寄存器，初值设为0.5
    SySyncBuffer.fSymbolSyncN = 0.9;                                //符号同步NCO寄存器，初值设为1
    SySyncBuffer.fSymbolSyncNTemp = 0.9;                            //符号同步NCO暂时的寄存器，初值设为1
    SySyncBuffer.fSymbolSyncU = 0.6;                                //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
    SySyncBuffer.nSymbolSyncKK = 0;                                 //符号同步用来表示Ti时间序号,指示u,yI,yQ
    SySyncBuffer.fSymbolSyncPastW = 0;                              //符号同步W的缓存值
    SySyncBuffer.fSymbolSyncPastN = 0;                              //符号同步N的缓存值
    SySyncBuffer.fSymbolSyncTimeError1 = 0;                         //符号同步time_error缓存值
    SySyncBuffer.fSymbolSyncTimeError2 = 0;                         //符号同步time_error缓存值
    
    //******************IQ路延迟***********
    cDelayBuff = new float[2];
    memset(cDelayBuff, 0x00, sizeof(float) * 2);
    nDelaySign = 2;                                                 //第一片需要延迟两个采样点
}


void Demo_OQPSK::InitSymbolSync() 
{
    symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823 / 10;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 1000;
    De_OQPSK.m_sAlgDemInit.nSampPerSymb = 4;                        //每个码元中的采样点数
    De_OQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_OQPSK::InitCostasPLL()
{
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;                                            //环路阻尼系数
    int Ko = 1;                                                     //压控振荡器增益
    int Kd = 1;                                                     //鉴相器增益
    //float fBLcoef = 0.0022;
	float fBLcoef = 0.001; //modify on 0729
    int K = Ko * Kd;
    BL = fBLcoef * m_Rb;
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (4 * (float)m_Rb);

    CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
    De_OQPSK.CostasPll = CostasPll;
}


// modify on 20240625
void Demo_OQPSK::Judgment(Complex buffin[], int bufflen, int buffout[])
{
	int* y = new int[bufflen];
	memset(y, 0x00, sizeof(int) * bufflen);

    for (int i = 0; i < bufflen; i++) {
        if (buffin[i].IData <0 && buffin[i].QData < 0) y[i] = 0;
        else if(buffin[i].IData >0 && buffin[i].QData < 0) y[i] = 1;
        else if (buffin[i].IData > 0 && buffin[i].QData > 0) y[i] = 2;
        else if (buffin[i].IData > 0 && buffin[i].QData < 0) y[i] = 3;
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

void Demo_OQPSK::DelayDelete(Complex buffin[], int bufflen, Complex buffout[], float DelayBuff[], int &nDelaySign)
{
    if (nDelaySign == 2) 
	{
        for (int i = 0; i < bufflen - 2; i++) 
		{
            buffout[i].IData = buffin[i + 2].IData;
            buffout[i].QData = buffin[i].QData;
        }
        DelayBuff[0]= buffin[bufflen-2].QData;
        DelayBuff[1]= buffin[bufflen - 1].QData;
        nDelaySign = 0;
    }
    else
	{
        buffout[0].QData = DelayBuff[0];
        buffout[1].QData = DelayBuff[1];
        for (int i = 0; i < bufflen - 2; i++)
		{
            buffout[i].IData = buffin[i].IData;
            buffout[i+2].QData = buffin[i].QData;
        }
        buffout[bufflen - 2].IData = buffin[bufflen - 2].IData;
        buffout[bufflen - 1].IData = buffin[bufflen - 1].IData;
        DelayBuff[0] = buffin[bufflen - 2].QData;
        DelayBuff[1] = buffin[bufflen - 1].QData;
    }

}

// add on 20240607
void Demo_OQPSK::AGCOQPSK(Complex* BuffIn, float target_power, int len)
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

float Demo_OQPSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}