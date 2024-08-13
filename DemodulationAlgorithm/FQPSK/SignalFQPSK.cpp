#include "SignalFQPSK.h"

Demo_FQPSK::Demo_FQPSK(int index) :m_index(index)
{

}

Demo_FQPSK::~Demo_FQPSK() {
	DELETE_ARR(Databuff0);

}

bool Demo_FQPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb);
	return 1;
}

void Demo_FQPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{

	memcpy(Databuff0 + flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;

	Demodulation(Databuff0, baseTime, demodulationResult);

	//Demodulation(dataInputRealSlice, baseTime, demodulationResult);

}

void Demo_FQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_FQPSK_Out)
{
	demodulation_flag++;   // 解调次数标志位

	//FILE*fp;
	//char file_path1[500];


	/*sprintf(file_path1, "D:/data/inputI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FQPSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", data_input_slice[i].IData);
	}
	fclose(fp);

	sprintf(file_path1, "D:/data/inputQ.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FQPSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", data_input_slice[i].QData);
	}
	fclose(fp);*/


	Complex* DataBuff = new Complex[FQPSK_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * FQPSK_Samplesize);

	De_FQPSK.BlockFilter(data_input_slice, FQPSK_Samplesize, DataBuff, cFilterBuff);

	// modify on 20240601
	//De_FQPSK.AGC(DataBuff, fAGCPastVc);
	/*sprintf(file_path1, "D:/data/BlockFilterI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FQPSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", DataBuff[i].QData);
	}
	fclose(fp);*/

	AGCFQPSK(DataBuff, fAGCPastVc, FQPSK_Samplesize);

	/*sprintf(file_path1, "D:/data/AGCFQPSKI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FQPSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", DataBuff[i].QData);
	}
	fclose(fp);*/


	Complex* DataPLLBuff = new Complex[FQPSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	Complex* DataFLLBuff = new Complex[FQPSK_Samplesize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * FQPSK_Samplesize);

	De_FQPSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, FQPSK_Samplesize, 4, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	/*sprintf(file_path1, "D:/data/CarrierSyncI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FQPSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", DataPLLBuff[i].IData);
	}
	fclose(fp);*/


	//int DelayOutLen = FQPSK_Samplesize - nDelaySign;
	//Complex* DataDelayOut = new Complex[DelayOutLen];
	//memset(DataDelayOut, 0x00, sizeof(Complex) * DelayOutLen);
	//DelayDelete(DataPLLBuff, FQPSK_Samplesize, DataDelayOut, cDelayBuff, nDelaySign);

	// modify on 20240604
	Complex* DataSymbolBuff = new Complex[FQPSK_Samplesize];
	memset(DataSymbolBuff, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	static float temp2[2] = { 0 };
	if (demodulation_flag == 1)
	{
		// 第一次数据解调
		int i = 0;
		for (i = 0; i < FQPSK_Samplesize - 2; i++)
		{
			DataSymbolBuff[i].IData = DataPLLBuff[i].IData;
			DataSymbolBuff[i + 2].QData = DataPLLBuff[i].QData;
		}

		temp2[0] = DataPLLBuff[i].QData;
		temp2[1] = DataPLLBuff[i + 1].QData;
		DataSymbolBuff[i].IData = DataPLLBuff[i].IData;
		DataSymbolBuff[i + 1].IData = DataPLLBuff[i + 1].IData;
		DataSymbolBuff[0].QData = 0;
		DataSymbolBuff[1].QData = 0;
	}
	else
	{
		// 第二次解调数据
		int i = 0;
		for (i = 0; i < FQPSK_Samplesize - 2; i++)
		{
			DataSymbolBuff[i].IData = DataPLLBuff[i].IData;
			DataSymbolBuff[i + 2].QData = DataPLLBuff[i].QData;
		}
		DataSymbolBuff[0].QData = temp2[0];
		DataSymbolBuff[1].QData = temp2[1];
		temp2[0] = DataPLLBuff[i].QData;
		temp2[1] = DataPLLBuff[i + 1].QData;
		DataSymbolBuff[i].IData = DataPLLBuff[i].IData;
		DataSymbolBuff[i + 1].IData = DataPLLBuff[i + 1].IData;
	}

	int nSymbolSyncSize = 0;
	Complex* DataSymbolSyncBuff1 = new Complex[FQPSK_Samplesize];
	memset(DataSymbolSyncBuff1, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	De_FQPSK.SymbolSync(DataSymbolBuff, DataSymbolSyncBuff1, FQPSK_Samplesize, &nSymbolSyncSize, &SySyncBuffer);


	Complex* DataSymbolSyncBuff = new Complex[nSymbolSyncSize];  // 码元同步输出
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);
	
	static float temp = 0;  // 缓存器
	int endline = 3;
	int nSymbolSynclen1 = nSymbolSyncSize - endline;

	if (demodulation_flag == 1)
	{
		// 第一次数据解调
		for (int i = 0; i < nSymbolSynclen1; i++)
		{
			DataSymbolSyncBuff[i].IData = -DataSymbolSyncBuff1[i + 3].IData;
			DataSymbolSyncBuff[i].QData = -DataSymbolSyncBuff1[i + 2].QData;
			if (i == (nSymbolSynclen1 - 1))
			{
				temp = -DataSymbolSyncBuff1[i + 3].QData;
			}
		}
	}
	else
	{
		// 第二次及之后数据解调
		for (int i = 0; i < nSymbolSyncSize; i++)
		{
			if (i < nSymbolSyncSize - 1)
			{
				DataSymbolSyncBuff[i].IData = -DataSymbolSyncBuff1[i].IData;
				DataSymbolSyncBuff[i + 1].QData = -DataSymbolSyncBuff1[i].QData;
			}
			else
			{
				DataSymbolSyncBuff[i].IData = -DataSymbolSyncBuff1[i].IData;
				DataSymbolSyncBuff[0].QData = temp;
				temp = -DataSymbolSyncBuff1[i].QData;
			}
		}
	}

	int nSymbolSynclen = 0;
	if (demodulation_flag == 1)
	{
		nSymbolSynclen = nSymbolSynclen1;
	}
	else
	{
		nSymbolSynclen = nSymbolSyncSize;
	}


	//// 码元同步输入
	//Complex* DataSymbolBuff = new Complex[FQPSK_Samplesize];
	//memset(DataSymbolBuff, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	//for (int i = 0; i < FQPSK_Samplesize; i++)
	//{
	//	DataSymbolBuff[i].IData = DataPLLBuff[i].IData;
	//	DataSymbolBuff[i].QData = DataPLLBuff[i].IData;
	//}

	//// I路码元同步
	//int nSymbolSyncSize = 0;
	//Complex* DataSymbolSyncBuffI = new Complex[FQPSK_Samplesize];
	//memset(DataSymbolSyncBuffI, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	//De_FQPSK.SymbolSync(DataSymbolBuff, DataSymbolSyncBuffI, FQPSK_Samplesize, &nSymbolSyncSize, &SySyncBuffer);

	//// Q路码元同步
	//for (int i = 0; i < FQPSK_Samplesize; i++)
	//{
	//	DataSymbolBuff[i].IData = DataPLLBuff[i].QData;
	//	DataSymbolBuff[i].QData = DataPLLBuff[i].QData;
	//}
	//int nSymbolSyncSizeQ = 0;
	//Complex* DataSymbolSyncBuffQ = new Complex[FQPSK_Samplesize];
	//memset(DataSymbolSyncBuffQ, 0x00, sizeof(Complex) * FQPSK_Samplesize);
	//De_FQPSK.SymbolSync(DataSymbolBuff, DataSymbolSyncBuffQ, FQPSK_Samplesize, &nSymbolSyncSizeQ, &SySyncBuffer);


	// modify on 20240604

	//// 码元同步输出
	//int endline = 3;
	//int nSymbolSynclen = nSymbolSyncSize - endline;
	//Complex* DataSymbolSyncBuff = new Complex[nSymbolSynclen];
	//memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * nSymbolSynclen);

	//for (int i = 0; i < nSymbolSynclen; i++)
	//{
	//	DataSymbolSyncBuff[i].IData = -DataSymbolSyncBuffI[i + 3].IData;
	//	DataSymbolSyncBuff[i].QData = -DataSymbolSyncBuffQ[i + 2].IData;
	//}

	/*sprintf(file_path1, "D:/data/SymbolSyncI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSynclen; i++)
	{
		fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	}
	fclose(fp);
*/

	//sprintf(file_path1, "D:/data/FQPSK_output_symbolI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSynclen; i++) {
	//	fprintf(fp, "%d\n", DataSymbolSyncBuff[i].IData > 0 ? 1 : -1);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/FQPSK_output_symbolQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSynclen; i++) {
	//	fprintf(fp, "%d\n", DataSymbolSyncBuff[i].QData > 0 ? 1 : -1);
	//}
	//fclose(fp);

	signal_FQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSynclen];
	memcpy(signal_FQPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSynclen);
	signal_FQPSK_Out->burstData->softDistinguishDataLen = nSymbolSynclen;

	char* DataResult = new char[nSymbolSynclen];
	memset(DataResult, 0x00, sizeof(char) * nSymbolSynclen);
	Judgment(DataSymbolSyncBuff, nSymbolSynclen, DataResult);

	/*sprintf(file_path1, "D:/data/FQPSK_output_symbol.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSynclen; i++) {
		fprintf(fp, "%d\n", DataResult[i]);
	}
	fclose(fp);*/

	signal_FQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSynclen];
	memcpy(signal_FQPSK_Out->burstData->nDemodulationByteI, DataResult, (sizeof(char))* nSymbolSynclen);
#pragma endregion 

#pragma region DELETE_ARR

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataResult);
	//DELETE_ARR(DataDelayOut);
	DELETE_ARR(DataSymbolBuff);
	//DELETE_ARR(DataSymbolSyncBuffI);
	//DELETE_ARR(DataSymbolSyncBuffQ);
	DELETE_ARR(DataSymbolSyncBuff1);
	DELETE_ARR(DataSymbolSyncBuff);


#pragma endregion

}

void Demo_FQPSK::InitBlockCFilter()
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

		// add on 20240420 

		CFilter.nFilterTaps = 17;
		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.00583374623082611	,0.0122805365181749	, 0.00968990653983874 , -0.0300173380010696,  -0.0547286102499980,
			0.0283416431463588,	0.247495724646726,	0.498059098399132,	0.609433757229813,	0.498059098399132,
			0.247495724646726,	0.0283416431463588, -0.0547286102499980, -0.0300173380010696,	0.00968990653983874,
			0.0122805365181749, -0.00583374623082611
		};

		//CFilter.nFilterTaps = 21;
		//CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		//float tempCoef[] =
		//{
		//	0.00643082701099247 ,-4.38566444872506e-18 ,-0.0101055853029882,	5.68512058168064e-18,	0.0181900535453787 ,
		//	-8.12160083097234e-18, -0.0424434582725503,	1.46188814957502e-17,	0.212217291362751,	0.500025212632458,
		//	0.636651874088254,	0.500025212632458,	0.212217291362751,	1.46188814957502e-17, -0.0424434582725503,
		//	-8.12160083097234e-18,	0.0181900535453787,	5.68512058168064e-18, -0.0101055853029882, -4.38566444872506e-18,
		//	0.00643082701099247
		//};

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i) {
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
	CFilter.bCoefEvenSym = false;
	De_FQPSK.CFilter = CFilter;
}

void Demo_FQPSK::InitialDemodulation(int Rb)
{
	m_Rb = Rb;
	InitBlockCFilter();
	FQPSKInit();
	InitCostasPLL();
	InitSymbolSync();
}

void Demo_FQPSK::FQPSKInit()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_FQPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                     //降采样时的时钟控制  
	fDownConversionPhase = 0;                               //下变频中的相位值
	fPLLNCO = 0;                                            //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                   //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;


	//add on 20240531
	flag = 0;
	slice = 3;
	FQPSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[FQPSK_Samplesize];

	// add on 20240604
	demodulation_flag = 0;

	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;                        //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.8;                        //符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.8;                    //符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;                        //符号同步NCO输出的定时分数间隔寄存器，初值设为0.


	SySyncBuffer.nSymbolSyncKK = 0;                         //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;                      //符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;                      //符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;                 //符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;                 //符号同步time_error缓存值

	cDelayBuff = new float[6];
	memset(cDelayBuff, 0x00, sizeof(float) * 6);
	nDelaySign = 6;                                         //第一片需要延迟七个采样点
	m_Judg = 0;                                             //判决缓存值

}

void Demo_FQPSK::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823/10;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05/100;
	De_FQPSK.m_sAlgDemInit.nSampPerSymb = 4;                    //每个码元中的采样点数
	De_FQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_FQPSK::InitCostasPLL() {
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;                                    //环路阻尼系数
	int Ko = 1;                                             //压控振荡器增益
	int Kd = 1;                                             //鉴相器增益
	float fBLcoef = 0.001;
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_FQPSK.CostasPll = CostasPll;

}

void Demo_FQPSK::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	for (int i = 0; i < bufflen; i++) {
		if (buffin[i].IData < 0 && buffin[i].QData < 0) buffout[i] = 2;
		else if (buffin[i].IData >= 0 && buffin[i].QData < 0) buffout[i] = 3;
		else if (buffin[i].IData >= 0 && buffin[i].QData >= 0) buffout[i] = 0;
		else if (buffin[i].IData >= 0 && buffin[i].QData < 0) buffout[i] = 1;
	}

	int temp = 0;
	for (int i = 0; i < bufflen; i++)
	{
		// modify on 20240531
		temp = buffout[i];
		if (buffout[i] - m_Judg < 0)
		{
			buffout[i] = (buffout[i] - m_Judg + 4) % 4;
		}
		else
		{
			buffout[i] = (buffout[i] - m_Judg) % 4;
		}
		m_Judg = temp;
	}
}

void Demo_FQPSK::DelayDelete(Complex buffin[], int bufflen, Complex buffout[], float DelayBuff[], int& nDelaySign)
{
	if (nDelaySign == 6)
	{
		for (int i = 0; i < bufflen - nDelaySign; i++)
		{
			buffout[i].IData = buffin[i + 6].IData;
			buffout[i].QData = buffin[i].QData;
		}
		for (int i = 0; i < 6; i++)
		{
			DelayBuff[i] = buffin[bufflen - 6 + i].QData;
		}
		nDelaySign = 0;
	}
	else
	{
		for (int i = 0; i < 6; i++) {
			buffout[i].QData = DelayBuff[i];
		}
		for (int i = 0; i < bufflen - 6; i++) {
			buffout[i].IData = buffin[i].IData;
			buffout[i + 7].QData = buffin[i].QData;
		}
		for (int i = 0; i < 6; i++) {
			buffout[bufflen - 6 + i].IData = buffin[bufflen - 6 + i].IData;
			DelayBuff[i] = buffin[bufflen - 6 + i].QData;
		}
	}

}

float Demo_FQPSK::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_FQPSK::AGCFQPSK(Complex* BuffIn, float target_power, int len)
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

float Demo_FQPSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}