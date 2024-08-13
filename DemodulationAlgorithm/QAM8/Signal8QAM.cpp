#include"Signal8QAM.h"

Demo_8QAM::Demo_8QAM(int index) :m_index(index)
{
}

Demo_8QAM::~Demo_8QAM()
{
	DELETE_ARR(Databuff0);
}

bool Demo_8QAM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb);
	return 1;
}

void Demo_8QAM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_8QAM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_8QAM_Out)
{

	Complex* DataBuff = new Complex[QAM8_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * QAM8_Samplesize);

	De_QAM.BlockFilter(data_input_slice, QAM8_Samplesize, DataBuff, cFilterBuff);

	//De_QAM.AGC(DataBuff, fAGCPastVc);

	AGCQAM8(DataBuff, fAGCPastVc, QAM8_Samplesize);

	//FILE* fp2;
	//char file_path2[500];

	/*sprintf(file_path2, "D:/data/inputI.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < QAM8_Samplesize; i++) {
		fprintf(fp2, "%f\n", data_input_slice[i].IData);
	}
	fclose(fp2);

	sprintf(file_path2, "D:/data/inputQ.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < QAM8_Samplesize; i++) {
		fprintf(fp2, "%f\n", data_input_slice[i].QData);
	}
	fclose(fp2);*/

	int nSymbolSyncSize = 0;										//突发在当前片中的长度

	Complex* DataSymbolSyncBuff = new Complex[QAM8_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * QAM8_Samplesize);

	De_QAM.SymbolSync(DataBuff, DataSymbolSyncBuff, QAM8_Samplesize, &nSymbolSyncSize, &SySyncBuffer);

	/*sprintf(file_path2, "D:/data/inputI1.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < nSymbolSyncSize; i++) {
		fprintf(fp2, "%f\n", DataSymbolSyncBuff[i].IData);
	}
	fclose(fp2);

	sprintf(file_path2, "D:/data/inputQ1.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < nSymbolSyncSize; i++) {
		fprintf(fp2, "%f\n", DataSymbolSyncBuff[i].QData);
	}
	fclose(fp2);*/

	Complex* DataPLLBuff = new Complex[nSymbolSyncSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);
	Complex* DataFLLBuff = new Complex[nSymbolSyncSize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);

	De_QAM.CarrierSync(DataSymbolSyncBuff, DataPLLBuff, DataFLLBuff, nSymbolSyncSize, 6, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	float* Signal_sort = new float[nSymbolSyncSize];
	memset(Signal_sort, 0x00, sizeof(float) * nSymbolSyncSize);

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		Signal_sort[i] = fabs(DataPLLBuff[i].IData);
	}

	De_QAM.Quick_Sort(Signal_sort, 0, nSymbolSyncSize - 1);

	//float mean_power = mean_function(Signal_sort, nSymbolSyncSize * 8 / 11, nSymbolSyncSize * 9 / 11);		   // 16QAM
	float mean_power = mean_function(Signal_sort, nSymbolSyncSize * 9.5 / 11, nSymbolSyncSize * 10.5 / 11);		   // 8QAM

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataPLLBuff[i].IData = DataPLLBuff[i].IData / mean_power;
		DataPLLBuff[i].QData = DataPLLBuff[i].QData / mean_power;
	}

	Complex* DataSymbol_m = new Complex[nSymbolSyncSize];
	memset(DataSymbol_m, 0x00, sizeof(Complex) * nSymbolSyncSize);

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataSymbol_m[i].IData = DataPLLBuff[i].IData;
		DataSymbol_m[i].QData = DataPLLBuff[i].QData;
	}


	Complex* Data_phase = new Complex[nSymbolSyncSize];
	memset(Data_phase, 0x00, sizeof(Complex) * nSymbolSyncSize);

	Phase_correction(DataSymbol_m, Data_phase, nSymbolSyncSize);

	signal_8QAM_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_8QAM_Out->burstData->softDistinguishData, Data_phase, (sizeof(Complex)) * nSymbolSyncSize);
	signal_8QAM_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	char* c_symbol = new char[nSymbolSyncSize];
	memset(c_symbol, 0x00, sizeof(char) * nSymbolSyncSize);

	//Judgment8qam(Data_phase, nSymbolSyncSize, c_symbol);
	Judgment(Data_phase, nSymbolSyncSize, c_symbol);  // modify on 20240624

	signal_8QAM_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_8QAM_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);

	//FILE* fp2;
	//char file_path2[500];

	//sprintf(file_path2, "D:/data/QAM8_output_symbol.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) 
	//{
	//	fprintf(fp2, "%d\n", c_symbol[i]);
	//}
	//fclose(fp2);



	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(Signal_sort);
	DELETE_ARR(DataSymbol_m);
	DELETE_ARR(Data_phase);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(c_symbol);
}

void Demo_8QAM::Judgment8qam(Complex buffin[], int bufflen, char buffout[])
{
	float th1 = 0.35;

	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData > th1)
		{
			buffout[i] = 1;
		}
		else if (buffin[i].IData < -th1)
		{
			buffout[i] = -1;
		}
		else
		{
			buffout[i] = 0;
		}
	}

}

void Demo_8QAM::InitBlockCFilter()
{
	if (m_Rb <= 2e6)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 9;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0424636059371430	,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
			0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430
		};


		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{

		CFilter.nFilterTaps = 9;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0424636059371430	,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
			0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430
		};


		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}



	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);

	CFilter.bCoefEvenSym = false;
	De_QAM.CFilter = CFilter;
}

void Demo_8QAM::InitialDemodulation(int Rb)
{
	m_Rb = Rb;

	InitBlockCFilter();

	QAM8Init();

	InitCostasPLL();

	InitSymbolSync();
}

void Demo_8QAM::QAM8Init()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_QAM.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;										//降采样时的时钟控制  
	fDownConversionPhase = 0;								//下变频中的相位值
	fPLLNCO = 0;											//锁相环中的本地NCO
	fPLLPastFreqPart = 0;									//锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;

	//add on 20240621
	flag = 0;
	slice = 20;
	QAM8_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[QAM8_Samplesize];

	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;						//符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;						//符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;					//符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;						//符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;							//符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;						//符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;						//符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;					//符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;					//符号同步time_error缓存值
}

void Demo_8QAM::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.001666498855103/2;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 1.389028703698034e-06/8;  // 16APSK 	
	//symbolsyncsactorInit.fSymbolSyncFactor1 = 6.665995420411309e-04;  
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 2.222445925916854e-07;  // 16QAM
	//symbolsyncsactorInit.fSymbolSyncFactor1 = 2.333098397143958e-04;
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 2.722496259248146e-08;  // 32QAM	

	//symbolsyncsactorInit.fSymbolSyncFactor1 = 6.665995420411310e-04;
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 2.222445925916854e-07;	// 8QAM

	De_QAM.m_sAlgDemInit.nSampPerSymb = 4;							//每个码元中的采样点数
	m_sAlgDemInit.nSampPerSymb = 4;										//每个码元中的采样点数

	De_QAM.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_8QAM::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;												//环路阻尼系数
	int Ko = 1;															//压控振荡器增益
	int Kd = 1;															//鉴相器增益

	float fBLcoef = 0.001;

	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_QAM.CostasPll = CostasPll;
}


void Demo_8QAM::AGCQAM8(Complex* BuffIn, float target_power, int len)
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

float Demo_8QAM::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}


float Demo_8QAM::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_8QAM::Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen)
{
	float* theta_0 = new float[bufflen];
	memset(theta_0, 0x00, sizeof(float) * bufflen);

	int n1 = 0;

	for (int i = 0; i < bufflen; i++)
	{
		float power_s = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
		// (sqrt(1.0 / 3 * 1.0 / 3 * 2) + sqrt(1.0 / 3 * 1.0 / 3 + 1 * 1)) / 2.0
		// =  0.7627
		if (power_s > (sqrt(1.0 * 1.0) * 0.95))
		{
			theta_0[n1] = atan2(Buffin[i].QData, Buffin[i].IData);

			if (theta_0[n1] > 0)
			{
				theta_0[n1] = theta_0[n1];
			}
			else
			{
				theta_0[n1] = theta_0[n1] + 2 * Pi;
			}

			n1 += 1;
		}
	}

	float* theta_1 = new float[n1];
	memset(theta_1, 0x00, sizeof(float) * n1);

	for (int i = 0; i < n1; i++)
	{
		if (theta_0[i] > 0 && theta_0[i] <= Pi / 2)
		{
			theta_1[i] = theta_0[i];
		}
		else if (theta_0[i] > Pi / 2 && theta_0[i] <= Pi)
		{
			theta_1[i] = theta_0[i] - Pi / 2;
		}
		else if (theta_0[i] > Pi / 2 * 2 && theta_0[i] <= Pi / 2 * 3)
		{
			theta_1[i] = theta_0[i] - Pi / 2 * 2;
		}
		else if (theta_0[i] > Pi / 2 * 3 && theta_0[i] <= Pi / 2 * 4)
		{
			theta_1[i] = theta_0[i] - Pi / 2 * 3;
		}
	}

	float* sort_theta_1 = new float[n1];
	memset(sort_theta_1, 0, sizeof(float) * n1);
	memcpy(sort_theta_1, theta_1, sizeof(float) * n1);

	De_QAM.Quick_Sort(sort_theta_1, 0, n1 - 1);

	float real_data = 0, imag_data = 0;
	//float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 3 / 8, (n1 - 1) * 5 / 8); //16APSK
	if (mean_sort_theta_1 == 0) {
		mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 6 / 16, (n1 - 1) * 10 / 16)- Pi / 4; //16QAM
		m_temp = new float[n1];
		memcpy(m_temp, sort_theta_1, sizeof(float) * n1);
	}
	for (int i = 0; i < bufflen; i++)
	{
		real_data = Buffin[i].IData * cos(-mean_sort_theta_1) - Buffin[i].QData * sin(-mean_sort_theta_1);
		imag_data = Buffin[i].IData * sin(-mean_sort_theta_1) + Buffin[i].QData * cos(-mean_sort_theta_1);

		Buffout[i].IData = real_data;

		Buffout[i].QData = imag_data;
	}


	DELETE_ARR(theta_0);
	DELETE_ARR(theta_1);
	DELETE_ARR(sort_theta_1);

}


void Demo_8QAM::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	float theta = 0.35;
	int* y = new int[bufflen];
	memset(y, 0x00, sizeof(int) * bufflen);
	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData >= theta)
		{
			if (buffin[i].QData >= theta)
			{
				y[i] = 1;
			}
			else if (buffin[i].QData <= -theta)
			{
				y[i] = 7;
			}
			else
			{
				y[i] = 0;
			}
		}
		else if (buffin[i].IData <= -theta)
		{
			if (buffin[i].QData >= theta)
			{
				y[i] = 3;
			}
			else if (buffin[i].QData <= -theta)
			{
				y[i] = 5;
			}
			else
			{
				y[i] = 4;
			}
		}
		else
		{
			if (buffin[i].QData >= 0)
			{
				y[i] = 2;
			}
			else
			{
				y[i] = 6;
			}
		}
	}



	for (int i = 0; i < bufflen; i++)
	{
		if (i == 0)
		{
			buffout[i] = ((y[i] - difftemp) + 8) % 8;
		}
		else
		{
			if (y[i] - y[i - 1] < 0)
			{
				buffout[i] = y[i] - y[i - 1] + 8;
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
