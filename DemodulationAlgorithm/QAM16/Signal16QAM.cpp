#define _CRT_SECURE_NO_WARNINGS 1
#include"Signal16QAM.h"

Demo_16QAM::Demo_16QAM(int index) :m_index(index)
{
}

Demo_16QAM::~Demo_16QAM()
{
}

bool Demo_16QAM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb); return 1;
}
void Demo_16QAM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_16QAM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16QAM_Out)
{

	//FILE*fp;
	//char file_path1[500];

	Complex* DataBuff = new Complex[m_SamleSize];

	memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_QAM.BlockFilter(data_input_slice, m_SamleSize, DataBuff, cFilterBuff);

	De_QAM.AGC(DataBuff, fAGCPastVc);


	int nSymbolSyncSize = 0;										//突发在当前片中的长度
	Complex* DataSymbolSyncBuff = new Complex[m_SamleSize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_QAM.SymbolSync(DataBuff, DataSymbolSyncBuff, m_SamleSize, &nSymbolSyncSize, &SySyncBuffer);


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

	float mean_power = mean_function(Signal_sort, nSymbolSyncSize * 8 / 11, nSymbolSyncSize * 9 / 11);		   // 16QAM

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



	// 符号同步数据传出
	signal_16QAM_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_16QAM_Out->burstData->softDistinguishData, Data_phase, (sizeof(Complex)) * nSymbolSyncSize);
	signal_16QAM_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	//int* c_symbol = new int[nSymbolSyncSize];
	//memset(c_symbol, 0x00, sizeof(int) * nSymbolSyncSize);

	//Judgment16qam(Data_phase, nSymbolSyncSize, c_symbol);


	//int* DataResult = new int[nSymbolSyncSize];
	//memset(DataResult, 0x00, sizeof(int) * nSymbolSyncSize);

	//diff_code_16qam(c_symbol, DataResult, nSymbolSyncSize);
	////sprintf(file_path1, "D:/data/QAM16_output_symbol.txt");
	////fp = fopen(file_path1, "at");
	////for (int i = 0; i < nSymbolSyncSize; i++)
	////{
	////	fprintf(fp, "%d\n", DataResult[i]);
	////}
	////fclose(fp);

	char* DataResultout = new char[nSymbolSyncSize];
	//memcpy(DataResultout, (char*)DataResult, sizeof(char) * nSymbolSyncSize);
	memset(DataResultout, 0x00, sizeof(char) * nSymbolSyncSize);

	signal_16QAM_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_16QAM_Out->burstData->nDemodulationByteI, DataResultout, (sizeof(char)) * nSymbolSyncSize);
	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(Signal_sort);
	DELETE_ARR(DataSymbol_m);
	DELETE_ARR(Data_phase);
	//DELETE_ARR(c_symbol);
	//DELETE_ARR(DataResult);
	DELETE_ARR(DataResultout);

}

float Demo_16QAM::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_16QAM::Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen)
{
	float* theta_0 = new float[bufflen];
	memset(theta_0, 0x00, sizeof(float) * bufflen);

	int n1 = 0;

	for (int i = 0; i < bufflen; i++)
	{
		float power_s = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
		// (sqrt(1.0 / 3 * 1.0 / 3 * 2) + sqrt(1.0 / 3 * 1.0 / 3 + 1 * 1)) / 2.0
		// =  0.7627
		if (power_s < (sqrt(1.0 / 3 * 1.0 / 3 * 2) + sqrt(1.0 / 3 * 1.0 / 3 + 1 * 1)) / 2.0)
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
	float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 7 / 16, (n1 - 1) * 9 / 16); //16QAM
	auto temp = Pi / 4 - mean_sort_theta_1;
	for (int i = 0; i < bufflen; i++)
	{
		real_data = Buffin[i].IData * cos(temp) - Buffin[i].QData * sin(temp);

		imag_data = Buffin[i].IData * sin(temp) + Buffin[i].QData * cos(temp);

		Buffout[i].IData = real_data;

		Buffout[i].QData = imag_data;
	}
	delete[]sort_theta_1;
}


void Demo_16QAM::Judgment16qam(Complex buffin[], int bufflen, int buffout[])
{

	float R = 2.0 / 3;
	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData <= -R && buffin[i].QData <= -R)
		{
			buffout[i] = 3;
		}
		else if (buffin[i].IData <= -R && buffin[i].QData > -R && buffin[i].QData <= 0)
		{
			buffout[i] = 1;
		}
		else if (buffin[i].IData <= -R && buffin[i].QData > 0 && buffin[i].QData <= R)
		{
			buffout[i] = 10;
		}
		else if (buffin[i].IData <= -R && buffin[i].QData > R)
		{
			buffout[i] = 11;
		}
		else if (buffin[i].IData > -R && buffin[i].IData <= 0 && buffin[i].QData <= -R)
		{
			buffout[i] = 2;
		}
		else if (buffin[i].IData > -R && buffin[i].IData <= 0 && buffin[i].QData > -R && buffin[i].QData <= 0)
		{
			buffout[i] = 0;
		}
		else if (buffin[i].IData > -R && buffin[i].IData <= 0 && buffin[i].QData > 0 && buffin[i].QData <= R)
		{
			buffout[i] = 8;
		}
		else if (buffin[i].IData > -R && buffin[i].IData <= 0 && buffin[i].QData > R)
		{
			buffout[i] = 9;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= R && buffin[i].QData <= -R)
		{
			buffout[i] = 5;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= R && buffin[i].QData > -R && buffin[i].QData <= 0)
		{
			buffout[i] = 4;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= R && buffin[i].QData > 0 && buffin[i].QData <= R)
		{
			buffout[i] = 12;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= R && buffin[i].QData > R)
		{
			buffout[i] = 14;
		}
		else if (buffin[i].IData > R && buffin[i].QData <= -R)
		{
			buffout[i] = 7;
		}
		else if (buffin[i].IData > R && buffin[i].QData > -R && buffin[i].QData <= 0)
		{
			buffout[i] = 6;
		}
		else if (buffin[i].IData > R && buffin[i].QData > 0 && buffin[i].QData <= R)
		{
			buffout[i] = 13;
		}
		else if (buffin[i].IData > R && buffin[i].QData > R)
		{
			buffout[i] = 15;
		}
	}
}

void Demo_16QAM::diff_code_16qam(int* y, int* c_symbol, int len)
{
	int c = 0, i = 0, j = 0;
	int n = 2;

	int** a = new int* [len];
	for (i = 0; i < len; ++i)
	{
		a[i] = new int[4];
	}

	for (i = 0; i < len; i++)
	{
		j = 0;
		while (y[i] > 0)
		{
			c = (y[i] % n);
			a[i][j] = c;
			y[i] = y[i] / n;
			j++;
		}
		if (j == 0)
		{
			a[i][0] = 0; a[i][1] = 0; a[i][2] = 0; a[i][3] = 0;
		}
		else if (j == 1)
		{
			a[i][1] = 0; a[i][2] = 0; a[i][3] = 0;
		}
		else if (j == 2)
		{
			a[i][2] = 0; a[i][3] = 0;
		}
		else if (j == 3)
		{
			a[i][3] = 0;
		}
	}

	int* c_symbol_2 = new int[len];
	memset(c_symbol_2, 0, sizeof(int) * len);

	for (i = 0; i < len; i++)
	{
		if (a[i][3] == 0 && a[i][2] == 0)
		{
			c_symbol_2[i] = 0;
		}
		else if (a[i][3] == 1 && a[i][2] == 0)
		{
			c_symbol_2[i] = 1;
		}
		else if (a[i][3] == 1 && a[i][2] == 1)
		{
			c_symbol_2[i] = 2;
		}
		else if (a[i][3] == 0 && a[i][2] == 1)
		{
			c_symbol_2[i] = 3;
		}
	}

	int** d_bit = new int* [len];
	for (i = 0; i < len; ++i)
	{
		d_bit[i] = new int[4];
	}

	for (i = 0; i < len; ++i)
	{
		d_bit[i][3] = a[i][0];
		d_bit[i][2] = a[i][1];
	}

	int d_bit_2;

	for (i = 0; i < len; i++)
	{
		if (i == 0)
		{
			d_bit_2 = ((c_symbol_2[i] - difftemp) + 4) % 4; //mod要加阶数如，QP要+4
		}
		else
		{
			d_bit_2 = ((c_symbol_2[i] - c_symbol_2[i - 1]) + 4) % 4;
		}
		if (d_bit_2 == 0)
		{
			d_bit[i][0] = 0;
			d_bit[i][1] = 0;
		}
		else if (d_bit_2 == 1)
		{
			d_bit[i][0] = 1;
			d_bit[i][1] = 0;
		}
		else if (d_bit_2 == 2)
		{
			d_bit[i][0] = 1;
			d_bit[i][1] = 1;
		}
		else if (d_bit_2 == 3)
		{
			d_bit[i][0] = 0;
			d_bit[i][1] = 1;
		}
	}



	difftemp = c_symbol_2[len - 1];



	int* d_bit_1 = new int[len];
	memset(d_bit_1, 0, sizeof(int) * len);

	for (i = 0; i < len; i++)
	{
		d_bit_1[i] = d_bit[i][0] * 1000 + d_bit[i][1] * 100 + d_bit[i][2] * 10 + d_bit[i][3] * 1;
	}

	int remainder = 0;
	for (j = 0; j < len; j++)
	{
		i = 0; c_symbol[j] = 0;
		while (d_bit_1[j] != 0)
		{
			remainder = d_bit_1[j] % 10;
			d_bit_1[j] /= 10;
			c_symbol[j] += remainder * std::pow(2, i);
			++i;
		}
	}

	// Memory cleanup
	for (i = 0; i < len; ++i)delete[] a[i];

	delete[] a;

	delete[] c_symbol_2;

	for (i = 0; i < len; ++i)delete[] d_bit[i];

	delete[] d_bit;

	delete[] d_bit_1;
}

void Demo_16QAM::InitBlockCFilter()
{
	if (m_Rb <= 23e3)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 9;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// 使用数组初始化的值
		float tempCoef[] =
		{
			-0.042463605937143, 1.462582099439346e-17, 0.212318029685715, 0.500262571713977, 0.636954089057146,
			0.500262571713977, 0.212318029685715, 1.462582099439346e-17, -0.042463605937143
		};

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{

		CFilter.nFilterTaps = 23;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// 使用数组初始化的值
		float tempCoef[] = { 1.29978452894083e-17,0.00643082701099247,-4.38566444872506e-18,-0.0101055853029882,5.68512058168064e-18,0.0181900535453787,-8.12160083097234e-18,-0.0424434582725503,1.46188814957502e-17,0.212217291362751,0.500025212632458,0.636651874088254,0.500025212632458,0.212217291362751,1.46188814957502e-17,-0.0424434582725503,-8.12160083097234e-18,0.0181900535453787,5.68512058168064e-18,-0.0101055853029882,-4.38566444872506e-18,0.00643082701099247,1.29978452894083e-17 };

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

void Demo_16QAM::InitialDemodulation(int Rb) {

	m_Rb = Rb;

	InitBlockCFilter();

	QAM16Init();

	InitCostasPLL();

	InitSymbolSync();
}

void Demo_16QAM::QAM16Init()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_QAM.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                                     //降采样时的时钟控制
	fDownConversionPhase = 0;                                               //下变频中的相位值
	fPLLNCO = 0;                                                            //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                                   //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;
	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;                                        //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;                                        //符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;                                    //符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;                                        //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;                                         //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;                                      //符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;                                      //符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;                                 //符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;                                 //符号同步time_error缓存值
}

void Demo_16QAM::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 6.665995420411309e-04;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 2.222445925916854e-07;
	De_QAM.m_sAlgDemInit.nSampPerSymb = 4;                                  //每个码元中的采样点数
	m_sAlgDemInit.nSampPerSymb = 4;                                         //每个码元中的采样点数

	De_QAM.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_16QAM::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;                                                    //环路阻尼系数
	int Ko = 1;                                                             //压控振荡器增益
	int Kd = 1;                                                             //鉴相器增益
	//float fBLcoef = 0.0001;
	float fBLcoef = 0.001;  // modify on 20240527
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_QAM.CostasPll = CostasPll;
}


