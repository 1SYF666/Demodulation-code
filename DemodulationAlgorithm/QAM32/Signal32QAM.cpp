#include"Signal32QAM.h"

Demo_32QAM::Demo_32QAM(int index) :m_index(index)
{
}

Demo_32QAM::~Demo_32QAM()
{
}

bool Demo_32QAM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb);
	return 1;
}

void Demo_32QAM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	//Demodulation(dataInputRealSlice, baseTime, demodulationResult);
	memcpy(Databuff0 + flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;

	Demodulation(Databuff0, baseTime, demodulationResult);
	//DELETE_ARR(Databuff0);
}

void Demo_32QAM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_32QAM_Out)
{

	//FILE* fp;
	//char file_path1[500];

	Complex* DataBuff = new Complex[QAM32_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * QAM32_Samplesize);

	De_QAM.BlockFilter(data_input_slice, QAM32_Samplesize, DataBuff, cFilterBuff);

	//De_QAM.AGC(DataBuff, fAGCPastVc);
	AGCQAM32(DataBuff, fAGCPastVc, QAM32_Samplesize);

	int nSymbolSyncSize = 0;
	Complex* DataSymbolSyncBuff = new Complex[QAM32_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * QAM32_Samplesize);

	De_QAM.SymbolSync(DataBuff, DataSymbolSyncBuff, QAM32_Samplesize, &nSymbolSyncSize, &SySyncBuffer);

	Complex* DataPLLBuff = new Complex[nSymbolSyncSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);
	Complex* DataFLLBuff = new Complex[nSymbolSyncSize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);

	De_QAM.CarrierSync_32_or_64qam(DataSymbolSyncBuff, DataPLLBuff, nSymbolSyncSize, 1, &fPLLNCO, &fPLLPastFreqPart);

	float* Signal_sort = new float[nSymbolSyncSize];
	memset(Signal_sort, 0x00, sizeof(float) * nSymbolSyncSize);

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		Signal_sort[i] = fabs(DataPLLBuff[i].IData);
	}

	De_QAM.Quick_Sort(Signal_sort, 0, nSymbolSyncSize - 1);

	float mean_power = mean_function(Signal_sort, nSymbolSyncSize * 9 / 11, nSymbolSyncSize * 9.5 / 11);		   // 32QAM

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

	Phase_correction(DataSymbol_m, Data_phase, nSymbolSyncSize);

	signal_32QAM_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_32QAM_Out->burstData->softDistinguishData, Data_phase, (sizeof(Complex)) * nSymbolSyncSize);
	signal_32QAM_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	int* c_symbol = new int[nSymbolSyncSize];
	Judgment32qam(Data_phase, nSymbolSyncSize, c_symbol);

	int* DataResult = new int[nSymbolSyncSize];
	diff_code_32qam(c_symbol, DataResult, nSymbolSyncSize);

	//sprintf(file_path1, "D:/data/QAM32_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);

	char* DataResultout = new char[nSymbolSyncSize];
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataResultout[i] = DataResult[i];
	}


	signal_32QAM_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_32QAM_Out->burstData->nDemodulationByteI, DataResultout, (sizeof(char)) * nSymbolSyncSize);

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(Signal_sort);
	DELETE_ARR(DataSymbol_m);
	DELETE_ARR(Data_phase);
	DELETE_ARR(c_symbol);
	DELETE_ARR(DataResult);
	DELETE_ARR(DataResultout);

}

float Demo_32QAM::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_32QAM::Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen)
{
	float* theta_0 = new float[bufflen];
	memset(theta_0, 0x00, sizeof(float) * bufflen);

	int n1 = 0;

	for (int i = 0; i < bufflen; i++)
	{
		float power_s = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
		// (sqrt(1.0 / 3 * 1.0 / 3 * 2) + sqrt(1.0 / 3 * 1.0 / 3 + 1 * 1)) / 2.0
		// =  0.7627
		if (power_s < (sqrt(1.0 / 5 * 1.0 / 5 * 2) + sqrt(1.0 / 5 * 1.0 / 5 + 3.0 / 5 * 3.0 / 5)) / 2.0)
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
	//float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 7 / 16, (n1 - 1) * 9 / 16); //16QAM
	float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 3 / 8, (n1 - 1) * 5 / 8); //32QAM
	for (int i = 0; i < bufflen; i++)
	{
		real_data = Buffin[i].IData * cos(Pi / 4 - mean_sort_theta_1) - Buffin[i].QData * sin(Pi / 4 - mean_sort_theta_1);

		imag_data = Buffin[i].IData * sin(Pi / 4 - mean_sort_theta_1) + Buffin[i].QData * cos(Pi / 4 - mean_sort_theta_1);

		Buffout[i].IData = real_data;

		Buffout[i].QData = imag_data;
	}

	DELETE_ARR(theta_0);
	DELETE_ARR(theta_1);
	DELETE_ARR(sort_theta_1);

}

void Demo_32QAM::Judgment32qam(Complex buffin[], int bufflen, int buffout[])
{
	float th1 = 2.0 / 5;
	float th2 = 4.0 / 5;

	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData <= -th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 24;
		}
		else if (buffin[i].IData <= -th2 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 25;
		}
		else if (buffin[i].IData <= -th2 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 19;
		}
		else if (buffin[i].IData <= -th2 && buffin[i].QData > th1)
		{
			buffout[i] = 18;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData <= -th2)
		{
			buffout[i] = 26;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData > -th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 28;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 29;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 22;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData > th1 && buffin[i].QData <= th2)
		{
			buffout[i] = 20;
		}
		else if (buffin[i].IData > -th2 && buffin[i].IData <= -th1 && buffin[i].QData > th2)
		{
			buffout[i] = 16;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData <= -th2)
		{
			buffout[i] = 27;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData > -th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 30;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 31;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 23;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData > th1 && buffin[i].QData <= th2)
		{
			buffout[i] = 21;
		}
		else if (buffin[i].IData > -th1 && buffin[i].IData <= 0 && buffin[i].QData > th2)
		{
			buffout[i] = 17;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData <= -th2)
		{
			buffout[i] = 9;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData > -th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 13;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 15;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 7;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData > th1 && buffin[i].QData <= th2)
		{
			buffout[i] = 6;
		}
		else if (buffin[i].IData > 0 && buffin[i].IData <= th1 && buffin[i].QData > th2)
		{
			buffout[i] = 3;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData <= -th2)
		{
			buffout[i] = 8;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData > -th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 12;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 14;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 5;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData > th1 && buffin[i].QData <= th2)
		{
			buffout[i] = 4;
		}
		else if (buffin[i].IData > th1 && buffin[i].IData <= th2 && buffin[i].QData > th2)
		{
			buffout[i] = 2;
		}
		else if (buffin[i].IData > th2 && buffin[i].QData <= -th1)
		{
			buffout[i] = 10;
		}
		else if (buffin[i].IData > th2 && buffin[i].QData > -th1 && buffin[i].QData <= 0)
		{
			buffout[i] = 11;
		}
		else if (buffin[i].IData > th2 && buffin[i].QData > 0 && buffin[i].QData <= th1)
		{
			buffout[i] = 1;
		}
		else if (buffin[i].IData > th2 && buffin[i].QData > th1)
		{
			buffout[i] = 0;
		}
	}

}

void Demo_32QAM::diff_code_32qam(int* y, int* c_symbol, int len)
{
	int c = 0, i = 0, j = 0;
	int n = 2;

	int** a = new int* [len];
	for (i = 0; i < len; ++i)
	{
		a[i] = new int[5];
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
			a[i][0] = 0; a[i][1] = 0; a[i][2] = 0; a[i][3] = 0; a[i][4] = 0;
		}
		else if (j == 1)
		{
			a[i][1] = 0; a[i][2] = 0; a[i][3] = 0; a[i][4] = 0;
		}
		else if (j == 2)
		{
			a[i][2] = 0; a[i][3] = 0; a[i][4] = 0;
		}
		else if (j == 3)
		{
			a[i][3] = 0; a[i][4] = 0;
		}
		else if (j == 4)
		{
			a[i][4] = 0;
		}
	}

	int* c_symbol_2 = new int[len];
	memset(c_symbol_2, 0, sizeof(int) * len);

	for (i = 0; i < len; i++)
	{
		if (a[i][4] == 0 && a[i][3] == 0)
		{
			c_symbol_2[i] = 0;
		}
		else if (a[i][4] == 1 && a[i][3] == 0)
		{
			c_symbol_2[i] = 1;
		}
		else if (a[i][4] == 1 && a[i][3] == 1)
		{
			c_symbol_2[i] = 2;
		}
		else if (a[i][4] == 0 && a[i][3] == 1)
		{
			c_symbol_2[i] = 3;
		}
	}

	int** d_bit = new int* [len];
	for (i = 0; i < len; ++i)
	{
		d_bit[i] = new int[5];
	}

	for (i = 0; i < len; ++i)
	{

		d_bit[i][4] = a[i][0];
		d_bit[i][3] = a[i][1];
		d_bit[i][2] = a[i][2];
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
		d_bit_1[i] = d_bit[i][0] * 10000 + d_bit[i][1] * 1000 + d_bit[i][2] * 100 + d_bit[i][3] * 10 + d_bit[i][4] * 1;
	}

	int remainder = 0;
	for (j = 0; j < len; j++)
	{
		i = 0;
		c_symbol[j] = 0;
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

void Demo_32QAM::InitBlockCFilter()
{
	if (m_Rb <= 2e6)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 33;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		//float tempCoef[] =
		//{
		//	-0.0424636059371430	,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
		//	0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430
		//};//32QAM

		float tempCoef[] =
		{
			-0.00249658235788025,7.65898743126971e-18,0.00326476154492033,3.71494110574812e-18,-0.00445194756125500,
			1.29973681111812e-17,0.00643059092181277,-4.38550344176276e-18,-0.0101052143057058,5.68491186895173e-18,
			0.0181893857502704,-8.12130266993104e-18,-0.0424419000839643,1.46183448058759e-17,0.212209500419822,
			0.500006855655654,0.636628501259465,0.500006855655654,0.212209500419822,1.46183448058759e-17,
			-0.0424419000839643,-8.12130266993104e-18,0.0181893857502704,5.68491186895173e-18,-0.0101052143057058,
			-4.38550344176276e-18,0.00643059092181277,1.29973681111812e-17,-0.00445194756125500,3.71494110574812e-18,
			0.00326476154492033,7.65898743126971e-18,-0.00249658235788025
		};//32QAM


		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 33;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.00249658235788025,7.65898743126971e-18,0.00326476154492033,3.71494110574812e-18,-0.00445194756125500,
			1.29973681111812e-17,0.00643059092181277,-4.38550344176276e-18,-0.0101052143057058,5.68491186895173e-18,
			0.0181893857502704,-8.12130266993104e-18,-0.0424419000839643,1.46183448058759e-17,0.212209500419822,
			0.500006855655654,0.636628501259465,0.500006855655654,0.212209500419822,1.46183448058759e-17,
			-0.0424419000839643,-8.12130266993104e-18,0.0181893857502704,5.68491186895173e-18,-0.0101052143057058,
			-4.38550344176276e-18,0.00643059092181277,1.29973681111812e-17,-0.00445194756125500,3.71494110574812e-18,
			0.00326476154492033,7.65898743126971e-18,-0.00249658235788025
		};//32QAM


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

void Demo_32QAM::InitialDemodulation(int Rb) {

	m_Rb = Rb;

	InitBlockCFilter();

	QAM32Init();

	InitCostasPLL();

	InitSymbolSync();
}

void Demo_32QAM::QAM32Init()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_QAM.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;												// 降采样时的时钟控制  
	fDownConversionPhase = 0;										// 下变频中的相位值
	fPLLNCO = 0;													// 锁相环中的本地NCO
	fPLLPastFreqPart = 0;											// 锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;

	// add on 20240621
	flag = 0;
	slice = 10;
	QAM32_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[QAM32_Samplesize];

	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;								//符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;								//符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;							//符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;								//符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;									//符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;								//符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;								//符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;							//符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;							//符号同步time_error缓存值
}

void Demo_32QAM::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 2.333098397143958e-04;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 2.722496259248146e-08;

	De_QAM.m_sAlgDemInit.nSampPerSymb = 4;							//每个码元中的采样点数
	m_sAlgDemInit.nSampPerSymb = 4;									//每个码元中的采样点数

	De_QAM.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_32QAM::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;											//环路阻尼系数
	int Ko = 1;														//压控振荡器增益
	int Kd = 1;														//鉴相器增益
	float fBLcoef = 0.0006;											//32QAM
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_QAM.CostasPll = CostasPll;
}

void Demo_32QAM::AGCQAM32(Complex* BuffIn, float target_power, int len)
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


float Demo_32QAM::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}