#include "Signal64QAM.h"

Demo_64QAM::Demo_64QAM(int index) :m_index(index)
{

}

Demo_64QAM::~Demo_64QAM() {

}

bool Demo_64QAM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter.rb);
	return 1;
}

void Demo_64QAM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_64QAM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_64QAM_Out)
{
	// add on 20240625
	if (nSilceCount==1)
	{
		// 第二次进入，修改PLL参数
		float fBLcoef = 0.0001;
		InitCostasPLL(&fBLcoef);
	}
	nSilceCount++;

	//FILE* fp2;
	//char file_path2[500];


	Complex* DataBuff = new Complex[QAM64_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * QAM64_Samplesize);

	De_QAM.BlockFilter(data_input_slice, QAM64_Samplesize, DataBuff, cFilterBuff);

	//De_QAM.AGC(DataBuff, fAGCPastVc);

	AGCQAM64(DataBuff, fAGCPastVc, QAM64_Samplesize);


	int nSymbolSyncSize = 0;//突发在当前片中的长度
	Complex* DataSymbolSyncBuff = new Complex[QAM64_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * QAM64_Samplesize);

	De_QAM.SymbolSync(DataBuff, DataSymbolSyncBuff, QAM64_Samplesize, &nSymbolSyncSize, &SySyncBuffer);

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataSymbolSyncBuff[i].IData = DataSymbolSyncBuff[i].IData * 7;
		DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].QData * 7;
	}

	Complex* DataPLLBuff = new Complex[nSymbolSyncSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);

	De_QAM.CarrierSync_32_or_64qam(DataSymbolSyncBuff, DataPLLBuff, nSymbolSyncSize, 2, &fPLLNCO, &fPLLPastFreqPart);


	float* II = new float[nSymbolSyncSize];
	memset(II, 0x00, sizeof(float) * nSymbolSyncSize);
	float* QQ = new float[nSymbolSyncSize];
	memset(QQ, 0x00, sizeof(float) * nSymbolSyncSize);
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		II[i] = DataPLLBuff[i].IData;
		QQ[i] = DataPLLBuff[i].QData;
	}
	float mean_I = mean_function(II, 0, nSymbolSyncSize);
	float mean_Q = mean_function(QQ, 0, nSymbolSyncSize);
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataPLLBuff[i].IData = DataPLLBuff[i].IData - mean_I;
		DataPLLBuff[i].QData = DataPLLBuff[i].QData - mean_Q;
	}

	float* Signal_sort = new float[nSymbolSyncSize];
	memset(Signal_sort, 0x00, sizeof(float) * nSymbolSyncSize);
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		Signal_sort[i] = fabs(DataPLLBuff[i].IData);
	}
	De_QAM.Quick_Sort(Signal_sort, 0, nSymbolSyncSize - 1);
	float mean_power = mean_function(Signal_sort, nSymbolSyncSize * 10.5 / 12, nSymbolSyncSize * 11.5 / 12);
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataPLLBuff[i].IData = DataPLLBuff[i].IData / mean_power;
		DataPLLBuff[i].QData = DataPLLBuff[i].QData / mean_power;
	}

	Complex* Data_phase = new Complex[nSymbolSyncSize];
	Phase_correction(DataPLLBuff, Data_phase, nSymbolSyncSize);


	int starmapcount = nSymbolSyncSize - 500;
	Complex* starData = new Complex[starmapcount];
	memset(starData, 0x00, sizeof(Complex) * starmapcount);
	for (int i = 0; i < starmapcount; i++)
	{
		starData[i] = Data_phase[i + 500];
	}

	signal_64QAM_Out->burstData->softDistinguishData = new Complex[starmapcount];
	memcpy(signal_64QAM_Out->burstData->softDistinguishData, DataPLLBuff, (sizeof(Complex)) * starmapcount);
	signal_64QAM_Out->burstData->softDistinguishDataLen = starmapcount;

	int* y1 = new int[nSymbolSyncSize];
	Judgment_64qam(Data_phase, nSymbolSyncSize, y1);

	char* c_symbol = new char[nSymbolSyncSize];
	memset(c_symbol, 0x00, sizeof(char) * nSymbolSyncSize);
	diff_code_64qam(y1, c_symbol, nSymbolSyncSize);

	signal_64QAM_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_64QAM_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);



	//sprintf(file_path2, "D:/data/QAM64_output_symbol.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp2, "%d\n", c_symbol[i]);
	//}
	//fclose(fp2);



	DELETE_ARR(DataBuff);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(II);
	DELETE_ARR(QQ);
	DELETE_ARR(Signal_sort);
	DELETE_ARR(Data_phase);
	DELETE_ARR(y1);
	DELETE_ARR(c_symbol);
	DELETE_ARR(starData);

#pragma endregion

}

void Demo_64QAM::InitialDemodulation(int Rb)
{
	float fBLcoef = 0.001;
	m_Rb = Rb;
	QAM64_Init();
	InitBlockCFilter();
	InitCostasPLL(&fBLcoef);
	InitSymbolSync();
}

void Demo_64QAM::InitBlockCFilter()
{
	if (m_Rb <= 2e6)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 17;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			//-0.0424636059371430	,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
			//0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430

			 /*- 0.0424636059371430,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146,
			 0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430 */

			//-0.00249658235788025,7.65898743126971e-18,0.00326476154492033,3.71494110574812e-18,-0.00445194756125500,
			//1.29973681111812e-17,0.00643059092181277,-4.38550344176276e-18,-0.0101052143057058,5.68491186895173e-18,
			//0.0181893857502704,-8.12130266993104e-18,-0.0424419000839643,1.46183448058759e-17,0.212209500419822,
			//0.500006855655654,0.636628501259465,0.500006855655654,0.212209500419822,1.46183448058759e-17,
			//-0.0424419000839643,-8.12130266993104e-18,0.0181893857502704,5.68491186895173e-18,-0.0101052143057058,
			//-4.38550344176276e-18,0.00643059092181277,1.29973681111812e-17,-0.00445194756125500,3.71494110574812e-18,
			//0.00326476154492033,7.65898743126971e-18,-0.00249658235788025

			- 0.0101060032508123,5.68535570747550e-18,0.0181908058514621,-8.12193672496501e-18,-0.0424452136534116,
			1.46194861049370e-17,0.212226068267058,0.500045892726026,0.636678204801174,0.500045892726026,0.212226068267058,
			1.46194861049370e-17,-0.0424452136534116,-8.12193672496501e-18,0.0181908058514621,5.68535570747550e-18,-0.0101060032508123

		};

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 17;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			//-0.0424636059371430	,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146	,
			//0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430

			 /*- 0.0424636059371430,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146,
			 0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430 */

			 //-0.00249658235788025,7.65898743126971e-18,0.00326476154492033,3.71494110574812e-18,-0.00445194756125500,
			 //1.29973681111812e-17,0.00643059092181277,-4.38550344176276e-18,-0.0101052143057058,5.68491186895173e-18,
			 //0.0181893857502704,-8.12130266993104e-18,-0.0424419000839643,1.46183448058759e-17,0.212209500419822,
			 //0.500006855655654,0.636628501259465,0.500006855655654,0.212209500419822,1.46183448058759e-17,
			 //-0.0424419000839643,-8.12130266993104e-18,0.0181893857502704,5.68491186895173e-18,-0.0101052143057058,
			 //-4.38550344176276e-18,0.00643059092181277,1.29973681111812e-17,-0.00445194756125500,3.71494110574812e-18,
			 //0.00326476154492033,7.65898743126971e-18,-0.00249658235788025

			 -0.0101060032508123,5.68535570747550e-18,0.0181908058514621,-8.12193672496501e-18,-0.0424452136534116,
			 1.46194861049370e-17,0.212226068267058,0.500045892726026,0.636678204801174,0.500045892726026,0.212226068267058,
			 1.46194861049370e-17,-0.0424452136534116,-8.12193672496501e-18,0.0181908058514621,5.68535570747550e-18,-0.0101060032508123

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

void Demo_64QAM::QAM64_Init() {
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
	slice = 10;
	QAM64_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[QAM64_Samplesize];

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

void Demo_64QAM::InitSymbolSync()
{
	//symbolsyncsactorInit.fSymbolSyncFactor1 = 0.001333199084082;
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667416e-07;
	symbolsyncsactorInit.fSymbolSyncFactor1 = 3.332997710205654e-04;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 5.556114814792135e-08;
	De_QAM.m_sAlgDemInit.nSampPerSymb = 4;					//每个码元中的采样点数
	De_QAM.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_64QAM::InitCostasPLL(float* fBLcoef) {
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;									//环路阻尼系数
	int Ko = 1;												//压控振荡器增益
	int Kd = 1;												//鉴相器增益
	//float fBLcoef = 0.0001;
	int K = Ko * Kd;
	BL = *fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);
	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_QAM.CostasPll = CostasPll;
}




void Demo_64QAM::Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen)
{
	float* theta_0 = new float[bufflen];
	memset(theta_0, 0x00, sizeof(float) * bufflen);

	int n1 = 0;

	for (int i = 0; i < bufflen; i++)
	{
		float power_s = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
		// (sqrt(1.0 / 3 * 1.0 / 3 * 2) + sqrt(1.0 / 3 * 1.0 / 3 + 1 * 1)) / 2.0
		// =  0.7627
		if (power_s < (sqrt(1.0 / 7 * 1.0 / 7 * 2) + sqrt(1.0 / 7 * 1.0 / 7 + 3.0 / 7 * 3.0 / 7)) / 2.0)
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
		else if (theta_0[i] > Pi / 2 && theta_0[i] <= Pi / 2 * 2)
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
	//float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 3 / 8, (n1 - 1) * 5 / 8); //32QAM
	//float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 3 / 8, (n1 - 1) * 5 / 8); //64QAM
	float mean_sort_theta_1 = mean_function(sort_theta_1, (n1 - 1) * 7 / 16, (n1 - 1) * 9 / 16); //64QAM
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

void Demo_64QAM::Judgment_64qam(Complex buffin[], int bufflen, int buffout[]) {
	float th1 = 6.0 / 7;
	float th2 = 4.0 / 7;
	float th3 = 2.0 / 7;
	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData >= th1 && buffin[i].QData >= th1)
		{
			buffout[i] = 0;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 1;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 2;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 3;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= th1)
		{
			buffout[i] = 4;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= th1)
		{
			buffout[i] = 5;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= th1)
		{
			buffout[i] = 6;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 7;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 8;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 9;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 10;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 11;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 12;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 13;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 14;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 15;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData < -th1)
		{
			buffout[i] = 16;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 20;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 21;
		}
		else if (buffin[i].IData >= th1 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 22;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData < -th1)
		{
			buffout[i] = 17;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData < -th1)
		{
			buffout[i] = 18;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData < -th1)
		{
			buffout[i] = 19;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 23;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 26;
		}
		else if (buffin[i].IData >= th2 && buffin[i].IData < th1 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 29;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 24;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 27;
		}
		else if (buffin[i].IData >= th3 && buffin[i].IData < th2 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 30;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 25;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 28;
		}
		else if (buffin[i].IData >= 0 && buffin[i].IData < th3 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 31;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= th1)
		{
			buffout[i] = 32;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= th1)
		{
			buffout[i] = 33;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= th1)
		{
			buffout[i] = 34;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= th1)
		{
			buffout[i] = 35;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 36;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 37;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 38;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 39;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 40;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= th2 && buffin[i].QData < th1)
		{
			buffout[i] = 41;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 42;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 43;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= th3 && buffin[i].QData < th2)
		{
			buffout[i] = 44;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 45;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 46;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= 0 && buffin[i].QData < th3)
		{
			buffout[i] = 47;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData < -th1)
		{
			buffout[i] = 48;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData < -th1)
		{
			buffout[i] = 52;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData < -th1)
		{
			buffout[i] = 53;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData < -th1)
		{
			buffout[i] = 54;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 49;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 50;
		}
		else if (buffin[i].IData < -th1 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 51;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 55;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 58;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= -th1 && buffin[i].QData < -th2)
		{
			buffout[i] = 61;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 56;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 59;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= -th2 && buffin[i].QData < -th3)
		{
			buffout[i] = 62;
		}
		else if (buffin[i].IData >= -th1 && buffin[i].IData < -th2 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 57;
		}
		else if (buffin[i].IData >= -th2 && buffin[i].IData < -th3 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 60;
		}
		else if (buffin[i].IData >= -th3 && buffin[i].IData < 0 && buffin[i].QData >= -th3 && buffin[i].QData < 0)
		{
			buffout[i] = 63;
		}
		else
		{
			buffout[i] = -1;
		}

	}
}

float Demo_64QAM::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_64QAM::diff_code_64qam(int* y, char* c_symbol, int len)
{
	int c = 0, i = 0, j = 0;
	int n = 2;

	int** a = new int* [len];
	for (i = 0; i < len; ++i)
	{
		a[i] = new int[6];
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
			a[i][0] = 0; a[i][1] = 0; a[i][2] = 0; a[i][3] = 0; a[i][4] = 0; a[i][5] = 0;
		}
		else if (j == 1)
		{
			a[i][1] = 0; a[i][2] = 0; a[i][3] = 0; a[i][4] = 0; a[i][5] = 0;
		}
		else if (j == 2)
		{
			a[i][2] = 0; a[i][3] = 0; a[i][4] = 0; a[i][5] = 0;
		}
		else if (j == 3)
		{
			a[i][3] = 0; a[i][4] = 0; a[i][5] = 0;
		}
		else if (j == 4)
		{
			a[i][4] = 0; a[i][5] = 0;
		}
		else if (j == 5)
		{
			a[i][5] = 0;
		}
	}

	int* c_symbol_2 = new int[len];
	memset(c_symbol_2, 0, sizeof(int) * len);

	for (i = 0; i < len; i++)
	{
		if (a[i][5] == 0 && a[i][4] == 0)
		{
			c_symbol_2[i] = 0;
		}
		else if (a[i][5] == 1 && a[i][4] == 0)
		{
			c_symbol_2[i] = 1;
		}
		else if (a[i][5] == 1 && a[i][4] == 1)
		{
			c_symbol_2[i] = 2;
		}
		else if (a[i][5] == 0 && a[i][4] == 1)
		{
			c_symbol_2[i] = 3;
		}
	}

	int** d_bit = new int* [len];
	for (i = 0; i < len; ++i)
	{
		d_bit[i] = new int[6];
	}

	for (i = 0; i < len; ++i)
	{

		d_bit[i][5] = a[i][0];
		d_bit[i][4] = a[i][1];
		d_bit[i][3] = a[i][2];
		d_bit[i][2] = a[i][3];
	}

	int d_bit_2;

	for (i = 0; i < len; i++)
	{
		if (i == 0)
		{
			d_bit_2 = ((c_symbol_2[i] - difftemp) + 4) % 4;//mod要加阶数如，QP要+4
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
		d_bit_1[i] = d_bit[i][0] * 100000 + d_bit[i][1] * 10000 + d_bit[i][2] * 1000 + d_bit[i][3] * 100 + d_bit[i][4] * 10 + d_bit[i][5] * 1;
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


void Demo_64QAM::AGCQAM64(Complex* BuffIn, float target_power, int len)
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

float Demo_64QAM::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}