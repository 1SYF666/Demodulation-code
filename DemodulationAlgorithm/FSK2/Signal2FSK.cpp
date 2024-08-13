#include"Signal2FSK.h"

Demo_2FSK::Demo_2FSK(int index) :m_index(index)
{

}

Demo_2FSK::~Demo_2FSK()
{

}

bool Demo_2FSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;

	InitialDemodulation(demoParameter);
	return 1;
}

void Demo_2FSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	memcpy(Databuff0 + flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;

	Demodulation(Databuff0, baseTime, demodulationResult);


}

void Demo_2FSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out)
{


	//FILE*fp;
	//char file_path1[500];
	/*sprintf(file_path1, "D:/data/input_I.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FSK_Samplesize; i++) {
		fprintf(fp, "%f\n", data_input_slice[i].IData);
	}
	fclose(fp);

	FILE*fp2;
	char file_path2[500];
	sprintf(file_path2, "D:/data/input_Q.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < FSK_Samplesize; i++) {
		fprintf(fp2, "%f\n", data_input_slice[i].QData);
	}
	fclose(fp2);*/


	//Complex* DataBuff = new Complex[FSK_Samplesize];
	//memset(DataBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	DataBuff[i] = data_input_slice[i];
	//}

	//float* Signal_sort = new float[FSK_Samplesize];
	//memset(Signal_sort, 0x00, sizeof(float) * FSK_Samplesize);

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	Signal_sort[i] = fabs(DataBuff[i].IData);
	//}

	//De_FSK.Quick_Sort(Signal_sort, 0, FSK_Samplesize - 1);

	//float mean_power = mean_function(Signal_sort, FSK_Samplesize - 70, FSK_Samplesize - 40);

	AGCFSK2(data_input_slice, fAGCPastVc, FSK_Samplesize); // add on 20240611 

	//Complex* DataBuff_normalize = new Complex[FSK_Samplesize];
	//memset(DataBuff_normalize, 0x00, sizeof(Complex) * FSK_Samplesize);

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	DataBuff_normalize[i].IData = DataBuff[i].IData / mean_power;
	//	DataBuff_normalize[i].QData = DataBuff[i].QData / mean_power;
	//}


	// De_FSK.AGC(DataBuff_normalize, fAGCPastVc);

	/*sprintf(file_path1, "D:/data/DataBuff_normalizeI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FSK_Samplesize; i++) {
		fprintf(fp, "%f\n", DataBuff_normalize[i].IData);
	}
	fclose(fp);


	sprintf(file_path2, "D:/data/DataBuff_normalizeQ.txt");
	fp2 = fopen(file_path2, "at");
	for (int i = 0; i < FSK_Samplesize; i++) {
		fprintf(fp2, "%f\n", DataBuff_normalize[i].QData);
	}
	fclose(fp2);*/

	if (0)
	{
		Complex* DataPLLBuff = new Complex[FSK_Samplesize];
		memset(DataPLLBuff, 0x00, sizeof(Complex) * FSK_Samplesize);

		float* fPLLPastFreqPart = new float[FSK_Samplesize];
		memset(fPLLPastFreqPart, 0x00, sizeof(float) * FSK_Samplesize);
		float* fPLLPastFreqPart_abs = new float[FSK_Samplesize];
		memset(fPLLPastFreqPart_abs, 0x00, sizeof(float) * FSK_Samplesize);

		De_FSK.CarrierSync_FSK(data_input_slice, DataPLLBuff, FSK_Samplesize, 1, &fPLLNCO, fPLLPastFreqPart, &fPLLPastFreqPartyemp);

		/*sprintf(file_path1, "D:/data/fPLLPastFreqPart.txt");
		fp = fopen(file_path1, "at");
		for (int i = 0; i < FSK_Samplesize; i++) {
			fprintf(fp, "%f\n", fPLLPastFreqPart[i]);
		}
		fclose(fp);*/


		float mean_fPLL = mean_function(fPLLPastFreqPart, 0, FSK_Samplesize);

		for (int i = 0; i < FSK_Samplesize; i++)
		{
			fPLLPastFreqPart[i] = fPLLPastFreqPart[i] - mean_fPLL;
		}

		/*sprintf(file_path1, "D:/data/fPLLPastFreqPartnormalization.txt");
		fp = fopen(file_path1, "at");
		for (int i = 0; i < FSK_Samplesize; i++)
		{
			fprintf(fp, "%f\n", fPLLPastFreqPart[i]);
		}
		fclose(fp);*/
	}

	// modify on 20240518
	Complex* DataPLLBuff = new Complex[FSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	Complex* DataFLLBuff = new Complex[FSK_Samplesize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * FSK_Samplesize);

	float* DataFLLBuff_normalize1 = new float[FSK_Samplesize];
	memset(DataFLLBuff_normalize1, 0x00, sizeof(float) * FSK_Samplesize);
	float* DataFLLBuff_normalize2 = new float[FSK_Samplesize];
	memset(DataFLLBuff_normalize2, 0x00, sizeof(float) * FSK_Samplesize);

	De_FSK.CarrierSync(data_input_slice, DataPLLBuff, DataFLLBuff, FSK_Samplesize, 9, &fPLLNCO, &fPLLPastFreqPart1, &FLLBuffer);

	for (int i = 0; i < FSK_Samplesize; i++)
	{
		DataFLLBuff_normalize1[i] = fabs(DataFLLBuff[i].IData);
		DataFLLBuff_normalize2[i] = DataFLLBuff[i].IData;
	}
	float max = max_function(DataFLLBuff_normalize1, 0, FSK_Samplesize);
	float mean_FLL = mean_function(DataFLLBuff_normalize2, 0, FSK_Samplesize);

	for (int i = 0; i < FSK_Samplesize; i++)
	{
		DataFLLBuff[i].IData = (DataFLLBuff[i].IData - mean_FLL) / max;
		DataFLLBuff[i].QData = (DataFLLBuff[i].QData - mean_FLL) / max;
	}

	/*sprintf(file_path1, "D:/data/DataFLLBuffI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", DataFLLBuff[i].IData);
	}
	fclose(fp);*/



	Complex* DataBuffin_Filter = new Complex[FSK_Samplesize];
	memset(DataBuffin_Filter, 0x00, sizeof(Complex) * FSK_Samplesize);
	Complex* DataBuffout_Filter = new Complex[FSK_Samplesize];
	memset(DataBuffout_Filter, 0x00, sizeof(Complex) * FSK_Samplesize);


	for (int i = 0; i < FSK_Samplesize; i++)
	{
		DataBuffin_Filter[i] = DataFLLBuff[i];
	}

	De_FSK.BlockFilter(DataBuffin_Filter, FSK_Samplesize, DataBuffout_Filter, cFilterBuff);

	/*sprintf(file_path1, "D:/data/DataBuffout_Filter.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < FSK_Samplesize; i++)
	{
		fprintf(fp, "%f\n", DataBuffout_Filter[i].IData);
	}
	fclose(fp);*/

	// 将采样率fs变换成4Rb,已进入符号同步
	//int nSymbolSyncSize = floor(m_SamleSize * ((4 * m_Rb) / m_fs));
	//int *temp_down = new int[nSymbolSyncSize];
	//memset(temp_down, 0x00, sizeof(int) * nSymbolSyncSize);
	//Complex* DataSymbolSyncBuff = new Complex[nSymbolSyncSize];
	//memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * nSymbolSyncSize);

	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	
	//	if (i == 0)
	//	{
	//		temp_down[i] = 0;
	//		DataSymbolSyncBuff[i].IData = DataBuffout_Filter[i].IData;
	//		DataSymbolSyncBuff[i].QData = DataBuffout_Filter[i].IData;
	//	}
	//	else
	//	{
	//		temp_down[i] = round(i * m_fs / (4 * m_Rb));
	//		DataSymbolSyncBuff[i].IData = DataBuffout_Filter[temp_down[i]].IData;
	//		DataSymbolSyncBuff[i].QData = DataBuffout_Filter[temp_down[i]].IData;
	//	}
	//}
	//Complex* DataSymbolSyncBuffout = new Complex[nSymbolSyncSize];
	//memset(DataSymbolSyncBuffout, 0x00, sizeof(Complex) * nSymbolSyncSize);
	//int nSymbolSync_count = 0;
	//De_FSK.SymbolSync(DataSymbolSyncBuff, DataSymbolSyncBuffout, nSymbolSyncSize, &nSymbolSync_count, &SySyncBuffer);

	// add on 20240511
	Complex* DataSymbolSyncBuff = new Complex[FSK_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	double step = ((double)(4 * m_Rb)) / ((double)m_fs);
	int nSymbolSyncSize = 0;
	Downfs(DataBuffout_Filter, DataSymbolSyncBuff, step, FSK_Samplesize, &nSymbolSyncSize);

	/*sprintf(file_path1, "D:/data/DataSymbolSyncBuff.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	}
	fclose(fp);*/

	Complex* DataSymbolSyncBuffout = new Complex[nSymbolSyncSize];
	memset(DataSymbolSyncBuffout, 0x00, sizeof(Complex) * nSymbolSyncSize);
	
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].IData;
	//}

	int nSymbolSync_count = 0;
	De_FSK.SymbolSync(DataSymbolSyncBuff, DataSymbolSyncBuffout, nSymbolSyncSize, &nSymbolSync_count, &SySyncBuffer);

	/*sprintf(file_path1, "D:/data/DataSymbolSyncBuffout.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < nSymbolSync_count; i++)
	{
		fprintf(fp, "%f\n", DataSymbolSyncBuffout[i].IData);
	}
	fclose(fp);*/



	// 画星座图

	//调用星座图绘制接口
	int star_len = nSymbolSync_count;

	Complex* DataPLLBuff_star = new Complex[star_len];
	memset(DataPLLBuff_star, 0x00, sizeof(Complex) * star_len);

	for (int i = 0; i < star_len; i++)
	{
		DataPLLBuff_star[i].IData = DataSymbolSyncBuffout[i].IData;
		DataPLLBuff_star[i].QData = DataSymbolSyncBuffout[i].IData;

	}
	signal_16APSK_Out->burstData->softDistinguishData = new Complex[star_len];
	memcpy(signal_16APSK_Out->burstData->softDistinguishData, DataPLLBuff_star, (sizeof(Complex))*star_len);
	signal_16APSK_Out->burstData->softDistinguishDataLen = star_len;

	// 判决
	char* c_symbol = new char[nSymbolSync_count];
	Judgment_2fsk(DataSymbolSyncBuffout, nSymbolSync_count, c_symbol);

	//sprintf(file_path1, "D:/data/FSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSync_count; i++)
	//{
	//	fprintf(fp, "%d\n", c_symbol[i]);
	//}
	//fclose(fp);

	signal_16APSK_Out->burstData->nDemodulationByteI = new char[nSymbolSync_count];
	memcpy(signal_16APSK_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSync_count);



	//DELETE_ARR(DataBuff);
	//DELETE_ARR(Signal_sort);
	//DELETE_ARR(DataBuff_normalize);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataPLLBuff_star);
	DELETE_ARR(DataFLLBuff_normalize1);
	DELETE_ARR(DataFLLBuff_normalize2);
	DELETE_ARR(DataBuffin_Filter);
	DELETE_ARR(DataBuffout_Filter);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(DataSymbolSyncBuffout);
	DELETE_ARR(c_symbol);
    DELETE_ARR(DataFLLBuff);
}

float Demo_2FSK::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

float Demo_2FSK::max_function(float* data, int start, int end)
{
	int i;
	float max1 = data[start];

	for (i = start + 1; i < end; i++)
	{
		if (max1 < data[i])
		{
			max1 = data[i];
		}
	}
	return max1;
}

int Demo_2FSK::max_index_function(float* data, int start, int end)
{
	int i;
	int index = start;
	float max1 = data[start];
	for (i = start + 1; i < end; i++)
	{
		if (max1 < data[i])
		{
			max1 = data[i];
			index = i;
		}
	}
	return index;
}

void Demo_2FSK::Judgment_2fsk(Complex buffin[], int bufflen, char buffout[])
{


	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData >0)
		{
			buffout[i] = 0;
		}
		else
		{
			buffout[i] = 1;
		}
	}

}

void Demo_2FSK::InitBlockCFilter()
{
	if (m_Rb <= 4e5)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 16 + 1;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0101060032508123, 5.68535570747550e-18, 0.0181908058514621, -8.12193672496501e-18 ,-0.0424452136534116,
			1.46194861049370e-17, 0.212226068267058, 0.500045892726026,	0.636678204801174, 0.500045892726026, 0.212226068267058,
			1.46194861049370e-17, -0.0424452136534116, -8.12193672496501e-18, 0.0181908058514621, 5.68535570747550e-18 ,-0.0101060032508123
		};//4FSK


		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{
		CFilter.nFilterTaps = 20 + 1;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0268661128107165, -0.0272597932734082, -0.0134774114144576, 0.0182063277002322	,0.0684930940231012
			,0.134330564053582,	0.208991748429463,	0.283025639703609	,0.345920226304412,	0.388127532797573
			,0.402991692160747,	0.388127532797573,0.345920226304412,0.283025639703609,0.208991748429463	,0.134330564053582
			,0.0684930940231012	,0.0182063277002322 ,-0.0134774114144576, -0.0272597932734082 ,-0.0268661128107165
		};//4FSK

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}

	}

	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);

	CFilter.bCoefEvenSym = false;

	De_FSK.CFilter = CFilter;
}

void Demo_2FSK::InitialDemodulation(const DemodulationInitParamater& signalinit)
{

	m_Rb = signalinit.rb;
	m_fs = signalinit.fs;
	memcpy(m_fc, signalinit.fc_demo, sizeof(float) * 8);

	InitBlockCFilter();

	FSK2Init();

	//InitCostasPLL();

	InitSymbolSync();

	//add on 20240518
	InitFLL();
}

void Demo_2FSK::FSK2Init()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_FSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;												//降采样时的时钟控制  
	fDownConversionPhase = 0;										//下变频中的相位值
	nPLLBuffSize = 0;
	
	//add on 20240514
	flag = 0;
	slice = 3;
	FSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[FSK_Samplesize];

	// add on 20240514
	//FLL初始化
	FLLBuffer.CFLLbuff = new Complex[4 / 2];
	memset(FLLBuffer.CFLLbuff, 0x00, sizeof(Complex) * (4 / 2));	               //nSampPerSymb/2长度  //FLL中间缓存区
	FLLBuffer.fFLLIaccum = 0;                                   //FLL中I_accum的缓存值
	FLLBuffer.fFLLQaccum = 0;                                   //FLL中Q_accum的缓存值
	FLLBuffer.fFLLfreqoutdlf = 0;                               //FLL中freqoutdlf的上一个缓存值
	FLLBuffer.fFLLPastfreqoutdlf = 0;                           //FLL中freqoutdlf的前两个缓存值
	FLLBuffer.fFLLfreqerror = 0;                                //FLL中freqerror的上一个缓存值
	FLLBuffer.fFLLNCOPhase = 0;                                 //FLL中NCOPhase的上一个缓存值

	fPLLNCO = 0;                                                //锁相环中的本地NCO
	fPLLPastFreqPart1 = 0;                                       //锁相环中的频率跟踪曲线


	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;								//符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.8;								//符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.8;							//符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;								//符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;									//符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;								//符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;								//符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;							//符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;							//符号同步time_error缓存值
}

void Demo_2FSK::InitSymbolSync()
{

	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.099989931306170/100;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 0.005000503333313/10000;	// 2FSK 202400518

	De_FSK.m_sAlgDemInit.nSampPerSymb = 4;							// 每个码元中的采样点数
	m_sAlgDemInit.nSampPerSymb = 4;									// 每个码元中的采样点数

	De_FSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_2FSK::InitCostasPLL()
{
	float BL = 0;
	float Wn = 0;
	float T_nco = 0;
	float sigma = 0.707;											//环路阻尼系数
	int Ko = 1;														//压控振荡器增益
	int Kd = 1;														//鉴相器增益
	float fBLcoef = 0.28;											//4FSK 
	int K = Ko * Kd;

	if ((m_fc[0] > 0 && m_fc[1] < 0 )||( m_fc[0] < 0 && m_fc[1]>0))
	{

		BL = fBLcoef * (m_Rb + fabs(m_fc[1]) + fabs(m_fc[0]));
	}
	else
	{
		BL = fBLcoef * (m_Rb + fabs(m_fc[1] - m_fc[0]));
	}

	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));

	T_nco = 1 / m_fs;

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_FSK.CostasPll = CostasPll;

	// modify on 2024/05/14
	fPLLPastFreqPartyemp = 0;										//锁相环中的频率跟踪中间变量
	fPLLNCO = 0;													//锁相环中的本地NCO											
}

// add on 2024/05/11

// 输入数据input
// 输出数据output
// 采样率变换倍数step = fs_after / fs_before
// 输入信号长度lengthin
// 输出信号长度lengthout
void Demo_2FSK::Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout)
{
	int i = 0;
	int j = 0;
	int k = 0;
	static double step_temp = 0;
	while(1)
	{
		while (1)
		{
			step_temp += step;

			if (step_temp >= 1)
			{
				step_temp -= 1;
				i += j;
				j = 1;
				break;
			}
			else
			{
				j++;
				if ((i + j) > lengthin)
				{
					break;
				}				
			}
		}

		if ((i + j) > lengthin)
		{
			break;
		}
		else
		{
			//output[k] = input[k];
			output[k].IData = input[i].IData;
			output[k].QData = input[i].QData;
			k++;
		}

	} 

	*lengthout = k;
}

//add on 2024/05/18
void Demo_2FSK::InitFLL()
{
	float fSymbolTime = 1.0 / (float)m_Rb;
	//float fDownSampleFs = m_Rb * (float)4;
	//Fll.fFLLDownSampleTs = 1.0 / ((float)(m_Rb * 4));
	float fDownSampleFs = m_fs;								//modify on 20240518
	Fll.fFLLDownSampleTs = 1.0 / float(fDownSampleFs);
	float fCarrierSyncFactor = 0.000002;
	float b_fll = fCarrierSyncFactor * m_Rb;
	float wn = 1.89 * b_fll;
	float gain_nco = 2 * Pi * fDownSampleFs;
	int gain_disctim = 1;
	int K = gain_nco * gain_disctim;
	Fll.fFLLPLLLoopFilterCoef1 = (sqrt(2) * wn * fSymbolTime + (wn * wn) * (fSymbolTime * fSymbolTime)) / K;
	Fll.fFLLPLLLoopFilterCoef2 = sqrt(2) * wn * fSymbolTime / K;
	De_FSK.Fll = Fll;
}


void Demo_2FSK::AGCFSK2(Complex* BuffIn, float target_power, int len)
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

float Demo_2FSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}