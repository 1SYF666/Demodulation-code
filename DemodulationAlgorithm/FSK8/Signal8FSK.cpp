#define _CRT_SECURE_NO_WARNINGS 1
#include"Signal8FSK.h"

Demo_8FSK::Demo_8FSK(int index) :m_index(index)
{

}

Demo_8FSK::~Demo_8FSK()
{
	DELETE_ARR(Databuff0);
}

bool Demo_8FSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter);
	return 1;
}

void Demo_8FSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_8FSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out)
{

	//FILE*fp;
	//char file_path1[500];
	//sprintf(file_path1, "D:/data/input_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FSK_Samplesize; i++) {
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//FILE*fp2;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/input_Q.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < FSK_Samplesize; i++) {
	//	fprintf(fp2, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp2);


	Complex* DataBuff = new Complex[FSK_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	DataBuff = data_input_slice;

	//float* Signal_sort = new float[FSK_Samplesize];
	//memset(Signal_sort, 0x00, sizeof(float) * FSK_Samplesize);

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	Signal_sort[i] = fabs(DataBuff[i].IData);
	//}

	//De_FSK.Quick_Sort(Signal_sort, 0, FSK_Samplesize - 1);

	//float mean_power = mean_function(Signal_sort, FSK_Samplesize - 100, FSK_Samplesize - 40);
	//Complex* DataBuff_normalize = new Complex[FSK_Samplesize];
	//memset(DataBuff_normalize, 0x00, sizeof(Complex) * FSK_Samplesize);

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	DataBuff_normalize[i].IData = DataBuff[i].IData / mean_power;
	//	DataBuff_normalize[i].QData = DataBuff[i].QData / mean_power;
	//}

	//De_FSK.AGC(DataBuff_normalize, fAGCPastVc);


	Complex* DataBuff_normalize = new Complex[FSK_Samplesize];
	memset(DataBuff_normalize, 0x00, sizeof(Complex) * FSK_Samplesize);
	for (int i = 0; i < FSK_Samplesize; i++)
	{
		DataBuff_normalize[i].IData = DataBuff[i].IData ;
		DataBuff_normalize[i].QData = DataBuff[i].QData ;
	}
	AGCFSK8(DataBuff_normalize, fAGCPastVc, FSK_Samplesize);  // add on 20240611


	//sprintf(file_path1, "D:/data/DataBuff_normalizeI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FSK_Samplesize; i++) {
	//	fprintf(fp, "%f\n", DataBuff_normalize[i].IData);
	//}
	//fclose(fp);


	//sprintf(file_path2, "D:/data/DataBuff_normalizeQ.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	fprintf(fp2, "%f\n", DataBuff_normalize[i].QData);
	//}
	//fclose(fp2);

	Complex* DataPLLBuff = new Complex[FSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	float* fPLLPastFreqPart = new float[FSK_Samplesize];
	memset(fPLLPastFreqPart, 0x00, sizeof(float) * FSK_Samplesize); //8FSK载波输出
	float* fPLLPastFreqPart_abs = new float[FSK_Samplesize];
	memset(fPLLPastFreqPart_abs, 0x00, sizeof(float) * FSK_Samplesize);

	De_FSK.CarrierSync_FSK(DataBuff_normalize, DataPLLBuff, FSK_Samplesize, 3, &fPLLNCO, fPLLPastFreqPart, &fPLLPastFreqPartemp); //8FSK

	//sprintf(file_path1, "D:/data/fPLLPastFreqPart.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FSK_Samplesize; i++) {
	//	fprintf(fp, "%f\n", fPLLPastFreqPart[i]);
	//}
	//fclose(fp);

	float mean_fPLL = mean_function(fPLLPastFreqPart, 0, FSK_Samplesize);

	for (int i = 0; i < FSK_Samplesize; i++)
	{
		fPLLPastFreqPart[i] = (fPLLPastFreqPart[i] - mean_fPLL);
	}

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	fPLLPastFreqPart_abs[i] = fabs(fPLLPastFreqPart[i]);
	//}

	//float max = max_function(fPLLPastFreqPart_abs, 0, FSK_Samplesize);

	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	fPLLPastFreqPart[i] = fPLLPastFreqPart[i] / max;
	//}

	//sprintf(file_path1, "D:/data/fPLLPastFreqPartnormalization.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", fPLLPastFreqPart[i]);
	//}
	//fclose(fp);

	Complex* DataBuffin_Filter = new Complex[FSK_Samplesize];
	memset(DataBuffin_Filter, 0x00, sizeof(Complex) * FSK_Samplesize);
	Complex* DataBuffout_Filter = new Complex[FSK_Samplesize];
	memset(DataBuffout_Filter, 0x00, sizeof(Complex) * FSK_Samplesize);


	for (int i = 0; i < FSK_Samplesize; i++)
	{
		DataBuffin_Filter[i].IData = fPLLPastFreqPart[i];
	}

	De_FSK.BlockFilter(DataBuffin_Filter, FSK_Samplesize, DataBuffout_Filter, cFilterBuff);

	//sprintf(file_path1, "D:/data/DataBuffout_Filter.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < FSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", DataBuffout_Filter[i].IData);
	//}
	//fclose(fp);

	//// 将采样率fs变换成4Rb,已进入符号同步
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

	// modify on 20240520
	Complex* DataSymbolSyncBuff = new Complex[FSK_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * FSK_Samplesize);
	double step = ((double)(4 * m_Rb)) / ((double)m_fs);
	int nSymbolSyncSize = 0;
	Downfs(DataBuffout_Filter, DataSymbolSyncBuff, step, FSK_Samplesize, &nSymbolSyncSize);

	//sprintf(file_path1, "D:/data/DataSymbolSyncBuff.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	//}
	//fclose(fp);

	Complex* DataSymbolSyncBuffout = new Complex[nSymbolSyncSize];
	memset(DataSymbolSyncBuffout, 0x00, sizeof(Complex) * nSymbolSyncSize);
	int nSymbolSync_count = 0;
	De_FSK.SymbolSync(DataSymbolSyncBuff, DataSymbolSyncBuffout, nSymbolSyncSize, &nSymbolSync_count, &SySyncBuffer);

	//sprintf(file_path1, "D:/data/DataSymbolSyncBuffout.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSync_count; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuffout[i].IData);
	//}
	//fclose(fp);

	// 画星座图

	// 调用星座图绘制接口
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

	// 阈值搜索
	float* DataI_sort = new float[nSymbolSync_count];
	memset(DataI_sort, 0x00, sizeof(float) * nSymbolSync_count);
	float* xielv = new float[nSymbolSync_count];
	memset(xielv, 0x00, sizeof(float) * nSymbolSync_count);

	for (int i = 0; i < nSymbolSync_count; i++)
	{
		DataI_sort[i] = DataSymbolSyncBuffout[i].IData;
	}

	De_FSK.Quick_Sort(DataI_sort, 0, nSymbolSync_count - 1);

	int aaa = 1;
	int count = 0;
	for (int i = 0 + aaa; i < nSymbolSync_count; i++)
	{
		xielv[i - aaa] = (DataI_sort[i] - DataI_sort[i - aaa]) / aaa; //差分
		count++;
	}

	int xielv_sort_end = nSymbolSync_count - 5;

	vector<float> xielv_sort;
	xielv_sort.resize(xielv_sort_end);   // 分配空间
	copy(xielv, xielv + xielv_sort_end, xielv_sort.begin());
	vector<size_t>idx;
	idx = sort_indexes_(xielv_sort);   // 从小到大排序，返回排序后元素在xielv_sort的index，xielv_sort数组不变
	int xielv_max1 = *(idx.end() - 1);

	int loc = 0;
	int xielv_max2 = 0;
	int th = 40; // modify on 20240525 将 th 从50设置到40
	while (++loc)
	{
		xielv_max2 = idx[idx.size() - 1 - loc];
		if (xielv_max1 - xielv_max2 > th || xielv_max1 - xielv_max2 < -th)
			break;
	}

	int xielv_max3 = 0;
	while (++loc)
	{
		xielv_max3 = idx[idx.size() - 1 - loc];
		if ((xielv_max2 - xielv_max3 > th || xielv_max2 - xielv_max3 < -th))
			break;
	}

	int xielv_max4 = 0;
	while (++loc)
	{
		xielv_max4 = idx[idx.size() - 1 - loc];
		if (
			(xielv_max3 - xielv_max4 > th || xielv_max3 - xielv_max4 < -th)
			&& (xielv_max2 - xielv_max4 > th || xielv_max2 - xielv_max4 < -th)
			&& (xielv_max1 - xielv_max4 > th || xielv_max1 - xielv_max4 < -th)
			)
			break;
	}

	int xielv_max5 = 0;
	while (++loc)
	{
		xielv_max5 = idx[idx.size() - 1 - loc];
		if (
			(xielv_max4 - xielv_max5 > th || xielv_max4 - xielv_max5 < -th)
			&& (xielv_max1 - xielv_max5 > th || xielv_max1 - xielv_max5 < -th)
			&& (xielv_max2 - xielv_max5 > th || xielv_max2 - xielv_max5 < -th)
			&& (xielv_max3 - xielv_max5 > th || xielv_max3 - xielv_max5 < -th)
			)
			break;
	}


	int xielv_max6 = 0;
	while (++loc)
	{
		xielv_max6 = idx[idx.size() - 1 - loc];
		if (
			(xielv_max5 - xielv_max6 > th || xielv_max5 - xielv_max6 < -th)
			&& (xielv_max1 - xielv_max6 > th || xielv_max1 - xielv_max6 < -th)
			&& (xielv_max2 - xielv_max6 > th || xielv_max2 - xielv_max6 < -th)
			&& (xielv_max3 - xielv_max6 > th || xielv_max3 - xielv_max6 < -th)
			&& (xielv_max4 - xielv_max6 > th || xielv_max4 - xielv_max6 < -th)
			)
			break;
	}

	int xielv_max7 = 0;
	while (++loc)
	{
		xielv_max7 = idx[idx.size() - 1 - loc];
		if (
			(xielv_max6 - xielv_max7 > th || xielv_max6 - xielv_max7 < -th)
			&& (xielv_max1 - xielv_max7 > th || xielv_max1 - xielv_max7 < -th)
			&& (xielv_max2 - xielv_max7 > th || xielv_max2 - xielv_max7 < -th)
			&& (xielv_max3 - xielv_max7 > th || xielv_max3 - xielv_max7 < -th)
			&& (xielv_max4 - xielv_max7 > th || xielv_max4 - xielv_max7 < -th)
			&& (xielv_max5 - xielv_max7 > th || xielv_max5 - xielv_max7 < -th)
			)
			break;
	}

	int xielc_index[7] = { xielv_max1,xielv_max2,xielv_max3,xielv_max4,xielv_max5,xielv_max6,xielv_max7 };
	sort(xielc_index, xielc_index + 7); //默认升序排序

	float threshold_test[7] = { 0 };
	float threshold_final[7] = { 0 };
	for (int i = 0; i < 7; i++)
	{
		threshold_test[i] = DataI_sort[xielc_index[i]]; // ？？？？？？
		threshold_final[7 - i - 1] = threshold_test[i];
	}

	char* c_symbol = new char[nSymbolSync_count];
	Judgment_8fsk(DataSymbolSyncBuffout, nSymbolSync_count, c_symbol, threshold_final);


	//sprintf(file_path1, "D:/data/FSK8_output_symbol.txt");
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
	DELETE_ARR(DataBuff_normalize);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(fPLLPastFreqPart);
	DELETE_ARR(fPLLPastFreqPart_abs);
	DELETE_ARR(DataBuffin_Filter);
	DELETE_ARR(DataBuffout_Filter);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(DataSymbolSyncBuffout);
	DELETE_ARR(DataI_sort);
	DELETE_ARR(xielv);
	//DELETE_ARR(c_symbol);

}

float Demo_8FSK::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

float Demo_8FSK::max_function(float* data, int start, int end)
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

int Demo_8FSK::max_index_function(float* data, int start, int end)
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

void Demo_8FSK::Judgment_8fsk(Complex buffin[], int bufflen, char buffout[], float* threshold)
{
	float th1 = threshold[0];
	float th2 = threshold[1];
	float th3 = threshold[2];
	float th4 = threshold[3];
	float th5 = threshold[4];
	float th6 = threshold[5];
	float th7 = threshold[6];

	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData >= th1)
		{
			buffout[i] = 7;
		}
		else if (buffin[i].IData < th1 && buffin[i].IData >= th2)
		{
			buffout[i] = 6;
		}
		else if (buffin[i].IData < th2 && buffin[i].IData >= th3)
		{
			buffout[i] = 5;
		}
		else if (buffin[i].IData < th3 && buffin[i].IData >= th4)
		{
			buffout[i] = 4;
		}
		else if (buffin[i].IData < th4 && buffin[i].IData >= th5)
		{
			buffout[i] = 3;
		}
		else if (buffin[i].IData < th5 && buffin[i].IData >= th6)
		{
			buffout[i] = 2;
		}
		else if (buffin[i].IData < th6 && buffin[i].IData >= th7)
		{
			buffout[i] = 1;
		}
		else if (buffin[i].IData < th7)
		{
			buffout[i] = 0;
		}
	}

}

void Demo_8FSK::InitBlockCFilter()
{
	if (m_Rb <= 4e5)
	{
		// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员

		CFilter.nFilterTaps = 21;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0268661128107165 ,-0.0272597932734082, -0.0134774114144576,	0.0182063277002322,	0.0684930940231012,	0.134330564053582,
			0.208991748429463,	0.283025639703609,	0.345920226304412,	0.388127532797573,	0.402991692160747,	0.388127532797573,
			0.345920226304412,	0.283025639703609,	0.208991748429463,	0.134330564053582,	0.0684930940231012,	0.0182063277002322,
			-0.0134774114144576, -0.0272597932734082, -0.0268661128107165
		};//8FSK

		// 将初始化值复制到动态分配的数组中
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{
		CFilter.nFilterTaps = 17;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		float tempCoef[] =
		{
			-0.0300349916229418, -0.0283172616659510,	1.03450095994087e-17,	0.0606798464270379,	0.150174958114709,	0.254855354993559,
			0.353841408874743,	0.424758924989266,	0.450524874344127,	0.424758924989266,	0.353841408874743,	0.254855354993559,
			0.150174958114709,	0.0606798464270379,	1.03450095994087e-17, -0.0283172616659510, -0.0300349916229418

		};//8FSK

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

void Demo_8FSK::InitialDemodulation(const DemodulationInitParamater& signalinit)
{

	m_fs = signalinit.fs;
	m_Rb = signalinit.rb;

	//for (int i = 0; i < 8; i++) 
	//{
	//	m_fc[i] = signalinit.fc_demo[i];
	//}
	// modify on 20240520
	memcpy(m_fc, signalinit.fc_demo, sizeof(float) * 8);

	InitBlockCFilter();

	FSK8Init();

	InitCostasPLL();

	InitSymbolSync();
}

void Demo_8FSK::FSK8Init()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_FSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;													//降采样时的时钟控制  
	fDownConversionPhase = 0;											//下变频中的相位值
	nPLLBuffSize = 0;

	// add on 20240520
	flag = 0;
	slice = 45;
	FSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[FSK_Samplesize];

	//******************符号同步初始化***********
	SySyncBuffer.CSymbolSyncBuff = new Complex[4];
	SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBuffer.fSymbolSyncW = 0.5;									//符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBuffer.fSymbolSyncN = 0.9;									//符号同步NCO寄存器，初值设为1
	SySyncBuffer.fSymbolSyncNTemp = 0.9;								//符号同步NCO暂时的寄存器，初值设为1
	SySyncBuffer.fSymbolSyncU = 0.6;									//符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBuffer.nSymbolSyncKK = 0;										//符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBuffer.fSymbolSyncPastW = 0;									//符号同步W的缓存值
	SySyncBuffer.fSymbolSyncPastN = 0;									//符号同步N的缓存值
	SySyncBuffer.fSymbolSyncTimeError1 = 0;								//符号同步time_error缓存值
	SySyncBuffer.fSymbolSyncTimeError2 = 0;								//符号同步time_error缓存值
}

void Demo_8FSK::InitSymbolSync()
{

	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.001999798626123;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 2.000201333325170e-06;	// 8FSK

	De_FSK.m_sAlgDemInit.nSampPerSymb = 4;								//每个码元中的采样点数

	m_sAlgDemInit.nSampPerSymb = 4;										//每个码元中的采样点数

	De_FSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_8FSK::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;												//环路阻尼系数
	int Ko = 1;															//压控振荡器增益
	int Kd = 1;															//鉴相器增益
	float fBLcoef1 = 0.33;												//8FSK 
	float fBLcoef2 = 0.28;												//8FSK 
	float fBLcoef = 0;													//8FSK 
	int K = Ko * Kd;

	fBLcoef = (m_Rb <= 4e5 ? fBLcoef1 : fBLcoef2);

	if (m_fc[0] > 0 && m_fc[1] < 0 || m_fc[0] < 0 && m_fc[1]>0)
	{

		BL = fBLcoef * (m_Rb + fabs(m_fc[1]) + fabs(m_fc[0]));
	}
	else
	{
		BL = fBLcoef * (m_Rb + fabs(m_fc[1] - m_fc[0]));
	}

	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = (float)1 / m_fs;

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_FSK.CostasPll = CostasPll;

	fPLLPastFreqPartemp = 0;										   //锁相环中的频率跟踪中间变量
	fPLLNCO = 0;														//锁相环中的本地NCO
	//fPLLPastFreqPart = 0;												//锁相环中的频率跟踪曲线
}

// add on 20240520
void Demo_8FSK::Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout)
{
	int i = 0;
	int j = 0;
	int k = 0;
	static double step_temp = 0;
	while (1)
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



void Demo_8FSK::AGCFSK8(Complex* BuffIn, float target_power, int len)
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


float Demo_8FSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}