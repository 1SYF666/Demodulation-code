#include "Signal16APSK.h"

Demo_16APSK::Demo_16APSK(int index):m_index(index)
{

}

Demo_16APSK::~Demo_16APSK() {

}

bool Demo_16APSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    m_demodulationInitParamater = demoParameter;
    subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb); return 1;
}

void Demo_16APSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_16APSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out) {

	//FILE*fp;
	//char file_path1[500];
	//sprintf(file_path1, "D:/data/input_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < 8192*5; i++) {
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//FILE*fp2;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/input_Q.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < 8192*5; i++) {
	//	fprintf(fp2, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp2);

	Complex* DataBuff = new Complex[APSK_Samplesize];
	memset(DataBuff, 0x00, sizeof(Complex) * APSK_Samplesize);

	De_APSK.BlockFilter(data_input_slice, APSK_Samplesize, DataBuff, cFilterBuff);


	De_APSK.AGC(DataBuff, fAGCPastVc);
	
	//sprintf(file_path1, "D:/data/DataBuff_normalizeI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < 8192*5; i++) {
	//	fprintf(fp, "%f\n", DataBuff[i].IData);
	//}
	//fclose(fp);


	//sprintf(file_path2, "D:/data/DataBuff_normalizeQ.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < 8192*5; i++)
	//{
	//	fprintf(fp2, "%f\n", DataBuff[i].QData);
	//}
	//fclose(fp2);

    Complex* DataPLLBuff = new Complex[APSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * APSK_Samplesize);
	Complex* DataFLLBuff = new Complex[APSK_Samplesize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * APSK_Samplesize);

	De_APSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, APSK_Samplesize, 6, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	//sprintf(file_path1, "D:/data/PLLI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < APSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", DataPLLBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/PLLQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < APSK_Samplesize; i++)
	//{
	//	fprintf(fp, "%f\n", DataPLLBuff[i].QData);
	//}
	//fclose(fp);

    int nSymbolSyncSize = 0;                                        //突发在当前片中的长度
	Complex* DataSymbolSyncBuff = new Complex[APSK_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * APSK_Samplesize);

	De_APSK.SymbolSync(DataPLLBuff, DataSymbolSyncBuff, APSK_Samplesize, &nSymbolSyncSize, &SySyncBuffer);

	//signal_16APSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	//memcpy(signal_16APSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	//signal_16APSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
	//sprintf(file_path1, "D:/data/DataSymbolSyncBuffoutI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/DataSymbolSyncBuffoutQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].QData);
	//}
	//fclose(fp);

	float* Signal_sort = new float[nSymbolSyncSize];
	memset(Signal_sort, 0x00, sizeof(float) * nSymbolSyncSize);
	for (int i = 0; i < nSymbolSyncSize; i++) {
		Signal_sort[i] = DataSymbolSyncBuff[i].IData;  
	}
	De_APSK.Quick_Sort(Signal_sort, 0, nSymbolSyncSize - 1);
	float mean_power = mean_function(Signal_sort,nSymbolSyncSize / 4 * 8 / 10, nSymbolSyncSize / 4 * 9 / 10);
	for (int i = 0; i < nSymbolSyncSize; i++) {
		DataSymbolSyncBuff[i].IData = DataSymbolSyncBuff[i].IData/ mean_power;
		DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].QData/ mean_power;
	}
    Complex* DataSymbol_m = new Complex[nSymbolSyncSize];
    for (int i = 0; i < nSymbolSyncSize; i++) {
        DataSymbol_m[i].IData = DataSymbolSyncBuff[i].IData;
        DataSymbol_m[i].QData = DataSymbolSyncBuff[i].QData;
    }

    
    Complex* Data_phase= new Complex[nSymbolSyncSize];
    Phase_correction(DataSymbol_m, Data_phase, nSymbolSyncSize);

	//sprintf(file_path1, "D:/data/DataSymbolSyncBuffoutI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", Data_phase[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/DataSymbolSyncBuffoutQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%f\n", Data_phase[i].QData);
	//}
	//fclose(fp);

	signal_16APSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_16APSK_Out->burstData->softDistinguishData, Data_phase, (sizeof(Complex))*nSymbolSyncSize);
	signal_16APSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	float* Data_power = new float[nSymbolSyncSize];
	Judgth1 = GetThreshold(Data_phase, nSymbolSyncSize, Data_power);
    
    
    int* c_symbol = new int[nSymbolSyncSize];
    Judgment(Data_phase, Data_power, nSymbolSyncSize, c_symbol);

	//sprintf(file_path1, "D:/data/csymbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", c_symbol[i]);
	//}
	//fclose(fp);

    char* DataResult = new char[nSymbolSyncSize];
	diff_code_16apsk(c_symbol, DataResult, nSymbolSyncSize);

	signal_16APSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_16APSK_Out->burstData->nDemodulationByteI, DataResult, (sizeof(char))*nSymbolSyncSize);

	//sprintf(file_path1, "D:/data/APSK16_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);

	nSilceCount++;

	DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataSymbolSyncBuff);
    DELETE_ARR(DataResult);
    DELETE_ARR(c_symbol);
    DELETE_ARR(Data_phase);
    DELETE_ARR(DataSymbol_m);

    // add on 20240627
    DELETE_ARR(Signal_sort);
    DELETE_ARR(Data_power);
    

}

void Demo_16APSK::InitBlockCFilter()
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
    De_APSK.CFilter = CFilter;
}

void Demo_16APSK::InitialDemodulation(int Rb) 
{

    m_Rb = Rb;
	InitBlockCFilter();
    APSKInit();
    InitCostasPLL();
    InitSymbolSync();
}

void Demo_16APSK::APSKInit()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	m_sAlgDemInit.llSliceLength = 8192*5;
    De_APSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	//add on 20240514
	flag = 0;
	slice = 5;
	APSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[APSK_Samplesize];

    nDownSampClock = 1;                                     //降采样时的时钟控制  
    fDownConversionPhase = 0;                               //下变频中的相位值
    fPLLNCO = 0;                                            //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                   //锁相环中的频率跟踪曲线
    nPLLBuffSize = 0;
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
}

void Demo_16APSK::InitSymbolSync() {
    symbolsyncsactorInit.fSymbolSyncFactor1 = 0.001666498855103;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 1.389028703698034e-06;
    De_APSK.m_sAlgDemInit.nSampPerSymb = 4;                 //每个码元中的采样点数

    De_APSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_16APSK::InitCostasPLL() {
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;                                    //环路阻尼系数
    int Ko = 1;                                             //压控振荡器增益
    int Kd = 1;                                             //鉴相器增益
    float fBLcoef = 0.0020;
    int K = Ko * Kd;
    BL = fBLcoef * m_Rb;
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (4 * (float)m_Rb);

    CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
    De_APSK.CostasPll = CostasPll;
}

void Demo_16APSK::Judgment(Complex buffin[], float power[], int bufflen, int buffout[]) {

    float theta = 0 ;
    for (int i = 0; i < bufflen; i++)
    {
        theta = atan2(buffin[i].QData, buffin[i].IData);
        if (theta > 0)
        {
            theta = theta;
        }
        else
        {
            theta = theta + 2 * Pi;
        }
        if (power[i] > Judgth1)
        {
            if (theta >= 0 && theta < 2 * Pi / 12)
            {
                buffout[i] = 13;
            }
            else if (theta >= 2 * Pi / 12 && theta < 4 * Pi / 12)
            {
                buffout[i] = 14;
            }
            else if (theta >= 4 * Pi / 12 && theta < 6 * Pi / 12)
            {
                buffout[i] = 15;
            }
            else if (theta >= 6 * Pi / 12 && theta < 8 * Pi / 12)
            {
                buffout[i] = 9;
            }
            else if (theta >= 8 * Pi / 12 && theta < 10 * Pi / 12)
            {
                buffout[i] = 10;
            }
            else if (theta >= 10 * Pi / 12 && theta < 12 * Pi / 12)
            {
                buffout[i] = 11;
            }
            else if (theta >= 12 * Pi / 12 && theta < 14 * Pi / 12)
            {
                buffout[i] = 1;
            }
            else if (theta >= 14 * Pi / 12 && theta < 16 * Pi / 12)
            {
                buffout[i] = 2;
            }
            else if (theta >= 16 * Pi / 12 && theta < 18 * Pi / 12)
            {
                buffout[i] = 3;
            }
            else if (theta >= 18 * Pi / 12 && theta < 20 * Pi / 12)
            {
                buffout[i] = 5;
            }
            else if (theta >= 20 * Pi / 12 && theta < 22 * Pi / 12)
            {
                buffout[i] = 6;
            }
            else if (theta >= 22 * Pi / 12 && theta <= 24.1 * Pi / 12)
            {
                buffout[i] = 7;
            }
        }
        else
        {
            if (theta >= 0 && theta < Pi / 2)
            {
                buffout[i] = 12;
            }
            else if (theta >= Pi / 2 && theta < 2 * Pi / 2)
            {
                buffout[i] = 8;
            }
            else if (theta >= 2 * Pi / 2 && theta < 3 * Pi / 2)
            {
                buffout[i] = 0;
            }
            else
            {
                buffout[i] = 4;
            }
        }
    }
}

float Demo_16APSK::mean_function(float* data, int start, int end) {
	float sum = 0;
	int i;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_16APSK::diff_code_16apsk(int* y, char* c_symbol,int len) {
    int c = 0, i = 0, j = 0;
    int n = 2;

    int** a = new int* [len];
    for (i = 0; i < len; ++i) {
        a[i] = new int[4];
    }

    for (i = 0; i < len; i++) {
        j = 0;
        while (y[i] > 0) {
            c = (y[i] % n);
            a[i][j] = c;
            y[i] = y[i] / n;
            j++;
        }
        if (j == 0) {
            a[i][0] = 0; a[i][1] = 0; a[i][2] = 0; a[i][3] = 0;
        }
        else if (j == 1) {
            a[i][1] = 0; a[i][2] = 0; a[i][3] = 0;
        }
        else if (j == 2) {
            a[i][2] = 0; a[i][3] = 0;
        }
        else if (j == 3) {
            a[i][3] = 0;
        }
    }

    int* c_symbol_2 = new int[len];
    memset(c_symbol_2, 0, sizeof(int) * len);

    for (i = 0; i < len; i++) {
        if (a[i][3] == 0 && a[i][2] == 0) {
            c_symbol_2[i] = 0;
        }
        else if (a[i][3] == 1 && a[i][2] == 0) {
            c_symbol_2[i] = 1;
        }
        else if (a[i][3] == 1 && a[i][2] == 1) {
            c_symbol_2[i] = 2;
        }
        else if (a[i][3] == 0 && a[i][2] == 1) {
            c_symbol_2[i] = 3;
        }
    }

    int** d_bit = new int* [len];
    for (i = 0; i < len; ++i) {
        d_bit[i] = new int[4];
    }

    for (i = 0; i < len; ++i) {
        d_bit[i][3] = a[i][0];
        d_bit[i][2] = a[i][1];
    }

    int d_bit_2;
    for (i = 0; i < len; i++) {
        if (i == 0) {
            d_bit_2 = ((c_symbol_2[i] - difftemp) + 4) % 4;
        }
        else {
            d_bit_2 = ((c_symbol_2[i] - c_symbol_2[i - 1]) + 4) % 4;
        }
        if (d_bit_2 == 0) {
            d_bit[i][0] = 0;
            d_bit[i][1] = 0;
        }
        else if (d_bit_2 == 1) {
            d_bit[i][0] = 1;
            d_bit[i][1] = 0;
        }
        else if (d_bit_2 == 2) {
            d_bit[i][0] = 1;
            d_bit[i][1] = 1;
        }
        else if (d_bit_2 == 3) {
            d_bit[i][0] = 0;
            d_bit[i][1] = 1;
        }
    }
	difftemp = c_symbol_2[len - 1];
    int* d_bit_1 = new int[len];
    memset(d_bit_1, 0, sizeof(int) * len);

    for (i = 0; i < len; i++) {
        d_bit_1[i] = d_bit[i][0] * 1000 + d_bit[i][1] * 100 + d_bit[i][2] * 10 + d_bit[i][3] * 1;
    }

    int remainder = 0;
    for (j = 0; j < len; j++) {
        i = 0; c_symbol[j] = 0;
        while (d_bit_1[j] != 0) {
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

void Demo_16APSK::Phase_correction(Complex Buffin[], Complex Buffout[], int bufflen) {
    float* theta_0 = new float[bufflen];
    memset(theta_0, 0x00, sizeof(float) * bufflen);
    int n1 = 0;
    for (int i = 0; i < bufflen; i++)
    {
        float power_s = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
        if (power_s < 0.8115)
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
    De_APSK.Quick_Sort(sort_theta_1, 0, n1 - 1);
    float real_data = 0, imag_data = 0;
	float x = mean_function(sort_theta_1, (n1 - 1) * 3 / 8, (n1 - 1) * 5 / 8);
	if (mean_sort_theta_1 != 0 && fabs(x) < fabs(mean_sort_theta_1)*1.05&&fabs(x) > fabs(mean_sort_theta_1)*0.95) {
		mean_sort_theta_1 = (x+ mean_sort_theta_1 * (nSilceCount - 1)) / nSilceCount;
	}
	else if(mean_sort_theta_1==0){
		mean_sort_theta_1 = x;
	}
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

float Demo_16APSK::GetThreshold(Complex Buffin[],int bufflen, float Buffout[])
{

    for (int i = 0; i < bufflen; i++)
    {
        Buffout[i] = sqrt(Buffin[i].IData * Buffin[i].IData + Buffin[i].QData * Buffin[i].QData);
    }
    float* sort_power = new float[bufflen];
    memset(sort_power, 0, bufflen);
	for (int i = 0; i < bufflen; i++) {
		sort_power[i] = Buffout[i];
	}
    De_APSK.Quick_Sort(sort_power, 0, bufflen - 1);

    int Window_width = 4;
    float* Slope = new float[bufflen];
    memset(Slope, 0, sizeof(float) * bufflen);
    for (int i = Window_width; i < bufflen; i++)
    {
        Slope[i - Window_width] = (sort_power[i] - sort_power[i - Window_width]) / Window_width;
    }

    float max1 = max_function(Slope, (bufflen - Window_width) * 1 / 16, (bufflen - Window_width) * 15 / 16);
    int loc1 = 0;
    for (int j = 0; j < (bufflen - Window_width + 1); j++)
    {
        if (Slope[j] == max1)
        {
            loc1 = j;
        }
    }
    float yuzhi = sort_power[loc1];

    DELETE_ARR(sort_power);
    DELETE_ARR(Slope);


    return  yuzhi;
}

float Demo_16APSK::max_function(float* data, int start, int end)
{
    float max1 = data[start];
    int i;
    for (i = start + 1; i < end; i++)
    {
        if (max1 < data[i])
        {
            max1 = data[i];
        }
    }
    return max1;
}
