#include "SignalPM_BP.h"

Demo_PM_BP::Demo_PM_BP(int index):m_index(index)
{

}

Demo_PM_BP::~Demo_PM_BP() 
{

}

bool Demo_PM_BP::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    m_demodulationInitParamater = demoParameter;
    subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter);
    return 1;
}

void Demo_PM_BP::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
    Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_PM_BP::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_PM_BP_Out)
{

    
	//FILE*fp;
	//char file_path1[500];

	//FILE*fp2;
	//char file_path2[500];


    Complex* DataBuff = new Complex[m_SamleSize];
    memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);
    for (int m = 0; m < m_SamleSize; m++)
    {
        DataBuff[m].IData = data_input_slice[m].IData;
        DataBuff[m].QData = data_input_slice[m].QData;
    }
	PM_BP.AGC(DataBuff, fAGCPastVc);
	//sprintf(file_path1, "D:/data/intputI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/intputQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataBuff[i].QData);
	//}
	//fclose(fp);

    Complex* DataPLLBuff_PM = new Complex[m_SamleSize];
    memset(DataPLLBuff_PM, 0x00, sizeof(Complex) * m_SamleSize);
    float* PLL_Freq_Part_2 = new float[m_SamleSize];
    memset(PLL_Freq_Part_2, 0x00, sizeof(float) * m_SamleSize);
    float* PLL_Phase_Part_2 = new float[m_SamleSize];
    memset(PLL_Phase_Part_2, 0x00, sizeof(float) * m_SamleSize);
    PLLCostas(DataBuff, DataPLLBuff_PM, m_SamleSize, &fPLLNCO,&fPLLPastFreqPart, PLL_Freq_Part_2, PLL_Phase_Part_2);
    
	//sprintf(file_path1, "D:/data/PLL_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", PLL_Phase_Part_2[i]);
	//}
	//fclose(fp);

    Complex* s_r1 = new Complex[m_SamleSize];
    memset(s_r1, 0x00, sizeof(Complex) * m_SamleSize);

	fc2 = -450000;
    downphase(PLL_Phase_Part_2, s_r1, m_SamleSize,&fc2,&fOrthDownConversionPhase); 
    
	//sprintf(file_path1, "D:/data/PLL_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", s_r1[i].IData);
	//}
	//fclose(fp);

    Complex* s_1 = new Complex[m_SamleSize];
    memset(s_1, 0x00, sizeof(Complex) * m_SamleSize); 
    PM_BP.BlockFilter(s_r1, m_SamleSize,s_1,cFilterBuff);
    

	//sprintf(file_path1, "D:/data/PLL_I1.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", s_1[i].IData);
	//}
	//fclose(fp);

	//FILE*fp;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/PLL_I.txt");
	//fp = fopen(file_path2, "at");
	//for (int i = 0; i < m_SamleSize; i++) {
	//	fprintf(fp, "%f\n", s_1[i].IData);
	//}
	//fclose(fp);

	//FILE*fp1;
	//char file_path21[500];
	//sprintf(file_path21, "D:/data/PLL_Q.txt");
	//fp1 = fopen(file_path21, "at");
	//for (int i = 0; i < m_SamleSize; i++) {
	//	fprintf(fp1, "%f\n", s_1[i].QData);
	//}
	//fclose(fp1);


    Complex* s_11 = new Complex[m_SamleSize/2];
    memset(s_11, 0x00, sizeof(Complex) * (m_SamleSize/2));
    for (int m = 0; m < m_SamleSize; m=m+2)
    {
        s_11[m / 2].IData = s_1[m].IData;
        s_11[m / 2].QData = s_1[m].QData;
    }

    float* Signal_sort = new float[m_SamleSize/2];
    memset(Signal_sort, 0x00, sizeof(float) * (m_SamleSize / 2));
    for (int i = 0; i < m_SamleSize / 2; i++)
    {
        Signal_sort[i] = fabs(s_11[i].IData);
    }
    PM_BP.Quick_Sort(Signal_sort, 0, m_SamleSize/2-1);
    float mean_power = mean_function(Signal_sort, m_SamleSize/2-70, m_SamleSize/2-40);
    Complex* Signal_Channel_3 = new Complex[m_SamleSize / 2];
    memset(Signal_Channel_3, 0x00, sizeof(Complex) * m_SamleSize / 2);
    for (int i = 0; i < m_SamleSize / 2; i++)
    {
        Signal_Channel_3[i].IData = s_11[i].IData / mean_power;
        Signal_Channel_3[i].QData = s_11[i].QData / mean_power;
    }

    
    Complex* DataPLLBuff_BP = new Complex[m_SamleSize/2];
    memset(DataPLLBuff_BP, 0x00, sizeof(Complex) * m_SamleSize/2);
    Complex* DataFLLBuff_BP = new Complex[m_SamleSize/2];
    memset(DataFLLBuff_BP, 0x00, sizeof(Complex) * m_SamleSize/2);
    PM_BP.CarrierSync(Signal_Channel_3, DataPLLBuff_BP, DataFLLBuff_BP, m_SamleSize/2, 1, &fPLLNCO_1, &fPLLPastFreqPart_1, &FLLBuffer);

	//FILE*fp;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/PLL_I.txt");
	//fp = fopen(file_path2, "at");
	//for (int i = 0; i < m_SamleSize; i++) {
	//	fprintf(fp, "%f\n", DataPLLBuff_BP[i].IData);
	//}
	//fclose(fp);

	//FILE*fp1;
	//char file_path21[500];
	//sprintf(file_path21, "D:/data/PLL_Q.txt");
	//fp1 = fopen(file_path21, "at");
	//for (int i = 0; i < m_SamleSize; i++) {
	//	fprintf(fp1, "%f\n", DataPLLBuff_BP[i].QData);
	//}
	//fclose(fp1);

    int nSymbolSyncSize = 0;
    Complex* DataSymbolSyncBuff = new Complex[m_SamleSize/2];
    memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize/2);
    PM_BP.SymbolSync(DataPLLBuff_BP, DataSymbolSyncBuff, m_SamleSize/2, &nSymbolSyncSize, &SySyncBuffer);

    
    char* c_symbol = new char[nSymbolSyncSize];
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbol);

	for (int i = 0; i < nSymbolSyncSize; i++)
	{
		DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].IData;
	}
	//sprintf(file_path1, "D:/data/PM_BPSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", c_symbol[i]);
	//}
	//fclose(fp);
	signal_PM_BP_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_PM_BP_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	signal_PM_BP_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
   
    signal_PM_BP_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
    memcpy(signal_PM_BP_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);

    DELETE_ARR(DataBuff);
    DELETE_ARR(DataPLLBuff_PM);
    DELETE_ARR(PLL_Freq_Part_2);
    DELETE_ARR(PLL_Phase_Part_2);
    DELETE_ARR(s_r1);
    DELETE_ARR(s_1);
    DELETE_ARR(s_11);
    DELETE_ARR(Signal_sort);
    DELETE_ARR(Signal_Channel_3);
    DELETE_ARR(DataPLLBuff_BP);
    DELETE_ARR(DataFLLBuff_BP);
    DELETE_ARR(DataSymbolSyncBuff);
    DELETE_ARR(c_symbol);

}
//全部初始化
void Demo_PM_BP::InitialDemodulation(const DemodulationInitParamater& info) {

    m_Rb = info.rb;
    fc2 = info.fc_demo[0];

    InitBlockCFilter();
    
    PM_BPInit();
    
    InitCostasPLL();
    
    InitSymbolSync();
    
}

void Demo_PM_BP::InitBlockCFilter() {
    // 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
    CFilter.nFilterTaps = 66;
    CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

    // 使用数组初始化的值
    float tempCoef[] = { 0.000300236790815230,- 0.000318083254178556,- 0.000855103154930395,- 0.000990268382954444,- 0.000487712148667362,
        0.000588523122424037,0.00172482575409843,0.00209416001114999,0.00104957199148652,- 0.00126316787160771,- 0.00364730576829181,
        - 0.00433407433765904,- 0.00211987302550528,0.00248853880247738,0.00701524733957493,0.00815335440937972,0.00390961798577038,
        - 0.00451151788578127,- 0.0125389944244636,- 0.0144144829075778,- 0.00686081390653839,0.00788960724885472,0.0219506333546247,
        0.0253950771878517,0.0122440038486131,- 0.0143814602156421,- 0.0413223332005430,- 0.0501396000048169,- 0.0259436042862064,
        0.0339401732374600,0.116211391424138,0.195362079644283,0.243811352622362,0.243811352622362,0.195362079644283,0.116211391424138,
        0.0339401732374600,- 0.0259436042862064,- 0.0501396000048169,- 0.0413223332005430,- 0.0143814602156421,0.0122440038486131,
        0.0253950771878517,0.0219506333546247,0.00788960724885472,- 0.00686081390653839,- 0.0144144829075778,- 0.0125389944244636,
        - 0.00451151788578127,0.00390961798577038,0.00815335440937972,0.00701524733957493,0.00248853880247738,- 0.00211987302550528,
        - 0.00433407433765904,- 0.00364730576829181,- 0.00126316787160771,0.00104957199148652,0.00209416001114999,0.00172482575409843,
        0.000588523122424037,- 0.000487712148667362,- 0.000990268382954444,- 0.000855103154930395,- 0.000318083254178556,0.000300236790815230 };
    // 将初始化值复制到动态分配的数组中
    for (int i = 0; i < CFilter.nFilterTaps; ++i) {
        CFilter.fFilterCoef[i] = tempCoef[i];
    }
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    PM_BP.CFilter = CFilter;
}

void Demo_PM_BP::PM_BPInit() 
{
    m_sAlgDemInit.llSliceLength = 8192;
    m_SamleSize = m_sAlgDemInit.llSliceLength;
    PM_BP.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;
    m_sAlgDemInit.nSigFs = 8 * m_Rb;
    nDownSampClock = 1;                                         //降采样时的时钟控制  
    fDownConversionPhase = 0;                                   //下变频中的相位值
    fPLLNCO = 0;                                                //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                       //锁相环中的频率跟踪曲线
    fPLLNCO_1 = 0;                                              //锁相环中的本地NCO
    fPLLPastFreqPart_1 = 0;                                     //锁相环中的频率跟踪曲线
    fOrthDownConversionPhase = 0;
    nPLLBuffSize = 0;
    //******************符号同步初始化***********
    SySyncBuffer.CSymbolSyncBuff = new Complex[4];
    SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
    memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
    memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
    SySyncBuffer.fSymbolSyncW = 0.5;                            //符号同步环路滤波器输出寄存器，初值设为0.5
    SySyncBuffer.fSymbolSyncN = 1;                              //符号同步NCO寄存器，初值设为1
    SySyncBuffer.fSymbolSyncNTemp = 1;                          //符号同步NCO暂时的寄存器，初值设为1
    SySyncBuffer.fSymbolSyncU = 0.6;                            //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
    SySyncBuffer.nSymbolSyncKK = 0;                             //符号同步用来表示Ti时间序号,指示u,yI,yQ
    SySyncBuffer.fSymbolSyncPastW = 0;                          //符号同步W的缓存值
    SySyncBuffer.fSymbolSyncPastN = 0;                          //符号同步N的缓存值
    SySyncBuffer.fSymbolSyncTimeError1 = 0;                     //符号同步time_error缓存值
    SySyncBuffer.fSymbolSyncTimeError2 = 0;                     //符号同步time_error缓存值
}

void Demo_PM_BP::InitSymbolSync()
{
    symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823/20;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05/400;
    PM_BP.m_sAlgDemInit.nSampPerSymb = 4;                       //每个码元中的采样点数

    PM_BP.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_PM_BP::InitCostasPLL() {
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;     //环路阻尼系数
    int Ko = 1;              //压控振荡器增益
    int Kd = 1;              //鉴相器增益
    float fBLcoef = 0.1;
    int K = Ko * Kd;
    BL = fBLcoef * m_Rb;
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (8 * (float)m_Rb);
    CostasPll.fPLLLoopFilterCoef1 = 2* sigma* Wn* T_nco;
    CostasPll.fPLLLoopFilterCoef2 = Wn * T_nco* Wn * T_nco;
    PM_BP.CostasPll = CostasPll;
}

void Demo_PM_BP::Judgment(Complex buffin[], int bufflen, char buffout[])
{
    for (int i=0;i< bufflen;i++)
    {
        if (buffin[i].IData>0)
        {
            buffout[i] = 1;
        }
        else
        {
            buffout[i] = 0;
        }
    }
}

void Demo_PM_BP::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart,float* PLL_Freq_Part, float* PLL_Phase_Part)
{
    float discriminator_out = 0;
    float pll_phase_part;
    float freq_control;
    float past_pll_freq_part;
    float past_nco_phase;
    float PLLCoef1 = 0.004011138433851;
    float PLLCoef2 = 7.047045975642138e-06;
    for (int i = 0; i < nBuffSize; i++)
    {
        COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fPLLNCO), sin(*fPLLNCO), BuffOut[i].IData, BuffOut[i].QData);
        discriminator_out = atan2(BuffOut[i].QData, BuffOut[i].IData);
        pll_phase_part = discriminator_out * PLLCoef1;
        PLL_Phase_Part[i] = pll_phase_part;

        past_pll_freq_part = *fPLLPastFreqPart;
        *fPLLPastFreqPart = discriminator_out * PLLCoef2 + past_pll_freq_part;
        PLL_Freq_Part[i] = *fPLLPastFreqPart;
        freq_control = pll_phase_part + *fPLLPastFreqPart;
        past_nco_phase = *fPLLNCO;
        *fPLLNCO = past_nco_phase + freq_control * 2 * Pi;    
    }
}

void Demo_PM_BP::downphase(float BuffIn[], Complex BuffOut[], int nBuffSize,float* fDownConversionFc, double* fOrthDownConversionPhase)
{
    float phase0 = 2 * Pi * (*fDownConversionFc) / m_sAlgDemInit.nSigFs;
    for (int i=0;i< nBuffSize;i++)
    {
        *fOrthDownConversionPhase += phase0;
        if (*fOrthDownConversionPhase>=2*Pi)
        {
            *fOrthDownConversionPhase -= 2 * Pi;
        }
        BuffOut[i].IData = BuffIn[i] *  cos(*fOrthDownConversionPhase);
        BuffOut[i].QData = BuffIn[i] * -sin(*fOrthDownConversionPhase);
    }

}

float Demo_PM_BP::mean_function(float* data, int start, int end)
{
    float sum = 0;
    int i = 0;
    for (i = start; i < end; i++)
    {
        sum += data[i];
    }
    return sum / (end - start);
}
