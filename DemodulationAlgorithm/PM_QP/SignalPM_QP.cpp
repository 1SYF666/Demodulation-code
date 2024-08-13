#include "SignalPM_QP.h"

Demo_PM_QP::Demo_PM_QP(int index):m_index(index)
{

}

Demo_PM_QP::~Demo_PM_QP() 
{

}

bool Demo_PM_QP::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    m_demodulationInitParamater = demoParameter;
    subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter);
    return 1;
}

void Demo_PM_QP::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
    Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}



void Demo_PM_QP::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_PM_QP_Out) {



    //FILE* fp;
    //char file_path[500];


    Complex* DataBuff = new Complex[m_SamleSize];
    memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);
    for (int m = 0; m < m_SamleSize; m++)
    {
        DataBuff[m].IData = data_input_slice[m].IData;
        DataBuff[m].QData = data_input_slice[m].QData;
    }
    PM_QP.AGC(DataBuff, fAGCPastVc);
    
   
    
    Complex* DataPLLBuff_PM = new Complex[m_SamleSize];
    memset(DataPLLBuff_PM, 0x00, sizeof(Complex) * m_SamleSize);
    float* PLL_Freq_Part_2 = new float[m_SamleSize];
    memset(PLL_Freq_Part_2, 0x00, sizeof(float) * m_SamleSize);
    float* PLL_Phase_Part_2 = new float[m_SamleSize];
    memset(PLL_Phase_Part_2, 0x00, sizeof(float) * m_SamleSize);
    PLLCostas(DataBuff, DataPLLBuff_PM, m_SamleSize, &fPLLNCO, &fPLLPastFreqPart, PLL_Freq_Part_2, PLL_Phase_Part_2);
        
	fc2 = -450000;
    Complex* s_r1 = new Complex[m_SamleSize];
    memset(s_r1, 0x00, sizeof(Complex) * m_SamleSize);
    downphase(PLL_Phase_Part_2, s_r1, m_SamleSize, &fc2, &fOrthDownConversionPhase);
    

    Complex* s_1 = new Complex[m_SamleSize];
    memset(s_1, 0x00, sizeof(Complex) * m_SamleSize);
    PM_QP.BlockFilter(s_r1, m_SamleSize, s_1, cFilterBuff);
    


    Complex* s_11 = new Complex[m_SamleSize / 2];
    memset(s_11, 0x00, sizeof(Complex) * m_SamleSize / 2);
    for (int m = 0; m < m_SamleSize; m = m + 2)
    {
        s_11[m / 2].IData = s_1[m].IData;
        s_11[m / 2].QData = s_1[m].QData;
    }
    float* s_11_real = new float[m_SamleSize / 2];
    memset(s_11_real, 0x00, sizeof(float) * m_SamleSize / 2);
    for (int m = 0; m < m_SamleSize/2; m++)
    {
        s_11_real[m] = s_11[m].IData;
    }
    double std_s=std(s_11_real,1, m_SamleSize / 2);
    Complex* Signal_Channel_3 = new Complex[m_SamleSize / 2];
    memset(Signal_Channel_3, 0x00, sizeof(Complex) * m_SamleSize / 2);
    for (int m = 0; m < m_SamleSize / 2; m++)
    {
        Signal_Channel_3[m].IData = s_11[m].IData*(0.015/ std_s);
        Signal_Channel_3[m].QData = s_11[m].QData * (0.015 / std_s);
    }
    

    
    Complex* DataPLLBuff_BP = new Complex[m_SamleSize / 2];
    memset(DataPLLBuff_BP, 0x00, sizeof(Complex) * m_SamleSize / 2);
    Complex* DataFLLBuff_BP = new Complex[m_SamleSize / 2];
    memset(DataFLLBuff_BP, 0x00, sizeof(Complex) * m_SamleSize / 2);
    PM_QP.CarrierSync(Signal_Channel_3, DataPLLBuff_BP, DataFLLBuff_BP, m_SamleSize / 2, 2, &fPLLNCO_1, &fPLLPastFreqPart_1, &FLLBuffer);
    

    for (int i = 0; i < m_SamleSize / 2; i++) {
        DataPLLBuff_BP[i].IData = DataPLLBuff_BP[i].IData * 30;
        DataPLLBuff_BP[i].QData = DataPLLBuff_BP[i].QData * 30;
    }
    int nSymbolSyncSize = 0;//突发在当前片中的长度
    Complex* DataSymbolSyncBuff = new Complex[m_SamleSize / 2];
    memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize / 2);
    PM_QP.SymbolSync(DataPLLBuff_BP, DataSymbolSyncBuff, m_SamleSize / 2, &nSymbolSyncSize, &SySyncBuffer);

    //sprintf(file_path, "D:/data/symbolI.txt");
    //fp = fopen(file_path, "at");
    //for (int i = 0; i < nSymbolSyncSize; i++) {
    //    fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
    //}
    //fclose(fp);

    //sprintf(file_path, "D:/data/symbolQ.txt");
    //fp = fopen(file_path, "at");
    //for (int i = 0; i < nSymbolSyncSize; i++) {
    //    fprintf(fp, "%f\n", DataSymbolSyncBuff[i].QData);
    //}
    //fclose(fp);

	signal_PM_QP_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_PM_QP_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	signal_PM_QP_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

    char* c_symbolI = new char[nSymbolSyncSize];
	char* c_symbolQ = new char[nSymbolSyncSize];
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbolI, c_symbolQ);
    
    //sprintf(file_path, "D:/data/PM_QPSK_output_symbol_I.txt");
    //fp = fopen(file_path, "at");
    //for (int i = 0; i < nSymbolSyncSize; i++) {
    //	fprintf(fp, "%d\n", c_symbolI[i]);
    //}
    //fclose(fp);

    //sprintf(file_path, "D:/data/PM_QPSK_output_symbol_Q.txt");
    //fp = fopen(file_path, "at");
    //for (int i = 0; i < nSymbolSyncSize; i++) {
    //	fprintf(fp, "%d\n", c_symbolQ[i]);
    //}
    //fclose(fp);

    signal_PM_QP_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
    memcpy(signal_PM_QP_Out->burstData->nDemodulationByteI, c_symbolI, (sizeof(char)) * nSymbolSyncSize);

    signal_PM_QP_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSize];
    memcpy(signal_PM_QP_Out->burstData->nDemodulationByteQ, c_symbolQ, (sizeof(char)) * nSymbolSyncSize);

    DELETE_ARR(DataBuff);
    DELETE_ARR(DataPLLBuff_PM);
    DELETE_ARR(PLL_Freq_Part_2);
    DELETE_ARR(PLL_Phase_Part_2);
    DELETE_ARR(s_r1);
    DELETE_ARR(s_1);
    DELETE_ARR(s_11);
    DELETE_ARR(s_11_real);
    DELETE_ARR(Signal_Channel_3);
    DELETE_ARR(DataPLLBuff_BP);
    DELETE_ARR(DataFLLBuff_BP);
    DELETE_ARR(DataSymbolSyncBuff);
    DELETE_ARR(c_symbolI);
	DELETE_ARR(c_symbolQ);

}

void Demo_PM_QP::InitialDemodulation(const DemodulationInitParamater& info) {

    m_Rb = info.rb;
    fc2 = info.fc_demo[0];
 
    InitBlockCFilter();
    PM_QPInit();
    InitCostasPLL();
    InitSymbolSync();

}

void Demo_PM_QP::InitBlockCFilter() 
{
    // 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
    CFilter.nFilterTaps = 17;
    CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

    // 使用数组初始化的值
    float tempCoef[] = { -0.0300349916229418,-0.0283172616659510,1.03450095994087e-17,0.0606798464270379,
        0.150174958114709,0.254855354993559,0.353841408874743,0.424758924989266,0.450524874344127,
        0.424758924989266,0.353841408874743,0.254855354993559,0.150174958114709,0.0606798464270379,
        1.03450095994087e-17,-0.0283172616659510,-0.0300349916229418 };
    // 将初始化值复制到动态分配的数组中
    for (int i = 0; i < CFilter.nFilterTaps; ++i) {
        CFilter.fFilterCoef[i] = tempCoef[i];
    }
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    PM_QP.CFilter = CFilter;
}

void Demo_PM_QP::PM_QPInit() {
    m_sAlgDemInit.llSliceLength = 8192;
    m_SamleSize = m_sAlgDemInit.llSliceLength;
    PM_QP.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;
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
    SySyncBuffer.fSymbolSyncN = 0.9;                            //符号同步NCO寄存器，初值设为1
    SySyncBuffer.fSymbolSyncNTemp = 0.9;                        //符号同步NCO暂时的寄存器，初值设为1
    SySyncBuffer.fSymbolSyncU = 0.6;                            //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
    SySyncBuffer.nSymbolSyncKK = 0;                             //符号同步用来表示Ti时间序号,指示u,yI,yQ
    SySyncBuffer.fSymbolSyncPastW = 0;                          //符号同步W的缓存值
    SySyncBuffer.fSymbolSyncPastN = 0;                          //符号同步N的缓存值
    SySyncBuffer.fSymbolSyncTimeError1 = 0;                     //符号同步time_error缓存值
    SySyncBuffer.fSymbolSyncTimeError2 = 0;                     //符号同步time_error缓存值
}


void Demo_PM_QP::InitSymbolSync() 
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823/4;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05/16 ;
    PM_QP.m_sAlgDemInit.nSampPerSymb = 4;                       //每个码元中的采样点数

    PM_QP.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_PM_QP::InitCostasPLL() {
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;                                        //环路阻尼系数
    int Ko = 1;                                                 //压控振荡器增益
    int Kd = 1;                                                 //鉴相器增益
    float fBLcoef = 0.04;
    int K = Ko * Kd;
    BL = fBLcoef * m_Rb;
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (8 * (float)m_Rb);
    CostasPll.fPLLLoopFilterCoef1 = 2 * sigma * Wn * T_nco;
    CostasPll.fPLLLoopFilterCoef2 = Wn * T_nco * Wn * T_nco;
    PM_QP.CostasPll = CostasPll;
}

//void Demo_PM_QP::Judgment(Complex buffin[], int bufflen, char buffoutI[], char buffoutQ[])
//{
//
//    //int* output_i = new int[bufflen];
//    //memset(output_i, 0x00, sizeof(int) * bufflen);
//    //int* output_q = new int[bufflen];
//    //memset(output_q, 0x00, sizeof(int) * bufflen);
//    //for (int i = 0; i < bufflen; i++)
//    //{
//    //    if (buffin[i].IData > 0)
//    //    {
//    //        output_i[i] = 1;
//    //    }
//    //    else
//    //    {
//    //        output_i[i] = -1;
//    //    }
//    //    if (buffin[i].QData > 0)
//    //    {
//    //        output_q[i] = 1;
//    //    }
//    //    else
//    //    {
//    //        output_q[i] = -1;
//    //    }
//    //}
//    //for (int i = 0; i < bufflen; i++)
//    //{
//    //    if ((output_i[i] == -1) && (output_q[i] == -1))
//    //    {
//    //        buffout[i] = 0;
//    //    }
//    //    else if ((output_i[i] == 1) && (output_q[i] == -1))
//    //    {
//    //        buffout[i] = 1;
//    //    }
//    //    else if ((output_i[i] == 1) && (output_q[i] == 1))
//    //    {
//    //        buffout[i] = 2;
//    //    }
//    //    else if ((output_i[i] == -1) && (output_q[i] == 1))
//    //    {
//    //        buffout[i] = 3;
//    //    }
//    //}
//}
void Demo_PM_QP::Judgment(Complex buffin[], int bufflen, char buffoutI[], char buffoutQ[]) {

	for (int i = 0; i < bufflen; i++) {
		if (buffin[i].IData >= 0) {
			buffoutI[i] = 1;
		}
		else {
			buffoutI[i] = 0;
		}
		if (buffin[i].QData >= 0) {
			buffoutQ[i] = 1;
		}
		else {
			buffoutQ[i] = 0;
		}
	}
}
void Demo_PM_QP::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart, float* PLL_Freq_Part, float* PLL_Phase_Part)
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

void Demo_PM_QP::downphase(float BuffIn[], Complex BuffOut[], int nBuffSize, float* fDownConversionFc, double* fOrthDownConversionPhase)
{
    float phase0 = 2 * Pi * (*fDownConversionFc) / m_sAlgDemInit.nSigFs;
    for (int i = 0; i < nBuffSize; i++)
    {
        *fOrthDownConversionPhase += phase0;
        if (*fOrthDownConversionPhase >= 2 * Pi)
        {
            *fOrthDownConversionPhase -= 2 * Pi;
        }
        BuffOut[i].IData = BuffIn[i] * cos(*fOrthDownConversionPhase);
        BuffOut[i].QData = BuffIn[i] * -sin(*fOrthDownConversionPhase);
    }
}

double Demo_PM_QP::std(float* data, int start, int end)
{
    int i;
    double average, var = 0, sum = 0;
    for (i = start - 1; i < end; i++)
    {
        sum += (double)data[i];
    }

    average = sum / (end - start + 1);

    for (i = start - 1; i < end; i++)
    {
        var += pow(((double)data[i] - average), 2);
    }

    return sqrt(var / (end - start + 1));
}