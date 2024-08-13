#include "SignalSOQPSK.h"

Demo_SOQPSK::Demo_SOQPSK(int index):m_index(index)
{

}

Demo_SOQPSK::~Demo_SOQPSK()
{

}

bool Demo_SOQPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    m_demodulationInitParamater = demoParameter;
    subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb);
    return 1;
}

void Demo_SOQPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
    Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_SOQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_SOQPSK_Out) {

	//FILE*fp;
	//char file_path1[500];

	//FILE*fp2;
	//char file_path2[500];

    Complex* DataBuff = new Complex[m_SamleSize];
    memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);
    De_SOQPSK.BlockFilter(data_input_slice, m_SamleSize, DataBuff, cFilterBuff);
    

    De_SOQPSK.AGC(DataBuff, fAGCPastVc);



    Complex* DataPLLBuff = new Complex[m_SamleSize];
    memset(DataPLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
    Complex* DataFLLBuff = new Complex[m_SamleSize];
    memset(DataFLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
    De_SOQPSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, m_SamleSize, 4, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);
    


    int DelayOutLen = m_SamleSize - nDelaySign;
    if (DelayOutLen <= 0) {
        DelayOutLen = 1;
    }
    Complex* DataDelayOut = new Complex[DelayOutLen];
    memset(DataDelayOut, 0x00, sizeof(Complex) * DelayOutLen);
    DelayDelete(DataPLLBuff, m_SamleSize, DataDelayOut, cDelayBuff, nDelaySign);
    
    Complex* aI_bQ = new Complex[(DelayOutLen)/2];
    memset(aI_bQ, 0x00, sizeof(Complex) *(DelayOutLen)/2 );
    for (int m = 0; m < DelayOutLen; m=m+2)
    {
        aI_bQ[m/2].IData = DataDelayOut[m].IData;
        aI_bQ[m/2].QData = DataDelayOut[m].QData;
    }



    int nSymbolSyncSize = 0;//突发在当前片中的长度
    Complex* DataSymbolSyncBuff = new Complex[m_SamleSize];
    memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);
    De_SOQPSK.SymbolSync(aI_bQ, DataSymbolSyncBuff, DelayOutLen /2, &nSymbolSyncSize, &SySyncBuffer);
    
	signal_SOQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_SOQPSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))*nSymbolSyncSize);
	signal_SOQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;


    char* c_symbolI = new char[nSymbolSyncSize];
	char* c_symbolQ = new char[nSymbolSyncSize];
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbolI, c_symbolQ);

	signal_SOQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_SOQPSK_Out->burstData->nDemodulationByteI, c_symbolI, (sizeof(char))*nSymbolSyncSize);
	signal_SOQPSK_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSize];
	memcpy(signal_SOQPSK_Out->burstData->nDemodulationByteQ, c_symbolQ, (sizeof(char))*nSymbolSyncSize);
	//sprintf(file_path1, "D:/data/SOQPSK_output_symbol_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", c_symbolI[i]);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/SOQPSK_output_symbol_Q.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", c_symbolQ[i]);
	//}
	//fclose(fp);
    
    DELETE_ARR(DataBuff);
    DELETE_ARR(DataPLLBuff);
    DELETE_ARR(DataFLLBuff);
    DELETE_ARR(DataDelayOut);
    DELETE_ARR(aI_bQ);
    DELETE_ARR(DataSymbolSyncBuff);
    DELETE_ARR(c_symbolI);
	DELETE_ARR(c_symbolQ);

}

void Demo_SOQPSK::InitBlockCFilter() 
{
        // 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
        CFilter.nFilterTaps = 5;
        CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

        // 使用数组初始化的值
        float tempCoef[] = { 0.438773933912439,0.451333483027147,0.455577703619451,0.451333483027147,0.438773933912439 };

        // 将初始化值复制到动态分配的数组中
        for (int i = 0; i < CFilter.nFilterTaps; ++i) {
            CFilter.fFilterCoef[i] = tempCoef[i];
        }
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    De_SOQPSK.CFilter = CFilter;
}

void Demo_SOQPSK::InitialDemodulation(int Rb)
{
    m_Rb = Rb;
    InitBlockCFilter();
    SOQPSKInit();
    InitCostasPLL();
    InitSymbolSync();
}

void Demo_SOQPSK::SOQPSKInit() 
{
    m_sAlgDemInit.llSliceLength = 8192;
    m_SamleSize = m_sAlgDemInit.llSliceLength;
    De_SOQPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

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
    SySyncBuffer.fSymbolSyncN = 0.8;                        //符号同步NCO寄存器，初值设为1
    SySyncBuffer.fSymbolSyncNTemp = 0.8;                    //符号同步NCO暂时的寄存器，初值设为1
    SySyncBuffer.fSymbolSyncU = 0.6;                        //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
    SySyncBuffer.nSymbolSyncKK = 0;                         //符号同步用来表示Ti时间序号,指示u,yI,yQ
    SySyncBuffer.fSymbolSyncPastW = 0;                      //符号同步W的缓存值
    SySyncBuffer.fSymbolSyncPastN = 0;                      //符号同步N的缓存值
    SySyncBuffer.fSymbolSyncTimeError1 = 0;                 //符号同步time_error缓存值
    SySyncBuffer.fSymbolSyncTimeError2 = 0;                 //符号同步time_error缓存值

    cDelayBuff = new float[2];
    memset(cDelayBuff, 0x00, sizeof(float) * 2);

    nDelaySign = 2;  //第一片需要延迟两个采样点
}

void Demo_SOQPSK::InitSymbolSync() 
{
    symbolsyncsactorInit.fSymbolSyncFactor1 = 0.0002666398168165;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 3.555913481466967e-08;
    De_SOQPSK.m_sAlgDemInit.nSampPerSymb = 4;//每个码元中的采样点数

    De_SOQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_SOQPSK::InitCostasPLL() {
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;                                    //环路阻尼系数
    int Ko = 1;                                             //压控振荡器增益
    int Kd = 1;                                             //鉴相器增益
    float fBLcoef = 0.002;
    int K = Ko * Kd;
    BL = fBLcoef * m_Rb;
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (float)m_Rb;

    CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
    De_SOQPSK.CostasPll = CostasPll;
}

//void Demo_SOQPSK::Judgment(Complex buffin[], int bufflen, int buffout[]) {
//    int* y1 = new int[bufflen];
//    memset(y1, 0x00, sizeof(int) * bufflen);
//    int* y2 = new int[bufflen];
//    memset(y2, 0x00, sizeof(int) * bufflen);
//    for (int i=0;i< bufflen;i++)
//    {
//        if (buffin[i].IData>0)
//        {
//            y1[i] = 1;
//        }
//        else
//        {
//            y1[i] = 0;
//        }
//        if (buffin[i].QData > 0)
//        {
//            y2[i] = 1;
//        }
//        else
//        {
//            y2[i] = 0;
//        }
//    }
//    for (int i = 0; i < bufflen*2; i++)
//    {
//        if (i % 2 == 0) 
//        { 
//            buffout[i] = y2[i / 2];
//        }
//        else
//        {
//            buffout[i] = y1[(i-1) / 2];
//        }
//    }
//    DELETE_ARR(y1);
//    DELETE_ARR(y2);
//}
void Demo_SOQPSK::Judgment(Complex buffin[], int bufflen, char buffoutI[], char buffoutQ[]) {

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
void Demo_SOQPSK::DelayDelete(Complex buffin[], int bufflen, Complex buffout[], float DelayBuff[], int& nDelaySign) 
{
    if (nDelaySign == 2) {
        for (int i = 0; i < bufflen - 2; i++) {
            buffout[i].IData = buffin[i + 2].IData;
            buffout[i].QData = buffin[i].QData;
        }
        DelayBuff[0] = buffin[bufflen - 2].QData;
        DelayBuff[1] = buffin[bufflen - 1].QData;
        nDelaySign = 0;
    }
    else {
        buffout[0].QData = DelayBuff[0];
        buffout[1].QData = DelayBuff[1];
        for (int i = 0; i < bufflen - 2; i++) {
            buffout[i].IData = buffin[i].IData;
            buffout[i + 2].QData = buffin[i].QData;
        }
        buffout[bufflen - 2].IData = buffin[bufflen - 2].IData;
        buffout[bufflen - 1].IData = buffin[bufflen - 1].IData;
        DelayBuff[0] = buffin[bufflen - 2].QData;
        DelayBuff[1] = buffin[bufflen - 1].QData;
    }

}