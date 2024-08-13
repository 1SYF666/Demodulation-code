#include "Signal8PSK.h"

Demo_8PSK::Demo_8PSK(int index):m_index(index)
{

}

Demo_8PSK::~Demo_8PSK()
{

}

bool Demo_8PSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    InitialDemodulation(demoParameter.rb);
    return 1;
}

#include  <QFile>
#include <QDataStream>
void Demo_8PSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

}

void Demo_8PSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_8PSK_Out)
{
    // add on 20240725
    if (nSilceCount == 1)
    {
        // 第二次进入该函数，载波系数保持不变
        De_8PSK.CostasPll = CostasPll3;   
        De_8PSK.CostasPll2 = CostasPll3;   
        De_8PSK.CostasPll3 = CostasPll3;   

    }
    nSilceCount++;

    int samplelen = PSK8_Samplesize;                        // add on 20240728

    Complex* DataBuff = new Complex[samplelen];
    memset(DataBuff, 0x00, sizeof(Complex) * samplelen);
    for (int m = 0; m < samplelen; m++)
    {
        DataBuff[m].IData = data_input_slice[m].IData;
        DataBuff[m].QData = data_input_slice[m].QData;
    }
    De_8PSK.AGC(DataBuff, fAGCPastVc, samplelen);

    //Complex* blockdata = new Complex[samplelen];
    //memset(blockdata, 0x00, sizeof(Complex) * samplelen);
    //De_8PSK.BlockFilter(DataBuff, samplelen, blockdata, cFilterBuff);

    De_8PSK.AGC(DataBuff, fAGCPastVc, samplelen);

    Complex* DataPLLBuff = new Complex[samplelen];
    memset(DataPLLBuff, 0x00, sizeof(Complex) * samplelen);
    Complex* DataFLLBuff = new Complex[samplelen];
    memset(DataFLLBuff, 0x00, sizeof(Complex) * samplelen);
    
    De_8PSK.CarrierSync(DataBuff, DataPLLBuff, DataFLLBuff, samplelen, 11, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

    Complex* aI_bQ = new Complex[samplelen];
    memset(aI_bQ, 0x00, sizeof(Complex) * samplelen);
    De_8PSK.BlockFilter(DataPLLBuff, samplelen, aI_bQ, cFilterBuff);

    int nSymbolSyncSize = 0;
    Complex* DataSymbolSyncBuff = new Complex[samplelen];
    memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * samplelen);
    De_8PSK.SymbolSync(aI_bQ, DataSymbolSyncBuff, samplelen, &nSymbolSyncSize, &SySyncBuffer);

    char* c_symbol = new char[nSymbolSyncSize];
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbol);

	signal_8PSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
	memcpy(signal_8PSK_Out->burstData->softDistinguishData, DataSymbolSyncBuff, (sizeof(Complex))* nSymbolSyncSize);
	signal_8PSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;

	signal_8PSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_8PSK_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char))*nSymbolSyncSize);


    DELETE_ARR(DataBuff);
    DELETE_ARR(DataPLLBuff);
    DELETE_ARR(DataFLLBuff);
    DELETE_ARR(aI_bQ);
    DELETE_ARR(DataSymbolSyncBuff);
    DELETE_ARR(c_symbol);

}

void Demo_8PSK::InitialDemodulation(int Rb) {

    m_Rb = Rb;
    PSK_8_Init();
    InitBlockCFilter();
    InitCostasPLL();
    InitSymbolSync();
}

void Demo_8PSK::InitBlockCFilter() 
{
    // 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
    CFilter.nFilterTaps = 9;
    CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

    // 使用数组初始化的值
    float tempCoef[] = { -0.0424636059371430,1.46258209943935e-17,0.212318029685715,0.500262571713977,0.636954089057146,
                          0.500262571713977,0.212318029685715,1.46258209943935e-17,-0.0424636059371430 };
    // 将初始化值复制到动态分配的数组中
    for (int i = 0; i < CFilter.nFilterTaps; ++i) {
        CFilter.fFilterCoef[i] = tempCoef[i];
    }
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    De_8PSK.CFilter = CFilter;
}

void Demo_8PSK::PSK_8_Init() {
    m_sAlgDemInit.llSliceLength = 8192;
    m_SamleSize = m_sAlgDemInit.llSliceLength;
    De_8PSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

    nDownSampClock = 1;                                         //降采样时的时钟控制  
    fDownConversionPhase = 0;                                   //下变频中的相位值
    fPLLNCO = 0;                                                //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                       //锁相环中的频率跟踪曲线
    nPLLBuffSize = 0;

    // add on 20240621
    flag = 0;
    slice = 10;
    PSK8_Samplesize = slice * m_SamleSize;
    Databuff0 = new Complex[PSK8_Samplesize];

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

void Demo_8PSK::InitSymbolSync() {
    //symbolsyncsactorInit.fSymbolSyncFactor1 = 0.002666398168165;
    //symbolsyncsactorInit.fSymbolSyncFactor2 = 3.555913481466967e-06; //改7-8-10:29
	symbolsyncsactorInit.fSymbolSyncFactor1 = 2.666398168164524e-4;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 3.555913481466967e-09; //改7-8-10:29

    De_8PSK.m_sAlgDemInit.nSampPerSymb = 4;                       //每个码元中的采样点数

    De_8PSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_8PSK::InitCostasPLL() {
    float BL;
    float Wn;
    float T_nco;
    float sigma = 0.707;                                        //环路阻尼系数
	int Ko = 1;                                                 //压控振荡器增益
    int Kd = 1;                                                 //鉴相器增益
    int K = Ko * Kd;

    // add on 20240724
    // 采用分段形式
    float temp1 = 0;
    float temp2 = 0;

	//float fBLcoef = 0.02;
    //float fBLcoef = 0.03;                                       // add on 20240724    
    //float fBLcoef = 0.0035;                                       // add on 20240724,效果比前两个好。
    //float fBLcoef = 0.0030;                                       // add on 20240724    
    float fBLcoef = 0.002;                                       // add on 20240724    
    BL = fBLcoef * m_Rb;                                        // add on 20240724
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (4 * (float)m_Rb);
    CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

    // 第二段
	fBLcoef = 0.0002;
    BL = fBLcoef * m_Rb;                                        // add on 20240724
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (4 * (float)m_Rb);
    CostasPll2.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll2.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

    // 第三段
    fBLcoef = 0.00001;
    BL = fBLcoef * m_Rb;                                        // add on 20240724
    Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
    T_nco = 1 / (4 * (float)m_Rb);
    CostasPll3.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
    CostasPll3.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	De_8PSK.CostasPll = CostasPll;   // 第一段
	De_8PSK.CostasPll2 = CostasPll2;   // 第二段
	De_8PSK.CostasPll3 = CostasPll3;   // 第三段
}

void Demo_8PSK::Judgment(Complex buffin[], int bufflen, char buffout[]) {
    float theta = 0;
    int* y = new int[bufflen];
    memset(y, 0x00, sizeof(int) * bufflen);
    for (int i=0;i< bufflen;i++)
    {
        theta = atan2(buffin[i].QData, buffin[i].IData);
        if (theta <= 0)
        {
            theta = theta+ 2 * Pi;
        }
        if (theta >= 0 && theta < Pi / 4)
            y[i] = 0;
        else if(theta >= Pi / 4 && theta < 4 * Pi / 8)
            y[i] = 1;
        else if(theta >= 4 * Pi / 8 && theta < 6 * Pi / 8)
            y[i] = 2;
        else if(theta >= 6 * Pi / 8 && theta < 8 * Pi / 8)
            y[i] = 3;
        else if(theta >= 8 * Pi / 8 && theta < 10 * Pi / 8)
            y[i] = 4;
        else if(theta >= 10 * Pi / 8 && theta < 12 * Pi / 8)
            y[i] = 5;
        else if(theta >= 12 * Pi / 8 && theta < 14 * Pi / 8)
            y[i] = 6;
        else if(theta >= 14 * Pi / 8)
            y[i] = 7;
    }

    for (int i = 0; i < bufflen; i++)
    {
        if (i==0)
        {
            buffout[i] = ((y[i]- difftemp)+8) % 8;
        }
        else 
        {
            if (y[i] - y[i - 1] < 0)
            {
                buffout[i] = y[i] - y[i - 1]+ 8;
            }
            else
            {
                buffout[i] = y[i] - y[i - 1];
            }
        }
    }
	difftemp = y[bufflen-1];
    DELETE_ARR(y);
}

