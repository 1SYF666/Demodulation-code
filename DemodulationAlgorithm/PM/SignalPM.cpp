#include "SignalPM.h"

Demo_PM::Demo_PM(int index):m_index(index)
{

}

Demo_PM::~Demo_PM() {

}

bool Demo_PM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb);
    return 1;
}

void Demo_PM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
    Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}


void Demo_PM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_AM_Out)
{

    Demo.AGC(data_input_slice, fAGCPastVc);

	//FILE*fp;
	//char file_path1[500];


	Complex* DataPLLBuff_PM = new Complex[m_SamleSize];
	memset(DataPLLBuff_PM, 0x00, sizeof(Complex) * m_SamleSize);
	float* PLL_Freq_Part_2 = new float[m_SamleSize];
	memset(PLL_Freq_Part_2, 0x00, sizeof(float) * m_SamleSize);
	float* PLL_Phase_Part_2 = new float[m_SamleSize];
	memset(PLL_Phase_Part_2, 0x00, sizeof(float) * m_SamleSize);

	PLLCostas(data_input_slice, DataPLLBuff_PM, m_SamleSize, &fPLLNCO, &fPLLPastFreqPart, PLL_Freq_Part_2, PLL_Phase_Part_2);

	Complex *PLL_Idata = new Complex[m_SamleSize];
	for (int i = 0; i < m_SamleSize; i++) {
		PLL_Idata[i].IData = PLL_Freq_Part_2[i];
		PLL_Idata[i].QData = PLL_Freq_Part_2[i];
	}

	//sprintf(file_path1, "D:/data/PLLi.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", DataPLLBuff[i].IData);
	//}
	//fclose(fp);


    Complex* DataBuff = new Complex[m_SamleSize];
    memset(DataBuff, 0x00, sizeof(Complex) * m_SamleSize);


    Demo.BlockFilter(PLL_Idata, m_SamleSize, DataBuff, cFilterBuff);

	//sprintf(file_path1, "D:/data/PM_output.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)	
	//{
	//	fprintf(fp, "%f\n", DataBuff[i].IData);
	//}
	//fclose(fp);
	
    DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff_PM);
	DELETE_ARR(PLL_Idata);
	DELETE_ARR(PLL_Freq_Part_2);
	DELETE_ARR(PLL_Phase_Part_2);
    

}

void Demo_PM::InitBlockCFilter() {

	CFilter.nFilterTaps = 132;
	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

	// 使用数组初始化的值
	float tempCoef[] = { 0.000105268961886467,0.000193705548487144,0.000278719387505542,0.000358016028082931,0.000428391818861667,0.000485569168137571,0.000524215020346765,0.000538171455812717,0.000520904737780625,0.000466154048671028,0.000368735858933526,0.000225436782635963,3.59092570513747e-05,-0.000196527469251344,-0.000464282117965790,-0.000755177286522161,-0.00105250657296325,-0.00133551878375984,-0.00158035325295247,-0.00176141297281933,-0.00185312213394798,-0.00183197493677846,-0.00167874649980388,-0.00138070761923682,-0.000933665969599617,-0.000343649432318946,0.000371945841405220,0.00118389869495142,0.00205151835142912,0.00292388508098020,0.00374199192386439,0.00444173960872347,0.00495767112252821,0.00522726927164842,0.00519558473760214,0.00481991805499349,0.00407425052936252,0.00295310939437777,0.00147456336608248,-0.000317923256744530,-0.00235499754653989,-0.00454243535513507,-0.00676373407419520,-0.00888425163769842,-0.0107567648030805,-0.0122282242275760,-0.0131473969660605,-0.0133730147527884,-0.0127819940756048,-0.0112772658017892,-0.00879475076000480,-0.00530904438023993,-0.000837427641650252,0.00455809919520673,0.0107709633269656,0.0176524340451940,0.0250162679264463,0.0326454208422367,0.0403005426705828,0.0477298636052160,0.0546799940892735,0.0609070992397736,0.0661878772749839,0.0703297722321766,0.0731798845494423,0.0746321072773874,0.0746321072773874,0.0731798845494423,0.0703297722321766,0.0661878772749839,0.0609070992397736,0.0546799940892735,0.0477298636052160,0.0403005426705828,0.0326454208422367,0.0250162679264463,0.0176524340451940,0.0107709633269656,0.00455809919520673,-0.000837427641650252,-0.00530904438023993,-0.00879475076000480,-0.0112772658017892,-0.0127819940756048,-0.0133730147527884,-0.0131473969660605,-0.0122282242275760,-0.0107567648030805,-0.00888425163769842,-0.00676373407419520,-0.00454243535513507,-0.00235499754653989,-0.000317923256744530,0.00147456336608248,0.00295310939437777,0.00407425052936252,0.00481991805499349,0.00519558473760214,0.00522726927164842,0.00495767112252821,0.00444173960872347,0.00374199192386439,0.00292388508098020,0.00205151835142912,0.00118389869495142,0.000371945841405220,-0.000343649432318946,-0.000933665969599617,-0.00138070761923682,-0.00167874649980388,-0.00183197493677846,-0.00185312213394798,-0.00176141297281933,-0.00158035325295247,-0.00133551878375984,-0.00105250657296325,-0.000755177286522161,-0.000464282117965790,-0.000196527469251344,3.59092570513747e-05,0.000225436782635963,0.000368735858933526,0.000466154048671028,0.000520904737780625,0.000538171455812717,0.000524215020346765,0.000485569168137571,0.000428391818861667,0.000358016028082931,0.000278719387505542,0.000193705548487144,0.000105268961886467 };

	// 将初始化值复制到动态分配的数组中
	for (int i = 0; i < CFilter.nFilterTaps; ++i) {
		CFilter.fFilterCoef[i] = tempCoef[i];
	}

	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
	CFilter.bCoefEvenSym = false;
	Demo.CFilter = CFilter;
}

void Demo_PM::InitialDemodulation(int Rb) {

    Rb = 1.25e6;
    m_Rb = Rb;
	InitBlockCFilter();
    AMInit();
    InitCostasPLL();
}

void Demo_PM::AMInit() {
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
    Demo.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

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


    cDelayBuff = new float[2];
    memset(cDelayBuff, 0x00, sizeof(float) * 2);

    nDelaySign = 2;                                         //第一片需要延迟两个采样点
}



void Demo_PM::InitCostasPLL() {

    CostasPll.fPLLLoopFilterCoef1 = 0.011665491985720;
    CostasPll.fPLLLoopFilterCoef2 = 6.806240648120366e-05;
    Demo.CostasPll = CostasPll;
}

float Demo_PM::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_PM::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart, float* PLL_Freq_Part, float* PLL_Phase_Part)
{
	float discriminator_out = 0;
	float pll_phase_part;
	float freq_control;
	float past_pll_freq_part;
	float past_nco_phase;
	float PLLCoef1 = 0.010580945111764;
	float PLLCoef2 = 5.599511025237728e-05;
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