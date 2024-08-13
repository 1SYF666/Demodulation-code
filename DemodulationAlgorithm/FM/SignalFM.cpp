#include "SignalFM.h"

Demo_FM::Demo_FM(int index):m_index(index)
{

}

Demo_FM::~Demo_FM() 
{

}

bool Demo_FM::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
    m_demodulationInitParamater = demoParameter;
    subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb);
    return 1;
}

void Demo_FM::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
    Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}



void Demo_FM::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_BPSK_Out)
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
	FM.AGC(DataBuff, fAGCPastVc);

	//sprintf(file_path1, "D:/data/input_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path2, "D:/data/input_Q.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < m_SamleSize; i++) {
	//	fprintf(fp2, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp2);

	float* signal_base = new float[m_SamleSize];
	memset(signal_base, 0x00, sizeof(float) * m_SamleSize);
	for (int i = 0; i < m_SamleSize; i++)
    {
		float x = DataBuff[i].IData * DataBuff[i].IData + DataBuff[i].QData * DataBuff[i].QData;
		if (x == 0) {
			signal_base[i] = (tempI * tempI - tempQ * tempQ);
		}
		else{
			signal_base[i] = (tempI * DataBuff[i].QData - DataBuff[i].IData * tempQ) / (DataBuff[i].IData * DataBuff[i].IData + DataBuff[i].QData * DataBuff[i].QData);
		}
		tempI = DataBuff[i].IData;
		tempQ = DataBuff[i].QData;
    }
    float signal_base_mean=mean_function(signal_base,0, m_SamleSize);
    for (int i = 0; i < m_SamleSize; i++)
    {
        signal_base[i] = signal_base[i]- signal_base_mean;
    }
   
    Complex* signal_base_I = new Complex[m_SamleSize];
    memset(signal_base_I, 0x00, sizeof(Complex) * m_SamleSize);
    for (int i = 1; i < m_SamleSize; i++)
    {
        signal_base_I[i].IData = signal_base[i];
    }
    
    Complex* yi = new Complex[m_SamleSize];
    memset(yi, 0x00, sizeof(Complex) * m_SamleSize);
    
    FM.BlockFilter(signal_base_I, m_SamleSize, yi, cFilterBuff);

	//sprintf(file_path1, "D:/data/FM_output.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < m_SamleSize; i++)
	//{
	//	fprintf(fp, "%f\n", yi[i].IData);
	//}
	//fclose(fp);

    DELETE_ARR(DataBuff);
    DELETE_ARR(signal_base);
    DELETE_ARR(signal_base_I);
    DELETE_ARR(yi);


}

void Demo_FM::InitBlockCFilter()
{
        // 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
        CFilter.nFilterTaps = 83;
        CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

        // 使用数组初始化的值
        float tempCoef[] = { 0.000614825385536176,0.000614308087911696,0.000609378904822616,0.000593117053905705,0.000555021333278777,0.000481751986983585,
            0.000358266387210140,0.000169290789833929,- 9.89582109658489e-05,- 0.000456912963951432,- 0.000909782909310685,- 0.00145580962305395,
            - 0.00208478791124812,- 0.00277698402543723,- 0.00350256492404075,- 0.00422162702388059,- 0.00488488035651076,- 0.00543500623967611,
            - 0.00580866565872870,- 0.00593909395397714,- 0.00575917767467252,- 0.00520487409927856,- 0.00421880525277775,- 0.00275383825270562,
            - 0.000776453994363081,0.00173029252241509,0.00476440422438322,0.00830332662966728,0.0123032117182940,0.0166991455574501,0.0214063760757837,
            0.0263225256889126,0.0313307211017101,0.0363035224990164,0.0411074893537395,0.0456081828499538,0.0496753776532034,0.0531882401355315,
            0.0560402272102006,0.0581434699724972,0.0594324289661824,0.0598666419723104,0.0594324289661824,0.0581434699724972,0.0560402272102006,
            0.0531882401355315,0.0496753776532034,0.0456081828499538,0.0411074893537395,0.0363035224990164,0.0313307211017101,0.0263225256889126,
            0.0214063760757837,0.0166991455574501,0.0123032117182940,0.00830332662966728,0.00476440422438322,0.00173029252241509,- 0.000776453994363081,
            - 0.00275383825270562,- 0.00421880525277775,- 0.00520487409927856,- 0.00575917767467252,- 0.00593909395397714,- 0.00580866565872870,
            - 0.00543500623967611,- 0.00488488035651076,- 0.00422162702388059,- 0.00350256492404075,- 0.00277698402543723,- 0.00208478791124812,
            - 0.00145580962305395,- 0.000909782909310685,- 0.000456912963951432,- 9.89582109658489e-05,0.000169290789833929,0.000358266387210140,
            0.000481751986983585,0.000555021333278777,0.000593117053905705,0.000609378904822616,0.000614308087911696,0.000614825385536176 };
        // 将初始化值复制到动态分配的数组中
        for (int i = 0; i < CFilter.nFilterTaps; ++i) {
            CFilter.fFilterCoef[i] = tempCoef[i];
        }
    cFilterBuff = new Complex[CFilter.nFilterTaps];
    memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    FM.CFilter = CFilter;
}

void Demo_FM::InitialDemodulation(int Rb) //全部初始化
{

    m_Rb = Rb; 

    FM_Init();

    InitBlockCFilter();

}

void Demo_FM::FM_Init()
{
    m_sAlgDemInit.llSliceLength = 8192;
    m_SamleSize = m_sAlgDemInit.llSliceLength;
    FM.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

    nDownSampClock = 1;                                             //降采样时的时钟控制  
    fDownConversionPhase = 0;                                       //下变频中的相位值
    fPLLNCO = 0;                                                    //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                           //锁相环中的频率跟踪曲线
    nPLLBuffSize = 0;
    //******************符号同步初始化***********
    SySyncBuffer.CSymbolSyncBuff = new Complex[4];
    SySyncBuffer.CSymbolSyncYBuff = new Complex[3];
    memset(SySyncBuffer.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
    memset(SySyncBuffer.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
    SySyncBuffer.fSymbolSyncW = 0.5;                                //符号同步环路滤波器输出寄存器，初值设为0.5
    SySyncBuffer.fSymbolSyncN = 0.9;                                //符号同步NCO寄存器，初值设为1
    SySyncBuffer.fSymbolSyncNTemp = 0.9;                            //符号同步NCO暂时的寄存器，初值设为1
    SySyncBuffer.fSymbolSyncU = 0.6;                                //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
    SySyncBuffer.nSymbolSyncKK = 0;                                 //符号同步用来表示Ti时间序号,指示u,yI,yQ
    SySyncBuffer.fSymbolSyncPastW = 0;                              //符号同步W的缓存值
    SySyncBuffer.fSymbolSyncPastN = 0;                              //符号同步N的缓存值
    SySyncBuffer.fSymbolSyncTimeError1 = 0;                         //符号同步time_error缓存值
    SySyncBuffer.fSymbolSyncTimeError2 = 0;                         //符号同步time_error缓存值
}

float Demo_FM::mean_function(float* data, int start, int end)
{
	vector<float>sum1(8192, 0);
    float sum = 0;
    int i = 0;
    for (i = start; i < end; i++)
    {
        sum += data[i];
		sum1[i] = sum;
    }
    return sum / (end - start);
}