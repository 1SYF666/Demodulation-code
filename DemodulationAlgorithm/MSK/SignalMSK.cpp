#include "SignalMSK.h"

Demo_MSK::Demo_MSK(int index):m_index(index)
{

}

Demo_MSK::~Demo_MSK() 
{
	DELETE_ARR(Databuff0);
}

bool Demo_MSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
    InitialDemodulation(demoParameter.rb);
    return 1;
}

void Demo_MSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_MSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_MSK_Out) 
{
	int datalength = MSK_Samplesize;


	//FILE*fp;
	//char file_path1[500];

	/*sprintf(file_path1, "D:/data/inputI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < datalength; i++)
	{
		fprintf(fp, "%f\n", data_input_slice[i].IData);
	}
	fclose(fp);

	sprintf(file_path1, "D:/data/inputQ.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < datalength; i++)
	{
		fprintf(fp, "%f\n", data_input_slice[i].QData);
	}
	fclose(fp);*/


    //De_MSK.AGC(data_input_slice, fAGCPastVc);
	AGCMSK(data_input_slice, fAGCPastVc, datalength);
	
    Complex* DataPLLBuff = new Complex[datalength];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * datalength);
	Complex* DataFLLBuff = new Complex[datalength];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * datalength);

    De_MSK.CarrierSync(data_input_slice, DataPLLBuff, DataFLLBuff, datalength, 9, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);
    

	/*sprintf(file_path1, "D:/data/MSK_output_symbolI.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < datalength; i++)
	{
		fprintf(fp, "%f\n", DataFLLBuff[i].IData);
	}
	fclose(fp);

	sprintf(file_path1, "D:/data/MSK_output_symbolQ.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < datalength; i++)
	{
		fprintf(fp, "%f\n", DataFLLBuff[i].QData);
	}
	fclose(fp);*/


    Complex* DataBuff = new Complex[datalength];
    memset(DataBuff, 0x00, sizeof(Complex) * datalength);

    De_MSK.BlockFilter(DataFLLBuff, datalength, DataBuff, cFilterBuff);


    Complex* DataNorm = new Complex[datalength];
    memset(DataNorm, 0x00, sizeof(Complex) * m_SamleSize);
    De_MSK.normolize(DataBuff, datalength, DataNorm);


	int nSymbolSyncSize = 0;
	Complex* DataSymbolSyncBuff = new Complex[datalength];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * datalength);

    De_MSK.SymbolSync(DataNorm, DataSymbolSyncBuff, datalength, &nSymbolSyncSize, &SySyncBuffer);


    signal_MSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSize];
    memcpy(signal_MSK_Out->burstData->softDistinguishData,DataSymbolSyncBuff,(sizeof (Complex))*nSymbolSyncSize);
    signal_MSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
    
    char* DataResult = new char[nSymbolSyncSize];
    memset(DataResult, 0x00, sizeof(char) * nSymbolSyncSize);
    Judgment(DataSymbolSyncBuff, nSymbolSyncSize, DataResult);


	//sprintf(file_path1, "D:/data/MSK_output_symbol.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	fprintf(fp, "%d\n", DataResult[i]);
	//}
	//fclose(fp);
	
	signal_MSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_MSK_Out->burstData->nDemodulationByteI, DataResult, (sizeof(char)) * nSymbolSyncSize);

    DELETE_ARR(DataBuff);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataResult);
    DELETE_ARR(DataNorm);
	DELETE_ARR(DataSymbolSyncBuff);
}

void Demo_MSK::InitBlockCFilter() 
{
	CFilter.nFilterTaps =33;
	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

	// 使用数组初始化的值
	float tempCoef[] = { 0.000932947996071962,0.00188009601734038,0.00153290572833484,-0.00118865119117817,-0.00540241983579544,-0.00664080691735891,4.24483894033689e-18,0.0128483025714876,0.0203788111596372,0.00882501430630615,-0.0222674279634046,-0.0505094850439060,-0.0403625545478313,0.0301637228733058,0.145668434287132,0.254548590542909,0.299185040033900,0.254548590542909,0.145668434287132,0.0301637228733058,-0.0403625545478313,-0.0505094850439060,-0.0222674279634046,0.00882501430630615,0.0203788111596372,0.0128483025714876,4.24483894033689e-18,-0.00664080691735891,-0.00540241983579544,-0.00118865119117817,0.00153290572833484,0.00188009601734038,0.000932947996071962 };

	// 将初始化值复制到动态分配的数组中
	for (int i = 0; i < CFilter.nFilterTaps; ++i) {
		CFilter.fFilterCoef[i] = tempCoef[i];
	}
	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);
    CFilter.bCoefEvenSym = false;
    De_MSK.CFilter = CFilter;
}

void Demo_MSK::InitialDemodulation(int Rb) 
{
    m_Rb = Rb;
	InitBlockCFilter();
    MSKInit();
    InitSymbolSync();
    InitFLL();
}

void Demo_MSK::MSKInit() 
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
    De_MSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;


    //FLL初始化
    FLLBuffer.CFLLbuff = new Complex[4/ 2];
    memset(FLLBuffer.CFLLbuff, 0x00, sizeof(Complex) * (4 / 2));	               //nSampPerSymb/2长度  //FLL中间缓存区
    FLLBuffer.fFLLIaccum = 0;                                   //FLL中I_accum的缓存值
    FLLBuffer.fFLLQaccum = 0;                                   //FLL中Q_accum的缓存值
    FLLBuffer.fFLLfreqoutdlf = 0;                               //FLL中freqoutdlf的上一个缓存值
    FLLBuffer.fFLLPastfreqoutdlf = 0;                           //FLL中freqoutdlf的前两个缓存值
    FLLBuffer.fFLLfreqerror = 0;                                //FLL中freqerror的上一个缓存值
    FLLBuffer.fFLLNCOPhase = 0;                                 //FLL中NCOPhase的上一个缓存值
    //
    nDownSampClock = 1;                                         //降采样时的时钟控制  
    fDownConversionPhase = 0;                                   //下变频中的相位值
    fPLLNCO = 0;                                                //锁相环中的本地NCO
    fPLLPastFreqPart = 0;                                       //锁相环中的频率跟踪曲线
    nPLLBuffSize = 0;
    
	
	// add on 20240607
	flag = 0;
	slice = 10;
	MSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[MSK_Samplesize];


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

void Demo_MSK::InitSymbolSync()
{
    symbolsyncsactorInit.fSymbolSyncFactor1 = 0.006665995420411/20;
    symbolsyncsactorInit.fSymbolSyncFactor2 = 2.222445925916854e-05/400;
    De_MSK.m_sAlgDemInit.nSampPerSymb = 4;                      //每个码元中的采样点数

    De_MSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_MSK::Judgment(Complex buffin[], int bufflen, char buffout[]) {

    for (int i = 0; i < bufflen; i++) 
	{
        if (buffin[i].IData <0 ) buffout[i] = 0;
        else if(buffin[i].IData >0) buffout[i] = 1;
    }

}

void Demo_MSK::InitFLL()
{
    float fSymbolTime = 1.0 / (float)m_Rb;
    float fDownSampleFs = m_Rb * (float)4;
    Fll.fFLLDownSampleTs = 1.0 / ((float)(m_Rb * 4));
    float fCarrierSyncFactor = 0.000002;
    float b_fll = fCarrierSyncFactor * m_Rb;
    float wn = 1.89 * b_fll;
    float gain_nco = 2 * Pi * fDownSampleFs;
    int gain_disctim = 1;
    int K = gain_nco * gain_disctim;
    Fll.fFLLPLLLoopFilterCoef1 = (sqrt(2) * wn * fSymbolTime + (wn * wn) * (fSymbolTime * fSymbolTime)) / K;
    Fll.fFLLPLLLoopFilterCoef2 = sqrt(2) * wn * fSymbolTime / K;
    De_MSK.Fll = Fll;
}


void Demo_MSK::AGCMSK(Complex* BuffIn, float target_power, int len)
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


float Demo_MSK::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}