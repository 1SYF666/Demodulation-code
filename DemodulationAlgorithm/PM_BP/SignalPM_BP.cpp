#include "SignalPM_BP.h"
//#include<QFile>
//#incldue<QTextStream>

Demo_PM_BP::Demo_PM_BP(int index) :m_index(index)
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
	memcpy(Databuff0 + flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;


	Demodulation(Databuff0, baseTime, demodulationResult);
}

void Demo_PM_BP::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_PM_BP_Out)
{

	FILE* fp;
	char file_path[500];

	//sprintf(file_path, "D:/data/data_input_sliceI.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < PMBP_Samplesize; i++) {
	//	fprintf(fp, "%f\n", data_input_slice[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path, "D:/data/data_input_sliceQ.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < PMBP_Samplesize; i++) {
	//	fprintf(fp, "%f\n", data_input_slice[i].QData);
	//}
	//fclose(fp);




	// add on 20240730
	if (nSilceCount == 1)
	{
		// 第二次进入该函数，载波系数保持不变
		CostasPll = CostasPll2;
		CostasPll_bpsk = CostasPll3_bpsk;
		CostasPll2_bpsk = CostasPll3_bpsk;

	}
	nSilceCount++;

	PM_BP.AGC(data_input_slice, fAGCPastVc, PMBP_Samplesize);

	Complex* DataPLLBuff_PM = new Complex[PMBP_Samplesize];
	memset(DataPLLBuff_PM, 0x00, sizeof(Complex) * PMBP_Samplesize);
	float* PLL_Freq_Part_2 = new float[PMBP_Samplesize];
	memset(PLL_Freq_Part_2, 0x00, sizeof(float) * PMBP_Samplesize);
	float* PLL_Phase_Part_2 = new float[PMBP_Samplesize];
	memset(PLL_Phase_Part_2, 0x00, sizeof(float) * PMBP_Samplesize);

	PLLCostas(data_input_slice, DataPLLBuff_PM, PMBP_Samplesize, &fPLLNCO, &fPLLPastFreqPart, PLL_Freq_Part_2, PLL_Phase_Part_2);

	//sprintf(file_path, "D:/data/PLL_Phase_Part_2.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < PMBP_Samplesize; i++) {
	//	fprintf(fp, "%f\n", PLL_Phase_Part_2[i]);
	//}
	//fclose(fp);

	for (int i = 0; i < PMBP_Samplesize; i++)
	{
		PLL_Phase_Part_2[i] = PLL_Phase_Part_2[i] * 100;
	}

	Complex* de_pm_bpsk = new Complex[PMBP_Samplesize];
	memset(de_pm_bpsk, 0x00, sizeof(Complex) * PMBP_Samplesize);
	downphase(PLL_Phase_Part_2, de_pm_bpsk, PMBP_Samplesize, &fcsub, &fOrthDownConversionPhase);

	Complex* de_bpsk = new Complex[PMBP_Samplesize];
	memset(de_bpsk, 0x00, sizeof(Complex) * PMBP_Samplesize);
	PM_BP.BlockFilter(de_pm_bpsk, PMBP_Samplesize, de_bpsk, cFilterBuff);

	// PM_BPSK码速估计

	if (flag1 == 0)
	{
		float est_rb = 5000;               
		int fftlen = 16384 * 64;
		est_pmbp_rb(de_bpsk, PMBP_Samplesize, fftlen, fs, &est_rb);

		// 假设码速是以K为单位的
		int rb_temp = (est_rb + 500) / 1000;
		m_Rb = rb_temp * 1000;

		qDebug() << "m_Rb" << m_Rb;

		//m_Rb = 16e3;
		InitCostasPLL_bpsk();

		flag1 = 1;
	}

	// 降采样操作
	Complex* bpsk_temp = new Complex[PMBP_Samplesize];
	memset(bpsk_temp, 0x00, sizeof(Complex) * PMBP_Samplesize);

	int rb_temp = m_Rb;
	double step = ((double)(4 * rb_temp)) / ((double)fs);
	int nSymbolSyncSize_temp = 0;	       //降采样之后长度
	Downfs(de_bpsk, bpsk_temp, step, PMBP_Samplesize, &nSymbolSyncSize_temp);

	// BPSK解调
	PM_BP.AGC(bpsk_temp, fAGCPastVc, nSymbolSyncSize_temp);

	Complex* DataPLLBuff_BP = new Complex[nSymbolSyncSize_temp];
	memset(DataPLLBuff_BP, 0x00, sizeof(Complex) * nSymbolSyncSize_temp);
	// PM_BP.CarrierSync(bpsk_temp, DataPLLBuff_BP, DataFLLBuff_BP, nSymbolSyncSize_temp, 1, &fPLLNCO_1, &fPLLPastFreqPart_1, &FLLBuffer);
	PLLCostas_bpsk(bpsk_temp, DataPLLBuff_BP, nSymbolSyncSize_temp, &fPLLNCO_1, &fPLLPastFreqPart_1);

	//sprintf(file_path, "D:/data/DataPLLBuff_BPI.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < nSymbolSyncSize_temp; i++) {
	//	fprintf(fp, "%f\n", DataPLLBuff_BP[i].IData);
	//}
	//fclose(fp);


	int nSymbolSyncSize = 0;
	Complex* DataSymbolSyncBuff = new Complex[nSymbolSyncSize_temp];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * nSymbolSyncSize_temp);
	PM_BP.SymbolSync(DataPLLBuff_BP, DataSymbolSyncBuff, nSymbolSyncSize_temp, &nSymbolSyncSize, &SySyncBuffer);

	char* c_symbol = new char[nSymbolSyncSize];
	Judgment(DataSymbolSyncBuff, nSymbolSyncSize, c_symbol);

	//for (int i = 0; i < nSymbolSyncSize; i++)
	//{
	//	DataSymbolSyncBuff[i].QData = DataSymbolSyncBuff[i].IData;
	//}

	//FILE* fp;
	//char file_path[500];
	//sprintf(file_path, "D:/data/DataSymbolSyncBuffI.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path, "D:/data/DataSymbolSyncBuffQ.txt");
	//fp = fopen(file_path, "at");
	//for (int i = 0; i < nSymbolSyncSize; i++) {
	//	fprintf(fp, "%f\n", DataSymbolSyncBuff[i].QData);
	//}
	//fclose(fp);

		// 画星座图
	int starlen = nSymbolSyncSize;

	Complex* StarBuff = new Complex[starlen];
	memset(StarBuff, 0x00, sizeof(Complex) * starlen);

	for (int i = 0; i < starlen; i++)
	{
		StarBuff[i].IData = DataSymbolSyncBuff[i].IData;
		StarBuff[i].QData = DataSymbolSyncBuff[i].QData;
	}
		// 移动相位pi/4
	Complex* StarBuffout = new Complex[starlen];
	memset(StarBuffout, 0x00, sizeof(Complex) * starlen);
	float phase = Pi / 4;
	downphaseFrequence(StarBuff, phase, StarBuffout, starlen);




	signal_PM_BP_Out->burstData->softDistinguishData = new Complex[starlen];
	memcpy(signal_PM_BP_Out->burstData->softDistinguishData, StarBuffout, (sizeof(Complex)) * starlen);
	signal_PM_BP_Out->burstData->softDistinguishDataLen = starlen;

	//sprintf(file_path, "D:/data/nSymbolSyncSize.txt");
	//fp = fopen(file_path, "at");
	//fprintf(fp, "%d\n", nSymbolSyncSize);
	//fclose(fp);


	signal_PM_BP_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_PM_BP_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * nSymbolSyncSize);


	DELETE_ARR(DataPLLBuff_PM);
	DELETE_ARR(PLL_Freq_Part_2);
	DELETE_ARR(PLL_Phase_Part_2);
	DELETE_ARR(de_pm_bpsk);
	DELETE_ARR(de_bpsk);
	DELETE_ARR(bpsk_temp);
	DELETE_ARR(DataPLLBuff_BP);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(c_symbol);

}
//全部初始化
void Demo_PM_BP::InitialDemodulation(const DemodulationInitParamater& info) {

	m_Rb = info.rb;
	fcsub = info.fc_demo[1];  // 副载波
	fs = info.fs;
	band = info.band;

	InitBlockCFilter();

	PM_BPInit();

	InitCostasPLL();

	//InitCostasPLL_bpsk();

	InitSymbolSync();

}

void Demo_PM_BP::InitBlockCFilter() {
	// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
	//CFilter.nFilterTaps = 66;
	//CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

	// 使用数组初始化的值
	//float tempCoef[] = { 0.000300236790815230,- 0.000318083254178556,- 0.000855103154930395,- 0.000990268382954444,- 0.000487712148667362,
	//    0.000588523122424037,0.00172482575409843,0.00209416001114999,0.00104957199148652,- 0.00126316787160771,- 0.00364730576829181,
	//    - 0.00433407433765904,- 0.00211987302550528,0.00248853880247738,0.00701524733957493,0.00815335440937972,0.00390961798577038,
	//    - 0.00451151788578127,- 0.0125389944244636,- 0.0144144829075778,- 0.00686081390653839,0.00788960724885472,0.0219506333546247,
	//    0.0253950771878517,0.0122440038486131,- 0.0143814602156421,- 0.0413223332005430,- 0.0501396000048169,- 0.0259436042862064,
	//    0.0339401732374600,0.116211391424138,0.195362079644283,0.243811352622362,0.243811352622362,0.195362079644283,0.116211391424138,
	//    0.0339401732374600,- 0.0259436042862064,- 0.0501396000048169,- 0.0413223332005430,- 0.0143814602156421,0.0122440038486131,
	//    0.0253950771878517,0.0219506333546247,0.00788960724885472,- 0.00686081390653839,- 0.0144144829075778,- 0.0125389944244636,
	//    - 0.00451151788578127,0.00390961798577038,0.00815335440937972,0.00701524733957493,0.00248853880247738,- 0.00211987302550528,
	//    - 0.00433407433765904,- 0.00364730576829181,- 0.00126316787160771,0.00104957199148652,0.00209416001114999,0.00172482575409843,
	//    0.000588523122424037,- 0.000487712148667362,- 0.000990268382954444,- 0.000855103154930395,- 0.000318083254178556,0.000300236790815230 };
	//// 将初始化值复制到动态分配的数组中


	//
	//CFilter.nFilterTaps = 67;
	//CFilter.fFilterCoef = new float[CFilter.nFilterTaps];
	//float tempCoef[] =
	//{
	//	3.13548025673911e-06, 6.99986464445268e-06, 1.45333581010432e-05, 2.70282696922728e-05, 4.65659590446208e-05, 7.57197831727885e-05,
	//	0.000117567551092813, 0.000175674700347455, 0.000254042522606023,
	//	0.000357017243475641, 0.000489157830921567, 0.000655062986543612, 0.000859160728618165, 0.00110546711993688, 0.00139732380280061,
	//	0.00173712683024656, 0.00212606157431009, 0.00256386001306700, 0.00304859725078647, 0.00357654357042533, 0.00414208659041524,
	//	0.00473773521989960, 0.00535421319080164, 0.00598064519200930, 0.00660483331891592, 0.00721361601904779, 0.00779329633655395,
	//	0.00833012141863541, 0.00881079130811533, 0.00922297232038642, 0.00955578902522939, 0.00980026916212181, 0.00994971773806806,
	//	0.0100000000000000, 0.00994971773806806, 0.00980026916212181, 0.00955578902522939, 0.00922297232038642,
	//	0.00881079130811533, 0.00833012141863541, 0.00779329633655395, 0.00721361601904779, 0.00660483331891592, 0.00598064519200930, 0.00535421319080164,
	//	0.00473773521989960, 0.00414208659041524, 0.00357654357042533, 0.00304859725078647, 0.00256386001306700, 0.00212606157431009,
	//	0.00173712683024656, 0.00139732380280061, 0.00110546711993688, 0.000859160728618165, 0.000655062986543612, 0.000489157830921567,
	//	0.000357017243475641, 0.000254042522606023, 0.000175674700347455, 0.000117567551092813, 7.57197831727885e-05,
	//	4.65659590446208e-05, 2.70282696922728e-05, 1.45333581010432e-05, 6.99986464445268e-06, 3.13548025673911e-06
	//};

	CFilter.nFilterTaps = 128;
	CFilter.fFilterCoef = new float[CFilter.nFilterTaps];
	float tempCoef[] =
	{ -5.52969580960963e-05, -9.54222028105188e-05, -0.000138381990948237, -0.000185504154776396, -0.000238053040240923, -0.000297174135478286, -0.000363838647224006,
		-0.000438789265857759, -0.000522488369819510, -0.000615069905353879, -0.000716296137845490, -0.000825520406731130, -0.000941656927929731, -0.00106315857720265,
		-0.00118800345657957, -0.00131369089611776, -0.00143724737734998, -0.00155524268570631, -0.00166381641014898, -0.00175871471265129, -0.00183533709157084,
		-0.00188879266510806, -0.00191396530763820, -0.00190558678645946, -0.00185831687301206, -0.00176682924432301, -0.00162590185051045, -0.00143051030553359,
		-0.00117592276353551, -0.000857794674224752, -0.000472261769454449, -1.60296206733429e-05, 0.000513541876084603, 0.00111836068439404, 0.00179953300170700,
		0.00255730324047379, 0.00339100490195894, 0.00429902347280596, 0.00527877229620072, 0.00632668217173297, 0.00743820522401304, 0.00860783335334502,
		0.00982913134613063, 0.0110947844822200, 0.0123966602352973, 0.0137258834248010, 0.0150729239480201, 0.0164276960029479, 0.0177796675101316, 0.0191179782587701,
		0.0204315651420279, 0.0217092927119032, 0.0229400871775431, 0.0241130718946691, 0.0252177023492959, 0.0262438986271556, 0.0271821733815850, 0.0280237533669100,
		0.0287606926908178, 0.0293859760565140, 0.0298936104117514, 0.0302787035947033, 0.0305375287632947, 0.0306675736117065, 0.0306675736117065, 0.0305375287632947,
		0.0302787035947033, 0.0298936104117514, 0.0293859760565140, 0.0287606926908178, 0.0280237533669100, 0.0271821733815850, 0.0262438986271556, 0.0252177023492959,
		0.0241130718946691, 0.0229400871775431, 0.0217092927119032, 0.0204315651420279, 0.0191179782587701, 0.0177796675101316, 0.0164276960029479, 0.0150729239480201,
		0.0137258834248010, 0.0123966602352973, 0.0110947844822200, 0.00982913134613063, 0.00860783335334502, 0.00743820522401304, 0.00632668217173297, 0.00527877229620072,
		0.00429902347280596, 0.00339100490195894, 0.00255730324047379, 0.00179953300170700, 0.00111836068439404, 0.000513541876084603, -1.60296206733429e-05, -0.000472261769454449,
		-0.000857794674224752, -0.00117592276353551, -0.00143051030553359, -0.00162590185051045, -0.00176682924432301, -0.00185831687301206, -0.00190558678645946,
		-0.00191396530763820, -0.00188879266510806, -0.00183533709157084, -0.00175871471265129, -0.00166381641014898, -0.00155524268570631, -0.00143724737734998,
		-0.00131369089611776, -0.00118800345657957, -0.00106315857720265, -0.000941656927929731, -0.000825520406731130, -0.000716296137845490, -0.000615069905353879,
		-0.000522488369819510, -0.000438789265857759, -0.000363838647224006, -0.000297174135478286, -0.000238053040240923, -0.000185504154776396, -0.000138381990948237,
		-9.54222028105188e-05, -5.52969580960963e-05 };

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
	m_sAlgDemInit.nSigFs = fs;
	nDownSampClock = 1;                                         //降采样时的时钟控制  
	fDownConversionPhase = 0;                                   //下变频中的相位值
	fPLLNCO = 0;                                                //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                       //锁相环中的频率跟踪曲线
	fPLLNCO_1 = 0;                                              //锁相环中的本地NCO
	fPLLPastFreqPart_1 = 0;                                     //锁相环中的频率跟踪曲线
	fOrthDownConversionPhase = 0;
	nPLLBuffSize = 0;

	// add on 20240730
	flag = 0;
	slice = 16 * 8;
	PMBP_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[PMBP_Samplesize];


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
	//symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823/20;
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05/400;

	//symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823 / 10;
	//symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 1000;

	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 100;

	PM_BP.m_sAlgDemInit.nSampPerSymb = 4;                       //每个码元中的采样点数

	PM_BP.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_PM_BP::InitCostasPLL() {
	// 用于PM解调
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;     //环路阻尼系数
	int Ko = 1;              //压控振荡器增益
	int Kd = 1;              //鉴相器增益
	int K = Ko * Kd;
	float fBLcoef = 0.01;

	// 第一段
	BL = fBLcoef * band;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第二段
	fBLcoef = 0.01;
	BL = fBLcoef * band;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll2.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll2.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第三段
	fBLcoef = 0.01;
	BL = fBLcoef * band;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll3.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll3.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
}

void Demo_PM_BP::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	for (int i = 0; i < bufflen; i++)
	{
		if (buffin[i].IData > 0)
		{
			buffout[i] = 1;
		}
		else
		{
			buffout[i] = 0;
		}
	}
}

void Demo_PM_BP::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart, float* PLL_Freq_Part, float* PLL_Phase_Part)
{
	float discriminator_out = 0;
	float pll_phase_part;
	float freq_control;
	float past_pll_freq_part;
	float past_nco_phase;
	float PLLCoef1 = CostasPll.fPLLLoopFilterCoef1;
	float PLLCoef2 = CostasPll.fPLLLoopFilterCoef2;
	float PLLCoef21 = CostasPll2.fPLLLoopFilterCoef1;
	float PLLCoef22 = CostasPll2.fPLLLoopFilterCoef2;
	float PLLCoef31 = CostasPll3.fPLLLoopFilterCoef1;
	float PLLCoef32 = CostasPll3.fPLLLoopFilterCoef2;

	for (int i = 0; i < nBuffSize; i++)
	{
		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fPLLNCO), sin(*fPLLNCO), BuffOut[i].IData, BuffOut[i].QData);
		discriminator_out = atan2(BuffOut[i].QData, BuffOut[i].IData);

		if (i < 1000)
		{
			pll_phase_part = discriminator_out * PLLCoef1;
			PLL_Phase_Part[i] = pll_phase_part;
			past_pll_freq_part = *fPLLPastFreqPart;
			*fPLLPastFreqPart = discriminator_out * PLLCoef2 + past_pll_freq_part;
			PLL_Freq_Part[i] = *fPLLPastFreqPart;
			freq_control = pll_phase_part + *fPLLPastFreqPart;
			past_nco_phase = *fPLLNCO;
			*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
			if (*fPLLNCO > 2 * Pi)
			{
				*fPLLNCO -= 2 * Pi;
			}
		}
		else
		{
			pll_phase_part = discriminator_out * PLLCoef21;
			PLL_Phase_Part[i] = pll_phase_part;
			past_pll_freq_part = *fPLLPastFreqPart;
			*fPLLPastFreqPart = discriminator_out * PLLCoef22 + past_pll_freq_part;
			PLL_Freq_Part[i] = *fPLLPastFreqPart;
			freq_control = pll_phase_part + *fPLLPastFreqPart;
			past_nco_phase = *fPLLNCO;
			*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
			if (*fPLLNCO > 2 * Pi)
			{
				*fPLLNCO -= 2 * Pi;
			}
		}
	}
}

void Demo_PM_BP::InitCostasPLL_bpsk()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;     //环路阻尼系数
	int Ko = 1;              //压控振荡器增益
	int Kd = 1;              //鉴相器增益
	int K = Ko * Kd;
	float fBLcoef = 0.02;

	// 第一段
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);
	CostasPll_bpsk.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll_bpsk.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第二段
	fBLcoef = 0.002;
	BL = fBLcoef * m_Rb;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);
	CostasPll2_bpsk.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll2_bpsk.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第三段
	fBLcoef = 0.002;
	BL = fBLcoef * m_Rb;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);
	CostasPll3_bpsk.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll3_bpsk.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

}

void Demo_PM_BP::PLLCostas_bpsk(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float* fPLLNCO, float* fPLLPastFreqPart)
{
	float discriminator_out = 0;
	float pll_phase_part;
	float freq_control;
	float past_pll_freq_part;
	float past_nco_phase;
	int nPLLSignI = 0;
	int nPLLSignQ = 0;
	float PLLCoef1 = CostasPll_bpsk.fPLLLoopFilterCoef1;
	float PLLCoef2 = CostasPll_bpsk.fPLLLoopFilterCoef2;
	float PLLCoef21 = CostasPll2_bpsk.fPLLLoopFilterCoef1;
	float PLLCoef22 = CostasPll2_bpsk.fPLLLoopFilterCoef2;
	float PLLCoef31 = CostasPll3_bpsk.fPLLLoopFilterCoef1;
	float PLLCoef32 = CostasPll3_bpsk.fPLLLoopFilterCoef2;

	//QString path = QString::fromLocal8Bit("C:\\outputfile\\PLLfreqpart.dat");
	//QFile newFile(path);
	//if (!newFile.open(QIODevice::ReadWrite | QIODevice::Append)) {

	//    return;
	//}

	for (int i = 0; i < nBuffSize; i++)
	{
		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fPLLNCO), sin(*fPLLNCO), BuffOut[i].IData, BuffOut[i].QData);
		SIGN(BuffOut[i].IData, nPLLSignI);
		SIGN(BuffOut[i].QData, nPLLSignQ);

		discriminator_out = nPLLSignI * BuffOut[i].QData;

		if (i < 1000)
		{
			pll_phase_part = discriminator_out * PLLCoef1;
			past_pll_freq_part = *fPLLPastFreqPart;
			*fPLLPastFreqPart = discriminator_out * PLLCoef2 + past_pll_freq_part;
			freq_control = pll_phase_part + *fPLLPastFreqPart;
			past_nco_phase = *fPLLNCO;
			*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
			if (*fPLLNCO > 2 * Pi)
			{
				*fPLLNCO -= 2 * Pi;
			}
		}
		else if (i < 5000)
		{
			pll_phase_part = discriminator_out * PLLCoef21;
			past_pll_freq_part = *fPLLPastFreqPart;
			*fPLLPastFreqPart = discriminator_out * PLLCoef22 + past_pll_freq_part;
			freq_control = pll_phase_part + *fPLLPastFreqPart;
			past_nco_phase = *fPLLNCO;
			*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
			if (*fPLLNCO > 2 * Pi)
			{
				*fPLLNCO -= 2 * Pi;
			}
		}
		else
		{
			pll_phase_part = discriminator_out * PLLCoef31;
			past_pll_freq_part = *fPLLPastFreqPart;
			*fPLLPastFreqPart = discriminator_out * PLLCoef32 + past_pll_freq_part;
			freq_control = pll_phase_part + *fPLLPastFreqPart;
			past_nco_phase = *fPLLNCO;
			*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
			if (*fPLLNCO > 2 * Pi)
			{
				*fPLLNCO -= 2 * Pi;
			}
		}



		//newFile.write((char*)fPLLPastFreqPart, 4);
	}

	//newFile.close();
}

void Demo_PM_BP::downphase(float BuffIn[], Complex BuffOut[], int nBuffSize, float* fDownConversionFc, double* fOrthDownConversionPhase)
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

// add on 20240730
void Demo_PM_BP::est_pmbp_rb(Complex* data, int len, int fftlen, float fs, float* est_rb)
{
	float* baoluo = new float[len];
	memset(baoluo, 0, sizeof(float) * len);
	for (int i = 0; i < len; i++)
	{
		baoluo[i] = data[i].IData * data[i].IData + data[i].QData * data[i].QData;
	}

	float* diff_baoluo = new float[len];
	memset(diff_baoluo, 0, sizeof(float) * len);

	for (int i = 1; i < len; i++)
	{
		diff_baoluo[i - 1] = baoluo[i] - baoluo[i - 1];
	}

	int datalen_temp = fftlen > len ? len : fftlen;

	float* data_I_fft = new float[fftlen];
	memset(data_I_fft, 0, sizeof(float) * fftlen);
	float* data_Q_fft = new float[fftlen];
	memset(data_Q_fft, 0, sizeof(float) * fftlen);

	memcpy(data_I_fft, diff_baoluo, sizeof(float) * datalen_temp);
	float* fft_spectrum = new float[fftlen];
	memset(fft_spectrum, 0, fftlen * sizeof(float));
	ippfft(data_I_fft, data_Q_fft, fftlen, fft_spectrum);

	// 消除无关谱线的影响
	for (int i = 0; i < 1000; i++)
	{
		fft_spectrum[i] = 0;
	}

	int max_index_rb = findmax1(fft_spectrum, 2, fftlen, 0);


	int datalen = max_index_rb + 5000;
	float* fftdata1 = new float[datalen];
	memcpy(fftdata1, fft_spectrum, sizeof(float) * datalen);

	// 排序找索引
	int* index = new int[datalen];
	dgz_sort_pmbpqp(fftdata1, index, datalen);

	// 最高峰
	int index_v[10] = { 0 };
	memcpy(index_v, index + datalen - 10, sizeof(int) * 10);
	int index_v_temp[10] = { 0 };
	int index_flag = 0;
	for (int i = 9; i > 0; i--)
	{
		index_v_temp[index_flag] = index_v[i] - index_v[i - 1];
		if (index_v_temp[index_flag] < 10)
		{
			index_v_temp[index_flag] = 0;
			continue;
		}
		index_flag++;
	}

	sort(index_v_temp, index_v_temp + index_flag);

	int difftemp = 0;

	difftemp = index_v_temp[0];  // 取最小的谱线之差

	*est_rb = difftemp * fs / fftlen;


	delete[] fftdata1;

}

void Demo_PM_BP::ippfft(float* data_real, float* data_imag, int fftLength, float* output)
{
	Ipp32fc* pDst = NULL;
	Ipp32fc* pSrc = NULL;
	int FFT_size = fftLength;
	int FFTOrder = (log(FFT_size) / log(2)); //add by zhuxue
	IppsFFTSpec_C_32fc* pSpec = 0;
	Ipp8u* pMemSpec = 0;
	Ipp8u* pMemInit = 0;
	Ipp8u* pMemBuffer = 0;
	int sizeSpec = 0;
	int sizeInit = 0;
	int sizeBuffer = 0;
	int flag = IPP_FFT_NODIV_BY_ANY;
	int sizeFft = (int)FFT_size;
	//add by zhuxue
	pSrc = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * sizeFft);
	pDst = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * sizeFft);
	std::memset(pSrc, 0, sizeof(Ipp32fc) * sizeFft);
	std::memset(pDst, 0, sizeof(Ipp32fc) * sizeFft);
	for (int i = 0; i < sizeFft; i++)
	{
		pSrc[i].re = data_real[i];
		pSrc[i].im = data_imag[i];
	}
	/// get sizes for required buffers
	ippsFFTGetSize_C_32fc(FFTOrder, flag, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuffer);
	//	printf("sizeSpec=%d,sizeInit=%d,sizeBuffer=%d\n", sizeSpec, sizeInit, sizeBuffer);
	/// allocate memory for required buffers
	pMemSpec = (Ipp8u*)ippMalloc(sizeSpec);
	if (sizeInit > 0)
	{
		pMemInit = (Ipp8u*)ippMalloc(sizeInit);
	}
	if (sizeBuffer > 0)
	{
		pMemBuffer = (Ipp8u*)ippMalloc(sizeBuffer);
	}
	/// initialize FFT specification structure
	ippsFFTInit_C_32fc(&pSpec, FFTOrder, flag, ippAlgHintNone, pMemSpec, pMemInit);
	/// free initialization buffer
	if (sizeInit > 0)
	{
		ippFree(pMemInit);
	}
	/// perform forward FFT
	ippsFFTFwd_CToC_32fc(pSrc, pDst, pSpec, pMemBuffer);
	for (int n = 0; n < sizeFft; n++)
	{
		output[n] = sqrt((float)pDst[n].re * (float)pDst[n].re + (float)pDst[n].im * (float)pDst[n].im);   //这个输出存放的是复信号的模值
		//Ioutput[n] = (float)pDst[n].re;
		//Qoutput[n] = (float)pDst[n].im;
		/*if (flag1 == -1)
		{
			Qoutput[n] = -Qoutput[n];
		}*/
	}
	/// ...
	/// free buffers
	if (sizeBuffer > 0)
	{
		ippFree(pMemBuffer);
	}
	ippFree(pMemSpec);
	ippFree(pDst);
	ippFree(pSrc);
}

int Demo_PM_BP::findmax1(float* a, int p, float FFT_size, int start)
{
	int i, max = start;
	for (i = start; i < FFT_size / p; i++)
	{
		if (a[max] < a[i])
		{
			max = i;
		}
	}
	return max;
}

// 输入数据input
// 输出数据output
// 采样率变换倍数step = fs_after / fs_before
// 输入信号长度lengthin
// 输出信号长度lengthout

void Demo_PM_BP::Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout)
{
	int i = 0;
	int j = 0;
	int k = 0;
	static double step_temp = 0;
	//double step_temp = 0;
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

void Demo_PM_BP::dgz_sort_pmbpqp(float arr[], int index[], int n)
{
	int i, j, min_idx;
	float temp;

	// 初始化索引数组
	for (i = 0; i < n; i++) {
		index[i] = i;
	}

	// 选择排序算法
	for (i = 0; i < n - 1; i++) {
		min_idx = i;
		for (j = i + 1; j < n; j++) {
			if (arr[j] < arr[min_idx]) {
				min_idx = j;
			}
		}

		// 交换最小元素到前面
		temp = arr[min_idx];
		arr[min_idx] = arr[i];
		arr[i] = temp;

		// 同时交换索引
		int temp_idx = index[min_idx];
		index[min_idx] = index[i];
		index[i] = temp_idx;
	}
}

void Demo_PM_BP::downphaseFrequence(Complex* data_input, float phase, Complex* downout, int dataLength)
{
	//两路
	for (int i = 0; i < dataLength; i++)
	{
		downout[i].IData = data_input[i].IData * cos(phase) - data_input[i].QData * sin(phase);
		downout[i].QData = data_input[i].QData * cos(phase) + data_input[i].IData * sin(phase);
	}
}
