#include "SignalUQPSK.h"
Demo_UQPSK::Demo_UQPSK(int index) :m_index(index)
{

}

Demo_UQPSK::~Demo_UQPSK() {

}

bool Demo_UQPSK::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter);
	return 1;
}

void Demo_UQPSK::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
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

void Demo_UQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_UQPSK_Out) {

	if (nSilceCount == 1)
	{
		// 第二次进入该函数，载波系数保持不变
		// 如果使用两段
		CostasPll = CostasPll2;

		//// 如果使用三段
		//CostasPll = CostasPll2;
		//CostasPll2 = CostasPll3;

	}
	nSilceCount++;

	De_UQPSK.AGC(data_input_slice, fAGCPastVc, UQPSK_Samplesize);


	//QString path4 = QString::fromLocal8Bit("C:\\outputfile\\Data.dat");
	//QFile newFile4(path4);
	//newFile4.open(QIODevice::ReadWrite | QIODevice::Append);
	//QDataStream outstm4(&newFile4);
	//outstm4.writeRawData((char*)data_input_slice, UQPSK_Samplesize * 2 * 4);
	//newFile4.close();

	Complex* DataPLLBuff = new Complex[UQPSK_Samplesize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * UQPSK_Samplesize);

	// add on 20240726
	PLLCostas(data_input_slice, DataPLLBuff, UQPSK_Samplesize, 2, &fPLLNCO, &fPLLPastFreqPart);

	//QString path = QString::fromLocal8Bit("C:\\outputfile\\DataPLLBuff.dat");
	//QFile newFile(path);
	//newFile.open(QIODevice::ReadWrite | QIODevice::Append);
	//QDataStream outstm(&newFile);
	//outstm.writeRawData((char*)DataPLLBuff, UQPSK_Samplesize * 2 * 4);
	//newFile.close();

	// 由于载波同步会相位翻转，所以要
	// 判断载波同步后的I路还是Q路对应的码速率、以进行码元同步
	int fftlen = 16384;


	//QString path3 = QString::fromLocal8Bit("C:\\outputfile\\ffti.dat");
	//QFile newFile3(path3);
	//newFile3.open(QIODevice::ReadWrite | QIODevice::Append);
	//QDataStream outstm3(&newFile3);

	if (flag1 == 0)
	{
		// UQPSK_Samplesize = 81920 = 16384*5
		float* datai = new float[fftlen];
		float* dataq = new float[fftlen];
		float* zero_temp = new float[fftlen];
		float* fftouti = new float[fftlen];
		float* fftoutq = new float[fftlen];
		float* dataiest = new float[fftlen];
		float* dataqest = new float[fftlen];
		float v_i[5] = { 0 };
		float v_q[5] = { 0 };
		int halffftlen = fftlen / 2;
		memset(zero_temp, 0x00, sizeof(float) * fftlen);
		for (int j = 0; j < 5; j++)
		{
			memset(datai, 0x00, sizeof(float) * fftlen);
			for(int i = 0;i< fftlen;i++)
			{
				datai[i] = DataPLLBuff[i + j * fftlen].IData* DataPLLBuff[i + j * fftlen].IData;
				dataq[i] = DataPLLBuff[i + j * fftlen].QData* DataPLLBuff[i + j * fftlen].QData;
			}

			ippfft(datai, zero_temp, fftlen, fftouti);
			memcpy(dataiest, fftouti,sizeof(float) * halffftlen);

			// 归一化，为了写dat文件
			int max_index = getFFTmax(dataiest, 0, halffftlen);
			int max_index2 = getFFTmax(dataqest, 0, halffftlen);
			int temp = dataiest[max_index];
			int temp2 = dataqest[max_index2];
			for(int i =0;i<halffftlen;i++)
			{
				dataiest[i] /= temp;
				dataqest[i] /= temp2;
			}
			//outstm3.writeRawData((char*)dataiest, halffftlen * 4);


			ippfft(dataq, zero_temp, fftlen, fftoutq);
			memcpy(dataqest, fftoutq, sizeof(float) * halffftlen);

			float rb1_I = 0;
			float rb2_I = 0;
			float rb1_Q = 0;
			float rb2_Q = 0;

			estuqpsk(dataiest, halffftlen, &rb1_I, &rb2_I, fs, fftlen);
			estuqpsk(dataqest, halffftlen, &rb1_Q, &rb2_Q, fs, fftlen);

			v_i[j] = rb1_I > rb2_I ? rb1_I : rb2_I;
			v_q[j] = rb1_Q > rb2_Q ? rb1_Q : rb2_Q;

		}

		// 选择
		sort(v_i, v_i + 5);
		sort(v_q, v_q + 5);
		// 去除最大值和最小值    add on 20240728
		float v_i_temp[3] = { v_i[1],v_i[2],v_i[3] };
		float v_q_temp[3] = { v_q[1],v_q[2],v_q[3] };

		int sum_i = accumulate(v_i_temp, v_i_temp + 3, 0);
		int sum_q = accumulate(v_q_temp, v_q_temp + 3, 0);
		

		if (sum_i > sum_q)
		{
			// 载波同步后的I路是大速率
			rbEst1 = m_Rb > m_RbQ ? m_Rb : m_RbQ;
			rbEst2 = (rbEst1 == m_Rb) ? m_RbQ : m_Rb;
		}
		else
		{
			// 载波同步后的I路是小速率
			rbEst1 = m_Rb < m_RbQ ? m_Rb : m_RbQ;
			rbEst2 = (rbEst1 == m_Rb) ? m_RbQ : m_Rb;
		}

		flag1 = 1;
		DELETE_ARR(datai);
		DELETE_ARR(dataq);
		DELETE_ARR(zero_temp);
		DELETE_ARR(fftouti);
		DELETE_ARR(fftoutq);
		DELETE_ARR(dataiest);
		DELETE_ARR(dataqest);
	}
	//newFile3.close();
	// 交换载波同步后I路和Q路数据
	// 只当作I路为大速率，Q路为小速率，方便后续处理
	if (rbEst1 < rbEst2)
	{
		float* temp = new float[UQPSK_Samplesize];
		for (int i = 0; i < UQPSK_Samplesize; i++) {
			temp[i] = DataPLLBuff[i].QData;
			DataPLLBuff[i].QData = DataPLLBuff[i].IData;
			DataPLLBuff[i].IData = temp[i];
		}
		DELETE_ARR(temp);
	}

	// I 路码元同步
	int nSymbolSyncSize = 0;//突发在当前片中的长度

	Complex* DataPLLBuff_Itemp = new Complex[UQPSK_Samplesize];
	memset(DataPLLBuff_Itemp, 0x00, sizeof(Complex) * UQPSK_Samplesize);
	for(int i =0;i< UQPSK_Samplesize;i++)
	{
		DataPLLBuff_Itemp[i].IData = DataPLLBuff[i].IData;
		DataPLLBuff_Itemp[i].QData = DataPLLBuff[i].IData;
	}

	Complex* DataSymbolSyncBuff = new Complex[UQPSK_Samplesize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * UQPSK_Samplesize);
	De_UQPSK.SymbolSync(DataPLLBuff_Itemp, DataSymbolSyncBuff, UQPSK_Samplesize, &nSymbolSyncSize, &SySyncBufferI);
	char* DataResultI = new char[nSymbolSyncSize];
	memset(DataResultI, 0x00, sizeof(char) * nSymbolSyncSize);
	Judgment(DataSymbolSyncBuff, nSymbolSyncSize, DataResultI);

	//QString path1 = QString::fromLocal8Bit("C:\\outputfile\\DataSymbolSyncBuffI.dat");
	//QFile newFile1(path1);
	//newFile1.open(QIODevice::ReadWrite | QIODevice::Append);
	//QDataStream outstm1(&newFile1);
	//outstm1.writeRawData((char*)DataSymbolSyncBuff, nSymbolSyncSize * 2 * 4);
	//newFile1.close();


	// Q 路码元同步 add on 20240728
	Complex* DataSymbolSyncBuffQ_temp = new Complex[UQPSK_Samplesize];
	memset(DataSymbolSyncBuffQ_temp, 0x00, sizeof(Complex) * UQPSK_Samplesize);

	for (int i = 0; i < UQPSK_Samplesize; i++)
	{
		DataSymbolSyncBuffQ_temp[i].IData = DataPLLBuff[i].QData;
		DataSymbolSyncBuffQ_temp[i].QData = DataPLLBuff[i].QData;
	}

	Complex* DataPLLBuff_Qtemp = new Complex[UQPSK_Samplesize];
	memset(DataPLLBuff_Qtemp, 0x00, sizeof(Complex) * UQPSK_Samplesize);

	int rb_temp = m_Rb > m_RbQ ? m_RbQ : m_Rb;
	double step = ((double)(4 * rb_temp)) / ((double)fs);
	int nSymbolSyncSizeQ_temp = 0;	       //降采样之后长度
	Downfs(DataSymbolSyncBuffQ_temp, DataPLLBuff_Qtemp, step, UQPSK_Samplesize, &nSymbolSyncSizeQ_temp);

	int nSymbolSyncSizeQ = 0;	       //Q路符号同步长度
	Complex* DataSymbolSyncBuffQ = new Complex[nSymbolSyncSizeQ_temp];
	memset(DataSymbolSyncBuffQ, 0x00, sizeof(Complex) * nSymbolSyncSizeQ_temp);

	De_UQPSK2.SymbolSync(DataPLLBuff_Qtemp, DataSymbolSyncBuffQ, nSymbolSyncSizeQ_temp, &nSymbolSyncSizeQ, &SySyncBufferQ);

	//QString path2 = QString::fromLocal8Bit("C:\\outputfile\\DataSymbolSyncBuffQ.dat");
	//QFile newFile2(path2);
	//newFile2.open(QIODevice::ReadWrite | QIODevice::Append);
	//QDataStream outstm2(&newFile2);
	//outstm2.writeRawData((char*)DataSymbolSyncBuffQ, nSymbolSyncSizeQ * 2 * 4);
	//newFile2.close();

	signal_UQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSizeQ];

	for (int i = 0; i < nSymbolSyncSizeQ; i++)
	{
		signal_UQPSK_Out->burstData->softDistinguishData[i].IData = DataSymbolSyncBuff[i].IData;
		signal_UQPSK_Out->burstData->softDistinguishData[i].QData = DataSymbolSyncBuffQ[i].IData;
	}

	signal_UQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
	//signal_UQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSizeQ;
	signal_UQPSK_Out->burstData->softDistinguishDataLenQ = nSymbolSyncSizeQ;
	char* DataResultQ = new char[nSymbolSyncSizeQ];
	memset(DataResultQ, 0x00, sizeof(char) * nSymbolSyncSizeQ);
	Judgment(DataSymbolSyncBuffQ, nSymbolSyncSizeQ, DataResultQ);

	signal_UQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_UQPSK_Out->burstData->nDemodulationByteI, DataResultI, (sizeof(char)) * nSymbolSyncSize);
	signal_UQPSK_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSizeQ];
	memcpy(signal_UQPSK_Out->burstData->nDemodulationByteQ, DataResultQ, (sizeof(char)) * nSymbolSyncSizeQ);

	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataResultI);
	DELETE_ARR(DataSymbolSyncBuffQ);
	DELETE_ARR(DataSymbolSyncBuffQ_temp);
	DELETE_ARR(DataResultQ);
	DELETE_ARR(DataSymbolSyncBuff);
	DELETE_ARR(DataPLLBuff_Itemp);
	DELETE_ARR(DataPLLBuff_Qtemp);
}

void Demo_UQPSK::InitBlockCFilter()
{

	// 假设CFilter是一个已定义的结构体或类，并且nFilterTaps和fFilterCoef是其成员
	CBlockFilter.nFilterTaps = 16;
	CBlockFilter.fFilterCoef = new float[CBlockFilter.nFilterTaps];

	// 使用数组初始化的值
	float tempCoef[] = { -0.0132465027355890, 8.32286573684456e-18, 0.0362191880550484, -1.37671571265597e-17, -0.0878600226999023, 1.80033486821016e-17, 0.313453567624625, 0.502867539511636, 0.313453567624625, 1.80033486821016e-17, -0.0878600226999023, -1.37671571265597e-17, 0.0362191880550484, 8.32286573684456e-18, -0.0132465027355890, -3.33169527851874e-18 };
	// 将初始化值复制到动态分配的数组中
	for (int i = 0; i < CFilter.nFilterTaps; ++i) {
		CBlockFilter.fFilterCoef[i] = tempCoef[i];
	}
	cBlockFilterBuff = new Complex[CBlockFilter.nFilterTaps];
	memset(cBlockFilterBuff, 0x00, sizeof(Complex) * CBlockFilter.nFilterTaps);
	CBlockFilter.bCoefEvenSym = true;
	De_UQPSK.CBlockFilter = CBlockFilter;

}

void Demo_UQPSK::InitialDemodulation(const DemodulationInitParamater& info) {

	m_Rb = info.rb;
	m_RbQ = info.rb1;
	rbEst1 = m_Rb;
	rbEst2 = m_RbQ;

	for (int i = 0; i < 2; i++)
	{
		fc[i] = info.fc_demo[i];
	}
	fs = info.fs;
	InitBlockCFilter();
	UQPSKInit();
	InitCostasPLL();
	InitSymbolSync();
	InitSRC();
}

void Demo_UQPSK::UQPSKInit()
{
	m_sAlgDemInit.llSliceLength = 8192;
	m_SamleSize = m_sAlgDemInit.llSliceLength;
	De_UQPSK.m_sAlgDemInit.llSliceLength = m_sAlgDemInit.llSliceLength;

	nDownSampClock = 1;                                         //降采样时的时钟控制  
	fDownConversionPhase = 0;                                   //下变频中的相位值
	fPLLNCO = 0;                                                //锁相环中的本地NCO
	fPLLPastFreqPart = 0;                                       //锁相环中的频率跟踪曲线
	nPLLBuffSize = 0;


	// add on 20240621
	flag = 0;
	slice = 10;
	UQPSK_Samplesize = slice * m_SamleSize;
	Databuff0 = new Complex[UQPSK_Samplesize];

	//******************符号同步I路初始化***********
	SySyncBufferI.CSymbolSyncBuff = new Complex[4];
	SySyncBufferI.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBufferI.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBufferI.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBufferI.fSymbolSyncW = 0.5;                           //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBufferI.fSymbolSyncN = 0.9;                           //符号同步NCO寄存器，初值设为1
	SySyncBufferI.fSymbolSyncNTemp = 0.9;                       //符号同步NCO暂时的寄存器，初值设为1
	SySyncBufferI.fSymbolSyncU = 0.6;                           //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBufferI.nSymbolSyncKK = 0;                            //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBufferI.fSymbolSyncPastW = 0;                         //符号同步W的缓存值
	SySyncBufferI.fSymbolSyncPastN = 0;                         //符号同步N的缓存值
	SySyncBufferI.fSymbolSyncTimeError1 = 0;                    //符号同步time_error缓存值
	SySyncBufferI.fSymbolSyncTimeError2 = 0;                    //符号同步time_error缓存值

	//******************符号同步Q路初始化***********
	SySyncBufferQ.CSymbolSyncBuff = new Complex[4];
	SySyncBufferQ.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBufferQ.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBufferQ.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBufferQ.fSymbolSyncW = 0.5;                           //符号同步环路滤波器输出寄存器，初值设为0.5
	SySyncBufferQ.fSymbolSyncN = 0.9;                           //符号同步NCO寄存器，初值设为1
	SySyncBufferQ.fSymbolSyncNTemp = 0.9;                       //符号同步NCO暂时的寄存器，初值设为1
	SySyncBufferQ.fSymbolSyncU = 0.6;                           //符号同步NCO输出的定时分数间隔寄存器，初值设为0.6
	SySyncBufferQ.nSymbolSyncKK = 0;                            //符号同步用来表示Ti时间序号,指示u,yI,yQ
	SySyncBufferQ.fSymbolSyncPastW = 0;                         //符号同步W的缓存值
	SySyncBufferQ.fSymbolSyncPastN = 0;                         //符号同步N的缓存值
	SySyncBufferQ.fSymbolSyncTimeError1 = 0;                    //符号同步time_error缓存值
	SySyncBufferQ.fSymbolSyncTimeError2 = 0;                    //符号同步time_error缓存值

}

void Demo_UQPSK::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823 / 10;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 1000;
	// add on 20240727
	symbolsyncsactorInit2.fSymbolSyncFactor1 = 0.013331990840823;
	symbolsyncsactorInit2.fSymbolSyncFactor2 = 8.889783703667417e-05/100;

	De_UQPSK.m_sAlgDemInit.nSampPerSymb = 4;
	De_UQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
	// add on 20240727
	De_UQPSK2.m_sAlgDemInit.nSampPerSymb = 4;
	De_UQPSK2.symbolsyncsactorInit = symbolsyncsactorInit2;


}

void Demo_UQPSK::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;     //环路阻尼系数
	int Ko = 1;              //压控振荡器增益
	int Kd = 1;              //鉴相器增益
	int K = Ko * Kd;
	float fBLcoef = 0.02;
	float band = m_Rb > m_RbQ ? m_Rb : m_RbQ;

	// 第一段
	BL = fBLcoef * band;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第二段
	fBLcoef = 0.002;
	BL = fBLcoef * band;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll2.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll2.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 第三段
	fBLcoef = 0.002;
	BL = fBLcoef * band;                                        // add on 20240724
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)band);
	CostasPll3.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll3.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;

	// 不调用其他类了，所以下面语句不需要
	//De_UQPSK.CostasPll = CostasPll;  
}

// add on 20240726
void Demo_UQPSK::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart)
{
	float discriminator_out = 0;
	float pll_phase_part;
	float freq_control;
	float past_pll_freq_part;
	float past_nco_phase;
	int nPLLSignI = 0;
	int nPLLSignQ = 0;
	float PLLCoef1 = CostasPll.fPLLLoopFilterCoef1;
	float PLLCoef2 = CostasPll.fPLLLoopFilterCoef2;
	float PLLCoef21 = CostasPll2.fPLLLoopFilterCoef1;
	float PLLCoef22 = CostasPll2.fPLLLoopFilterCoef2;
	float PLLCoef31 = CostasPll3.fPLLLoopFilterCoef1;
	float PLLCoef32 = CostasPll3.fPLLLoopFilterCoef2;

	//QString path = QString::fromLocal8Bit("C:\\outputfile\\PLLfreqpart.dat");
	//QFile newFile(path);
	//if (!newFile.open(QIODevice::ReadWrite | QIODevice::Append)) {

	//	return;
	//}

	for (int i = 0; i < nBuffSize; i++)
	{
		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fPLLNCO), sin(*fPLLNCO), BuffOut[i].IData, BuffOut[i].QData);
		SIGN(BuffOut[i].IData, nPLLSignI);
		SIGN(BuffOut[i].QData, nPLLSignQ);
		switch (nSignalType)
		{
		case 1: 											//BPSK/PM_BP
		{
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
			else
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
			break;
		}
		case 2:
		{
			discriminator_out = nPLLSignI * BuffOut[i].QData - nPLLSignQ * BuffOut[i].IData;             //QPSK/PM_QP/UQPSK
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
			else
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
			//newFile.write((char*)fPLLPastFreqPart, 4);
			break;
		}

		}
	}
	//newFile.close();
}


void Demo_UQPSK::Judgment(Complex buffin[], int bufflen, char buffout[])
{
	Complex* Judgm = new Complex[bufflen];
	for (int i = 0; i < bufflen; i++) {
		if (buffin[i].IData > 0) {
			buffout[i] = 1;
		}
		else {
			buffout[i] = 0;
		}
	}
	delete(Judgm);
}

void Demo_UQPSK::InitSRC()
{
	// 该降采样函数貌似只适用整数倍
	SRCParam.bUpSampFilterValid = false;

	SRCParam.nUpSampRate = 1;
	SRCParam.nDownSampRate = m_Rb / m_RbQ;
	//SRCParam.nDownSampRate = ceil((float)(m_Rb > m_RbQ) ? (m_Rb / m_RbQ) : (m_RbQ / m_Rb));   // modify on 20240727
	//SRCParam.nDownSampRate =2;   // modify on 20240727
	De_UQPSK.SRCParam = SRCParam;
}

void Demo_UQPSK::ippfft(float* data_real, float* data_imag, int fftLength, float* output)
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

void Demo_UQPSK::fft_shift(float* fft_s, int FFT_size)
{
	int i;
	float temp;
	for (i = 0; i < FFT_size / 2; i++)
	{
		temp = fft_s[i];
		fft_s[i] = fft_s[i + FFT_size / 2];
		fft_s[i + FFT_size / 2] = temp;
	}
	return;
}


int Demo_UQPSK::getFFTmax(float* a, int start, int end)
{
	int i, max = start;
	for (i = start; i < end; i++)
	{
		if (a[max] < a[i])
		{
			max = i;
		}
	}
	return max;
}

void Demo_UQPSK::estuqpsk(float* s_2_real_temp, int len_temp, float* rb1, float* rb2, int fs, int fftLength)
{
	int max_index = getFFTmax(s_2_real_temp, 0, len_temp);
	// 去除最高峰周围影响
	for (int i = 0; i < max_index + 20; i++)              // 20 可随调试修改
	{
		s_2_real_temp[i] = 0;
	}

	// 发现后20个点左右会有峰值影响 add on 20240728
	for (int i = len_temp-20; i < len_temp; i++)              // 20 可随调试修改
	{
		s_2_real_temp[i] = 0;
	}

	// 查找大速率
	int max_index2 = getFFTmax(s_2_real_temp, 0, len_temp);
	//estInfo.rbEst1 = (float)(max_index2) * fs / fftLength;
	*rb1 = (float)(max_index2)*fs / fftLength;

	// 查找小速率
	int len_temp2 = max_index2 * 2 - 20;				// 20 可随调试修改
	len_temp2 = len_temp2 > len_temp ? len_temp : len_temp2;

	float* s_2_real_temp2 = new float[len_temp2];
	memcpy(s_2_real_temp2, s_2_real_temp, sizeof(float) * len_temp2);

	// 去除一些影响
	for (int i = max_index2 - 20; i < max_index2 + 20; i++)
	{
		s_2_real_temp2[i] = 0;
	}

	int max_index3 = getFFTmax(s_2_real_temp2, 0, len_temp2);

	// 求均值
	float mean_ave = mean(s_2_real_temp2, 0, len_temp2);

	// 设置阈值
	float thread3 = mean_ave * 4.5;  // 6 可随调试修改

	// 做选择
	*rb2 = *rb1;								// 相同码速
	if (s_2_real_temp2[max_index3] > thread3)
	{
		*rb2 = (float)max_index3 * fs / fftLength;  // 不同码速
	}

	DELETE_ARR(s_2_real_temp2);
}

float  Demo_UQPSK::mean(float* data, int start, int end)
{
	int i = 0;
	float mean, sum = 0;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	mean = sum / (end - start);
	return mean;
}

// 输入数据input
// 输出数据output
// 采样率变换倍数step = fs_after / fs_before
// 输入信号长度lengthin
// 输出信号长度lengthout
void Demo_UQPSK::Downfs(Complex* input, Complex* output, double step, int lengthin, int* lengthout)
{
	int i = 0;
	int j = 0;
	int k = 0;
	static double step_temp = 0;   // 适用在线版本
	//double step_temp = 0;		     // 适用离线版本
	
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