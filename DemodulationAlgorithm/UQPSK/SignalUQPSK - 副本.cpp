#include "SignalUQPSK.h"
#include  <QDataStream>
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
	Demodulation(dataInputRealSlice, baseTime, demodulationResult);
}

void Demo_UQPSK::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_UQPSK_Out) {


	De_UQPSK.AGC(data_input_slice, fAGCPastVc);

	QString path = QString::fromLocal8Bit("D:\\data\\data_input_slice.dat");
	QFile newFile(path);
	newFile.open(QIODevice::ReadWrite | QIODevice::Append);
	QDataStream outstm(&newFile);
	outstm.writeRawData((char*)data_input_slice, m_SamleSize * 2 * 4);
	newFile.close();



	Complex* DataPLLBuff = new Complex[m_SamleSize];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * m_SamleSize);
	Complex* DataFLLBuff = new Complex[m_SamleSize];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_UQPSK.CarrierSync(data_input_slice, DataPLLBuff, DataFLLBuff, m_SamleSize, 1, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);



	if (flag == 0) {
		float* signCiFang_I = new float[m_SamleSize];
		memset(signCiFang_I, 0, sizeof(float) * m_SamleSize);
		float* signCiFang_Q = new float[m_SamleSize];
		memset(signCiFang_Q, 0, sizeof(float) * m_SamleSize);
		float* zeroPart = new float[m_SamleSize];
		memset(zeroPart, 0, sizeof(float) * m_SamleSize);

		for (int i = 0; i < m_SamleSize; i++)
		{
			signCiFang_I[i] = DataPLLBuff[i].IData * DataPLLBuff[i].IData;
			signCiFang_Q[i] = DataPLLBuff[i].QData * DataPLLBuff[i].QData;
		}

		float* ciFangModel_I = new float[m_SamleSize];
		memset(ciFangModel_I, 0, sizeof(float) * m_SamleSize);
		float* ciFangModel_Q = new float[m_SamleSize];
		memset(ciFangModel_Q, 0, sizeof(float) * m_SamleSize);
		ippfft(signCiFang_I, zeroPart, m_SamleSize, ciFangModel_I);
		ippfft(zeroPart, signCiFang_Q, m_SamleSize, ciFangModel_Q);


		fft_shift(ciFangModel_I, m_SamleSize);
		fft_shift(ciFangModel_Q, m_SamleSize);
		int ciFangMax1_I = getFFTmax(ciFangModel_I, 0, m_SamleSize);
		int ciFangMax2_I = getFFTmax(ciFangModel_I, 0, m_SamleSize / 2 - 20);
		int ciFangMax1_Q = getFFTmax(ciFangModel_Q, 0, m_SamleSize);
		int ciFangMax2_Q = getFFTmax(ciFangModel_Q, 3700 / 2, m_SamleSize / 2 - 20);
		rbEst1 = (float)(ciFangMax1_I - ciFangMax2_I) * fs / m_SamleSize;//ǿ������ת������ֹ��������ͣ�
		rbEst2 = (float)(ciFangMax1_Q - ciFangMax2_Q) * fs / m_SamleSize;
		flag = 1;

		DELETE_ARR(signCiFang_I);
		DELETE_ARR(signCiFang_Q);
		DELETE_ARR(zeroPart);
		DELETE_ARR(ciFangModel_I);
		DELETE_ARR(ciFangModel_Q);
	}


	if (rbEst1 < rbEst2) {
		float* temp = new float[m_SamleSize];
		for (int i = 0; i < m_SamleSize; i++) {
			temp[i] = DataPLLBuff[i].QData;
			DataPLLBuff[i].QData = DataPLLBuff[i].IData;
			DataPLLBuff[i].IData = temp[i];
		}
		DELETE_ARR(temp);
	}

	int nSymbolSyncSize = 0;//ͻ���ڵ�ǰƬ�еĳ���
	Complex* DataSymbolSyncBuff = new Complex[m_SamleSize];
	memset(DataSymbolSyncBuff, 0x00, sizeof(Complex) * m_SamleSize);

	De_UQPSK.SymbolSync(DataPLLBuff, DataSymbolSyncBuff, m_SamleSize, &nSymbolSyncSize, &SySyncBufferI);


	char* DataResultI = new char[nSymbolSyncSize];
	memset(DataResultI, 0x00, sizeof(char) * nSymbolSyncSize);
	Judgment(DataSymbolSyncBuff, nSymbolSyncSize, DataResultI);

	int nBuffSRCSize;
	Complex* DataSRCBuff = new Complex[(m_SamleSize + nDownSampClock - 1) * SRCParam.nUpSampRate / SRCParam.nDownSampRate];
	memset(DataSRCBuff, 0x00, sizeof(Complex) * ((m_SamleSize + nDownSampClock - 1) * SRCParam.nUpSampRate / SRCParam.nDownSampRate));
	De_UQPSK.SRC(DataPLLBuff, m_SamleSize, DataSRCBuff, &nBuffSRCSize, &nDownSampClock, cBlockFilterBuff);


	for (int i = 0; i < nBuffSRCSize; i++) {
		DataSRCBuff[i].IData = DataSRCBuff[i].QData;
	}
	int nSymbolSyncSizeQ = 0;                                       //Q·����ͬ������
	Complex* DataSymbolSyncBuffQ = new Complex[nBuffSRCSize];
	memset(DataSymbolSyncBuffQ, 0x00, sizeof(Complex) * nBuffSRCSize);

	De_UQPSK.SymbolSync(DataSRCBuff, DataSymbolSyncBuffQ, nBuffSRCSize, &nSymbolSyncSizeQ, &SySyncBufferQ);

	signal_UQPSK_Out->burstData->softDistinguishData = new Complex[nSymbolSyncSizeQ];
	for (int i = 0; i < nSymbolSyncSizeQ; i++) {
		signal_UQPSK_Out->burstData->softDistinguishData[i].IData = DataSymbolSyncBuff[i].IData;
		signal_UQPSK_Out->burstData->softDistinguishData[i].QData = DataSymbolSyncBuffQ[i].IData;
	}
	signal_UQPSK_Out->burstData->softDistinguishDataLen = nSymbolSyncSize;
	signal_UQPSK_Out->burstData->softDistinguishDataLenQ = nSymbolSyncSizeQ;
	char* DataResultQ = new char[nSymbolSyncSizeQ];
	memset(DataResultQ, 0x00, sizeof(char) * nSymbolSyncSizeQ);
	Judgment(DataSymbolSyncBuffQ, nSymbolSyncSizeQ, DataResultQ);

	signal_UQPSK_Out->burstData->nDemodulationByteI = new char[nSymbolSyncSize];
	memcpy(signal_UQPSK_Out->burstData->nDemodulationByteI, DataResultI, (sizeof(char)) * nSymbolSyncSize);
	signal_UQPSK_Out->burstData->nDemodulationByteQ = new char[nSymbolSyncSizeQ];
	memcpy(signal_UQPSK_Out->burstData->nDemodulationByteQ, DataResultQ, (sizeof(char)) * nSymbolSyncSizeQ);


	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(DataResultI);
	DELETE_ARR(DataSRCBuff);
	DELETE_ARR(DataSymbolSyncBuffQ);
	DELETE_ARR(DataResultQ);
	//
	DELETE_ARR(DataSymbolSyncBuff);
}

void Demo_UQPSK::InitBlockCFilter()
{

	// ����CFilter��һ���Ѷ���Ľṹ����࣬����nFilterTaps��fFilterCoef�����Ա
	CBlockFilter.nFilterTaps = 16;
	CBlockFilter.fFilterCoef = new float[CBlockFilter.nFilterTaps];

	// ʹ�������ʼ����ֵ
	float tempCoef[] = { -0.0132465027355890, 8.32286573684456e-18, 0.0362191880550484, -1.37671571265597e-17, -0.0878600226999023, 1.80033486821016e-17, 0.313453567624625, 0.502867539511636, 0.313453567624625, 1.80033486821016e-17, -0.0878600226999023, -1.37671571265597e-17, 0.0362191880550484, 8.32286573684456e-18, -0.0132465027355890, -3.33169527851874e-18 };
	// ����ʼ��ֵ���Ƶ���̬�����������
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
	//m_RbQ = m_Rb / 4;
	for (int i = 0; i < 2; i++) {
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

	nDownSampClock = 1;                                         //������ʱ��ʱ�ӿ���  
	fDownConversionPhase = 0;                                   //�±�Ƶ�е���λֵ
	fPLLNCO = 0;                                                //���໷�еı���NCO
	fPLLPastFreqPart = 0;                                       //���໷�е�Ƶ�ʸ�������
	nPLLBuffSize = 0;
	//******************����ͬ��I·��ʼ��***********
	SySyncBufferI.CSymbolSyncBuff = new Complex[4];
	SySyncBufferI.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBufferI.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBufferI.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBufferI.fSymbolSyncW = 0.5;                           //����ͬ����·�˲�������Ĵ�������ֵ��Ϊ0.5
	SySyncBufferI.fSymbolSyncN = 0.9;                           //����ͬ��NCO�Ĵ�������ֵ��Ϊ1
	SySyncBufferI.fSymbolSyncNTemp = 0.9;                       //����ͬ��NCO��ʱ�ļĴ�������ֵ��Ϊ1
	SySyncBufferI.fSymbolSyncU = 0.6;                           //����ͬ��NCO����Ķ�ʱ��������Ĵ�������ֵ��Ϊ0.6
	SySyncBufferI.nSymbolSyncKK = 0;                            //����ͬ��������ʾTiʱ�����,ָʾu,yI,yQ
	SySyncBufferI.fSymbolSyncPastW = 0;                         //����ͬ��W�Ļ���ֵ
	SySyncBufferI.fSymbolSyncPastN = 0;                         //����ͬ��N�Ļ���ֵ
	SySyncBufferI.fSymbolSyncTimeError1 = 0;                    //����ͬ��time_error����ֵ
	SySyncBufferI.fSymbolSyncTimeError2 = 0;                    //����ͬ��time_error����ֵ

	//******************����ͬ��Q·��ʼ��***********
	SySyncBufferQ.CSymbolSyncBuff = new Complex[4];
	SySyncBufferQ.CSymbolSyncYBuff = new Complex[3];
	memset(SySyncBufferQ.CSymbolSyncBuff, 0x00, sizeof(Complex) * 4);
	memset(SySyncBufferQ.CSymbolSyncYBuff, 0x00, sizeof(Complex) * 3);
	SySyncBufferQ.fSymbolSyncW = 0.5;                           //����ͬ����·�˲�������Ĵ�������ֵ��Ϊ0.5
	SySyncBufferQ.fSymbolSyncN = 0.9;                           //����ͬ��NCO�Ĵ�������ֵ��Ϊ1
	SySyncBufferQ.fSymbolSyncNTemp = 0.9;                       //����ͬ��NCO��ʱ�ļĴ�������ֵ��Ϊ1
	SySyncBufferQ.fSymbolSyncU = 0.6;                           //����ͬ��NCO����Ķ�ʱ��������Ĵ�������ֵ��Ϊ0.6
	SySyncBufferQ.nSymbolSyncKK = 0;                            //����ͬ��������ʾTiʱ�����,ָʾu,yI,yQ
	SySyncBufferQ.fSymbolSyncPastW = 0;                         //����ͬ��W�Ļ���ֵ
	SySyncBufferQ.fSymbolSyncPastN = 0;                         //����ͬ��N�Ļ���ֵ
	SySyncBufferQ.fSymbolSyncTimeError1 = 0;                    //����ͬ��time_error����ֵ
	SySyncBufferQ.fSymbolSyncTimeError2 = 0;                    //����ͬ��time_error����ֵ

}

void Demo_UQPSK::InitSymbolSync()
{
	symbolsyncsactorInit.fSymbolSyncFactor1 = 0.013331990840823 / 10;
	symbolsyncsactorInit.fSymbolSyncFactor2 = 8.889783703667417e-05 / 1000;
	De_UQPSK.m_sAlgDemInit.nSampPerSymb = 4;//ÿ����Ԫ�еĲ�������

	De_UQPSK.symbolsyncsactorInit = symbolsyncsactorInit;
}

void Demo_UQPSK::InitCostasPLL()
{
	float BL;
	float Wn;
	float T_nco;
	float sigma = 0.707;     //��·����ϵ��
	int Ko = 1;              //ѹ����������
	int Kd = 1;              //����������
	float fBLcoef = 0.004;
	int K = Ko * Kd;
	BL = fBLcoef * m_Rb;
	Wn = (8 * sigma * BL / (1 + 4 * sigma * sigma));
	T_nco = 1 / (4 * (float)m_Rb);

	CostasPll.fPLLLoopFilterCoef1 = (2 * sigma * Wn * T_nco) / K;
	CostasPll.fPLLLoopFilterCoef2 = ((Wn * T_nco) * (Wn * T_nco)) / K;
	De_UQPSK.CostasPll = CostasPll;
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

	SRCParam.bUpSampFilterValid = false;

	SRCParam.nUpSampRate = 1;
	SRCParam.nDownSampRate = m_Rb / m_RbQ;
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
		output[n] = sqrt((float)pDst[n].re * (float)pDst[n].re + (float)pDst[n].im * (float)pDst[n].im);   //��������ŵ��Ǹ��źŵ�ģֵ
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
