//#include "Algorithm_Demodulation.h"

//float Algorithm_Demodulation::meanArray(float* data, int len) {
//	float sum = 0;
//	int i;
//	for (i = 0; i < len; i++)
//	{
//		sum += data[i];
//	}
//	return sum / len;
//}

//void Algorithm_Demodulation::AGC(Complex* BuffIn, float target_power) {
//	float* data_pow2 = new float[m_sAlgDemInit.llSliceLength];
//	memset(data_pow2, 0, sizeof(float) * m_sAlgDemInit.llSliceLength);
//	for (int i = 0; i < m_sAlgDemInit.llSliceLength; i++)
//	{
//		data_pow2[i] = BuffIn[i].IData * BuffIn[i].IData;
//	}
//	float input_power = meanArray(data_pow2, m_sAlgDemInit.llSliceLength);
//	float gain_factor = sqrt(target_power / input_power);
//	for (int i = 0; i < m_sAlgDemInit.llSliceLength; i++)
//	{
//		BuffIn[i].IData = BuffIn[i].IData * gain_factor;
//		BuffIn[i].QData = BuffIn[i].QData * gain_factor;
//	}
//	DELETE_ARR(data_pow2);
//}

//void Algorithm_Demodulation::down_orth_conversion(float Buffin[], Complex BuffOut[], double& fOrthDownConversionPhase)        //已检测
//{
//	double phase0 = 2 * Pi * m_sAlgDemInit.nOrthDownConversionFc / m_sAlgDemInit.nSigFs;
//	for (int i = 0; i < m_sAlgDemInit.llSliceLength; i++)
//	{
//		fOrthDownConversionPhase += phase0;
//		//新增 2022年2月22日  保护
//		if (fOrthDownConversionPhase > (2 * Pi))
//		{
//			fOrthDownConversionPhase = fOrthDownConversionPhase - 2 * Pi;
//		}
//		BuffOut[i].IData = Buffin[i] * cos(fOrthDownConversionPhase);
//		BuffOut[i].QData = -Buffin[i] * sin(fOrthDownConversionPhase);
//	}
//}

//void Algorithm_Demodulation::Quick_Sort(float BufferIn[], int nQuickSortStart, int nQuickSortEnd)   //已检测
//{
//	int i, j;
//	i = nQuickSortStart;
//	j = nQuickSortEnd;
//	float fQuickSortkey;
//	fQuickSortkey = BufferIn[nQuickSortStart];
//	while (i < j)
//	{
//		while (i < j && fQuickSortkey < BufferIn[j])
//			j--;
//		if (i < j)
//		{
//			BufferIn[i] = BufferIn[j];
//			i++;
//		}
//		while (i < j && BufferIn[i] <= fQuickSortkey)
//			i++;
//		if (i < j)
//		{
//			BufferIn[j] = BufferIn[i];
//			j--;
//		}
//	}
//	BufferIn[i] = fQuickSortkey;
//	if (nQuickSortStart < i)
//		Quick_Sort(BufferIn, nQuickSortStart, j - 1);
//	if (i < nQuickSortEnd)
//		Quick_Sort(BufferIn, j + 1, nQuickSortEnd);
//	return;
//}

//int Algorithm_Demodulation::FindPlace(float BuffIn[], int nBuffSize, float FindNum)    //查找峰值频点与次峰值频点
//{
//	int flag = 0;
//	for (int i = 0; i < nBuffSize; i++)
//	{
//		if (BuffIn[i] == FindNum)
//		{
//			flag = i;
//			break;
//		}
//	}
//	return flag;
//}


//bool Algorithm_Demodulation::Fft(Complex m_SlipFFTBuff[], int FFT_size, float* fftdata)         //已检测
//{
//	//int  p, pow_n, lenNum, endNum, l, j_l, s_l, i_s, j, k, n;
//	//float R1, R2, T1, T2;
//	//Complex* Outputfft;
//	//Outputfft = new Complex[FFT_size];        //32点滑动DFT长度
//	////计算log2（length）
//	//int stages = (log(FFT_size) / log(2));
//	//memset(Outputfft, 0x00, sizeof(Complex) * FFT_size);
//	////进行位置变换
//	//exechangNum(m_SlipFFTBuff, stages, Outputfft, FFT_size);
//	////exechangNum(m_SlipFFTBuff.QData, stages, data_imag_fft, FFT_size);
//	////fft 数据由主函数开辟
//	//for (i_s = 0; i_s <= stages - 1; i_s++)
//	//{
//	//	lenNum = 0;
//	//	l = stages - (i_s + 1);
//	//	j_l = pow(2, l);
//	//	s_l = pow(2, i_s);
//	//	for (j = 0; j <= j_l - 1; j++)
//	//	{
//	//		for (k = 0; k <= s_l - 1; k++)
//	//		{
//	//			p = (int)(lenNum + pow(2, i_s));
//	//			pow_n = stages - i_s - 1;
//	//			R1 = (float)(Outputfft[p].IData * cos(2 * pi * k * pow(2, pow_n) / FFT_size));
//	//			R2 = (float)(Outputfft[p].QData * sin(2 * pi * k * pow(2, pow_n) / FFT_size));
//	//			T1 = (float)(Outputfft[p].IData * sin(2 * pi * k * pow(2, pow_n) / FFT_size));
//	//			T2 = (float)(Outputfft[p].QData * cos(2 * pi * k * pow(2, pow_n) / FFT_size));
//	//			Outputfft[p].IData = Outputfft[lenNum].IData - R1 - R2;
//	//			Outputfft[p].QData = Outputfft[lenNum].QData + T1 - T2;
//	//			Outputfft[lenNum].IData = Outputfft[lenNum].IData + R1 + R2;
//	//			Outputfft[lenNum].QData = Outputfft[lenNum].QData - T1 + T2;
//	//			lenNum = lenNum + 1;
//	//			endNum = (int)(lenNum + pow(2, i_s));
//	//		}
//	//		lenNum = endNum;
//	//	}
//	//}
//	//for (n = 0; n < FFT_size; n++)
//	//{
//	//	fftdata[n] = Outputfft[n].IData * Outputfft[n].IData + Outputfft[n].QData * Outputfft[n].QData;;
//	//}
//	//DELETE_ARR(Outputfft);

//	//Ipp32fc Dst[32];
//	Ipp32fc* pDst = NULL;
//	Ipp32fc* pSrc = NULL;
//	int FFTOrder = (log(FFT_size) / log(2)); //add by zhuxue
//	IppsFFTSpec_C_32fc* pSpec = 0;
//	Ipp8u* pMemSpec = 0;
//	Ipp8u* pMemInit = 0;
//	Ipp8u* pMemBuffer = 0;
//	int sizeSpec = 0;
//	int sizeInit = 0;
//	int sizeBuffer = 0;
//	int flag = IPP_FFT_NODIV_BY_ANY;
//	int sizeFft = (int)FFT_size;
//	//add by zhuxue
//	pSrc = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * sizeFft);
//	pDst = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * sizeFft);
//	memset(pSrc, 0, sizeof(Ipp32fc) * sizeFft);
//	memset(pDst, 0, sizeof(Ipp32fc) * sizeFft);
//	for (int i = 0; i < sizeFft; i++)
//	{
//		pSrc[i].re = m_SlipFFTBuff[i].IData;
//		pSrc[i].im = m_SlipFFTBuff[i].QData;
//	}
//	/// get sizes for required buffers
//	ippsFFTGetSize_C_32fc(FFTOrder, flag, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuffer);
//	//	printf("sizeSpec=%d,sizeInit=%d,sizeBuffer=%d\n", sizeSpec, sizeInit, sizeBuffer);
//		/// allocate memory for required buffers
//	pMemSpec = (Ipp8u*)ippMalloc(sizeSpec);
//	if (sizeInit > 0)
//	{
//		pMemInit = (Ipp8u*)ippMalloc(sizeInit);
//	}
//	if (sizeBuffer > 0)
//	{
//		pMemBuffer = (Ipp8u*)ippMalloc(sizeBuffer);
//	}
//	/// initialize FFT specification structure
//	ippsFFTInit_C_32fc(&pSpec, FFTOrder, flag, ippAlgHintNone, pMemSpec, pMemInit);
//	/// free initialization buffer
//	if (sizeInit > 0)
//	{
//		ippFree(pMemInit);
//	}
//	/// perform forward FFT
//	ippsFFTFwd_CToC_32fc(pSrc, pDst, pSpec, pMemBuffer);
//	for (int n = 0; n < sizeFft; n++)
//	{
//		fftdata[n] = (pDst[n].re * pDst[n].re + pDst[n].im * pDst[n].im);   //高精FFT用
//	}
//	/// ...
//	/// free buffers
//	if (sizeBuffer > 0)
//	{
//		ippFree(pMemBuffer);
//	}
//	ippFree(pMemSpec);
//	ippFree(pDst);
//	ippFree(pSrc);

//	return true;
//}

//// 改变索引值
//void Algorithm_Demodulation::exechangNum(Complex* Inputexechang, int pos, Complex* Outputexechang, int FFT_size)   //已检测
//{
//	//输出的二进制 索引值
//	int* new_index = new int[pos];
//	memset(new_index, 0, pos * sizeof(int));
//	//新索引下表
//	int* out_put = new int[FFT_size];
//	memset(out_put, 0, FFT_size * sizeof(int));

//	int i, j, cnt, k;

//	for (i = 0; i < FFT_size; i++)
//	{
//		for (j = pos - 1; j >= 0; j--)
//		{
//			new_index[j] = (i >> (pos - j - 1)) & 1;
//		}
//		cnt = 1;
//		for (j = 0; j < pos; j++)
//		{
//			out_put[i] += cnt * new_index[j];
//			cnt *= 2;
//		}
//	}
//	for (k = 0; k < FFT_size; k++)
//	{
//		Outputexechang[k] = Inputexechang[out_put[k]];
//	}
//	DELETE_ARR(new_index);
//	DELETE_ARR(out_put);
//}

//int Algorithm_Demodulation::CountOne(unsigned char Buffin[], int nBuffSize)
//{
//	int nCouter = 0;
//	for (int i = 0; i < nBuffSize; i++)
//		nCouter += Buffin[i];

//	return nCouter;
//}


////下变频
//void Algorithm_Demodulation::down_conversion(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float fDownConversionFc, float* fDownConversionPhase)
//{
//	float phase0 = 2 * Pi * fDownConversionFc / m_sAlgDemInit.nSigFs;
//	for (int i = 0; i < nBuffSize; i++)
//	{
//		*fDownConversionPhase += phase0;
//		//新增 2022年2月22日
//		if (*fDownConversionPhase > (2 * Pi))
//		{
//			*fDownConversionPhase = *fDownConversionPhase - 2 * Pi;
//		}
//		if (*fDownConversionPhase < -(2 * Pi))
//		{
//			*fDownConversionPhase = *fDownConversionPhase + 2 * Pi;
//		}
//		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fDownConversionPhase), sin(*fDownConversionPhase), BuffOut[i].IData, BuffOut[i].QData);
//	}
//}

///*******************************CarrierSync************************************
//	BuffIn[]:          //输入及其输出
//	BuffOutPLL[]:      //PLL输出
//	BuffOutFLL[];      //FLL输出
//	fCarrierSyncPLL;   //PLL结构体
//	fCarrierSyncFLL;   //FLL结构体
//	nBuffSize;         //长度
//	NSignalType；      //信号类型
//**************************************************************************/
//void Algorithm_Demodulation::CarrierSync(Complex BuffIn[], Complex BuffOutPLL[], Complex BuffOutFLL[], int nBuffSize, int NSignalType, float* fPLLNCO, float* fPLLPastFreqPart, FreqLLBuffer* FLLBuffer)
//{
//	switch (NSignalType)
//	{
//	case 1: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //BPSK
//	case 2: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //QPSK
//	case 3: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //DEQPSK
//	case 4: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //OQPSK
//	case 5: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //DPSK
//	case 6: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //APSK
//	case 10: PLLCostas(BuffIn, BuffOutPLL, nBuffSize, NSignalType, fPLLNCO, fPLLPastFreqPart); break;              //64QAM
//	//case 6: FLL(BuffIn, BuffOutFLL, nBuffSize, NSignalType, FLLBuffer); break;              //FSK
//	case 7: FLL(BuffIn, BuffOutFLL, nBuffSize, NSignalType, FLLBuffer); break;              //CPM
//	case 8: FLL(BuffIn, BuffOutFLL, nBuffSize, NSignalType, FLLBuffer); break;              //SBPSK
//	case 9: FLL(BuffIn, BuffOutFLL, nBuffSize, NSignalType, FLLBuffer); break;              //SOQPSK
//	default:break;
//	}
//}

///*******************************PLLCostas*********************************
//	BuffIn[]:          //输入及其输出
//	BuffOut[]:          //输入及其输出
//	CPLLCostas:       //PLL结构体
//	nBuffSize;         //长度
//	nSignalType;      //信号类型
//**************************************************************************/
//void Algorithm_Demodulation::PLLCostas(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, float* fPLLNCO, float* fPLLPastFreqPart)
//{
//	float discriminator_out = 0;
//	float pll_phase_part;
//	float freq_control;
//	float past_pll_freq_part;
//	float past_nco_phase;
//	int nPLLSignI = 0;
//	int nPLLSignQ = 0;
//	float PLLCoef1 = CostasPll.fPLLLoopFilterCoef1;
//	float PLLCoef2 = CostasPll.fPLLLoopFilterCoef2;
//	if (nSignalType==10)
//	{
//		for (int i = 0; i < nBuffSize; i++)
//		{

//		}
//	}
//	for (int i = 0; i < nBuffSize; i++)
//	{
//		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(*fPLLNCO), sin(*fPLLNCO), BuffOut[i].IData, BuffOut[i].QData);
//		SIGN(BuffOut[i].IData, nPLLSignI);
//		SIGN(BuffOut[i].QData, nPLLSignQ);
//		switch (nSignalType)
//		{
//		case 1: discriminator_out = nPLLSignI* BuffOut[i].QData; break;              //BPSK
//		case 2: discriminator_out = nPLLSignI * BuffOut[i].QData - nPLLSignQ * BuffOut[i].IData; break;             //QPSK
//		case 3: discriminator_out = nPLLSignI * BuffOut[i].QData - nPLLSignQ * BuffOut[i].IData; break;             //DEQPSK
//		case 4: discriminator_out = nPLLSignI * BuffOut[i].QData - nPLLSignQ * BuffOut[i].IData; break;             //OQPSK
//		case 5: discriminator_out = nPLLSignI * BuffOut[i].QData; break;              //DPSK
//		case 6: {
//			float ccc = BuffIn[i].IData;
//			float I_PLL_PE = pow(BuffOut[i].IData, 3) - 3.0 * BuffOut[i].IData * pow(BuffOut[i].QData, 2);
//			float Q_PLL_PE = 3.0 * BuffOut[i].QData * pow(BuffOut[i].IData, 2) - pow(BuffOut[i].QData, 3);
//			SIGN(I_PLL_PE, nPLLSignI);
//			SIGN(Q_PLL_PE, nPLLSignQ);
//			discriminator_out = nPLLSignI * Q_PLL_PE - nPLLSignQ * I_PLL_PE;
//			break;
//		}//APSK
//		case 10: {
//		}
//		default:break;
//		}
//		//pll_phase_part = discriminator_out * PLLCoef1;
//		//*fPLLPastFreqPart = discriminator_out * PLLCoef1 + *fPLLPastFreqPart;
//		//freq_control = pll_phase_part + *fPLLPastFreqPart;
//		//*fPLLNCO = *fPLLNCO + freq_control * 2 * pi;

//		//新增 2022年2月22日
//		pll_phase_part = discriminator_out * PLLCoef1;
//		past_pll_freq_part = *fPLLPastFreqPart;
//		*fPLLPastFreqPart = discriminator_out * PLLCoef2 + past_pll_freq_part;
//		freq_control = pll_phase_part + *fPLLPastFreqPart;
//		past_nco_phase = *fPLLNCO;
//		*fPLLNCO = past_nco_phase + freq_control * 2 * Pi;
//	}
//}

///*******************************FLL************************************
//	BuffIn[]:          //输入及其输出
//	fCarrierSyncFactor[]:  //载波同步环路系数
//	StEndFlag;         //起始点与结束点
//	nSignalType;      //信号类型
//	Rb;               //码速
//**************************************************************************/
//void Algorithm_Demodulation::FLL(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int nSignalType, FreqLLBuffer* FLLBuffer)
//{
//	float dot = 0;
//	float cross = 0;
//	float freqerror = 0;
//	float freqoutdlf = 0;
//	float NCO_Phase = 0;
//	int NC1 = m_sAlgDemInit.nSampPerSymb / 2;
//	float FLLCoef1 = Fll.fFLLPLLLoopFilterCoef1;
//	float FLLCoef2 = Fll.fFLLPLLLoopFilterCoef2;
//	float FLLDownSampleTs = Fll.fFLLDownSampleTs;
//	for (int i = 0; i < nBuffSize; i++)
//	{
//		float fI_accum = 0;
//		float fQ_accum = 0;
//		COMPLEX_MULTIPLY(BuffIn[i].IData, BuffIn[i].QData, cos(FLLBuffer->fFLLNCOPhase), sin(FLLBuffer->fFLLNCOPhase), FLLBuffer->CFLLbuff[NC1 - 1].IData, FLLBuffer->CFLLbuff[NC1 - 1].QData);
//		for (int j = 0; j < NC1; j++)
//		{
//			fI_accum = fI_accum + FLLBuffer->CFLLbuff[j].IData;
//			fQ_accum = fQ_accum + FLLBuffer->CFLLbuff[j].QData;
//		}
//		memmove(FLLBuffer->CFLLbuff, FLLBuffer->CFLLbuff + 1, (NC1 - 1) * sizeof(Complex));	// 观察窗数组移位
//		dot = FLLBuffer->fFLLIaccum * fI_accum + FLLBuffer->fFLLQaccum * fQ_accum;
//		cross = FLLBuffer->fFLLIaccum * fQ_accum - FLLBuffer->fFLLQaccum * fI_accum;
//		//freqerror = (atan2(cross, dot)) / (2 * pi * FLLDownSampleTs);
//		if (fabs(dot) == 0 && fabs(cross) == 0)
//		{
//			freqerror = 0;
//		}
//		else
//		{
//			freqerror = (atan2(cross, dot)) / (2 * Pi * FLLDownSampleTs);
//		}
//		BuffOut[i].IData = freqerror;
//		BuffOut[i].QData = freqerror;
//		FLLBuffer->fFLLIaccum = fI_accum;
//		FLLBuffer->fFLLQaccum = fQ_accum;
//		freqoutdlf = 2 * FLLBuffer->fFLLfreqoutdlf - FLLBuffer->fFLLPastfreqoutdlf + FLLCoef1 * freqerror - FLLCoef2 * FLLBuffer->fFLLfreqerror;
//		FLLBuffer->fFLLPastfreqoutdlf = FLLBuffer->fFLLfreqoutdlf;
//		FLLBuffer->fFLLfreqoutdlf = freqoutdlf;
//		FLLBuffer->fFLLfreqerror = freqerror;
//		NCO_Phase = FLLBuffer->fFLLNCOPhase + freqoutdlf * 2 * Pi;
//		FLLBuffer->fFLLNCOPhase = NCO_Phase;
//	}
//}

//void Algorithm_Demodulation::SymbolSync(Complex BuffIn[], Complex BuffOut[], int nBuffSize, int* nBuffSizeOut, SymbolSyncBuffer* SySyncBuffer)
//{
//	float fSymbolSyncFI1 = 0, fSymbolSyncFI2 = 0, fSymbolSyncFI3 = 0, fSymbolSyncFI4 = 0;
//	float fSymbolSyncFQ1 = 0, fSymbolSyncFQ2 = 0, fSymbolSyncFQ3 = 0, fSymbolSyncFQ4 = 0;
//	float fSymbolSyncIa = 0;
//	float fSymbolSyncQa = 0;
//	int nSymbolSyncstrobe = 0;
//	float fSymbolSyncFactor1 = symbolsyncsactorInit.fSymbolSyncFactor1;
//	float fSymbolSyncFactor2 = symbolsyncsactorInit.fSymbolSyncFactor2;
//	int nSampPerSymb = m_sAlgDemInit.nSampPerSymb;

//	for (int i = 0; i < nBuffSize; i++)
//	{
//		memmove(SySyncBuffer->CSymbolSyncBuff, SySyncBuffer->CSymbolSyncBuff + 1, (nSampPerSymb - 1) * sizeof(Complex));
//		SySyncBuffer->CSymbolSyncBuff[nSampPerSymb - 1] = BuffIn[i];
//		SySyncBuffer->fSymbolSyncPastN = SySyncBuffer->fSymbolSyncN;
//		SySyncBuffer->fSymbolSyncPastW = SySyncBuffer->fSymbolSyncW;
//		SySyncBuffer->fSymbolSyncNTemp = SySyncBuffer->fSymbolSyncN - SySyncBuffer->fSymbolSyncW;
//		if (SySyncBuffer->fSymbolSyncNTemp > 0)
//		{
//			SySyncBuffer->fSymbolSyncN = SySyncBuffer->fSymbolSyncNTemp;
//		}
//		else
//		{
//			if (SySyncBuffer->fSymbolSyncNTemp < 0)
//			{
//				SySyncBuffer->fSymbolSyncN = (float)SySyncBuffer->fSymbolSyncNTemp - (int)SySyncBuffer->fSymbolSyncNTemp + 1;
//			}
//			else
//			{
//				SySyncBuffer->fSymbolSyncN = SySyncBuffer->fSymbolSyncNTemp;
//			}
//			//插值滤波器一
//			//fSymbolSyncFI1 = 0.5 * SySyncBuffer->CSymbolSyncBuff[3].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[2].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].IData + 0.5 * SySyncBuffer->CSymbolSyncBuff[0].IData;
//			//fSymbolSyncFI2 = (-0.5) * SySyncBuffer->CSymbolSyncBuff[3].IData + 1.5 * SySyncBuffer->CSymbolSyncBuff[2].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[0].IData;
//			//fSymbolSyncFI3 = SySyncBuffer->CSymbolSyncBuff[1].IData;


//			//fSymbolSyncFQ1 = 0.5 * SySyncBuffer->CSymbolSyncBuff[3].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[2].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].QData + 0.5 * SySyncBuffer->CSymbolSyncBuff[0].QData;
//			//fSymbolSyncFQ2 = (-0.5) * SySyncBuffer->CSymbolSyncBuff[3].QData + 1.5 * SySyncBuffer->CSymbolSyncBuff[2].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[0].QData;
//			//fSymbolSyncFQ3 = SySyncBuffer->CSymbolSyncBuff[1].QData;
//			//插值滤波器二
//			fSymbolSyncFI1 = 0 * SySyncBuffer->CSymbolSyncBuff[3].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[2].IData + 0.5 * SySyncBuffer->CSymbolSyncBuff[1].IData - SySyncBuffer->CSymbolSyncBuff[0].IData / 6;
//			fSymbolSyncFI2 = SySyncBuffer->CSymbolSyncBuff[3].IData / 6 + 0.5 * SySyncBuffer->CSymbolSyncBuff[2].IData - 1 * SySyncBuffer->CSymbolSyncBuff[1].IData + 0.5 * SySyncBuffer->CSymbolSyncBuff[0].IData;
//			fSymbolSyncFI3 = -SySyncBuffer->CSymbolSyncBuff[3].IData / 6 + 1 * SySyncBuffer->CSymbolSyncBuff[2].IData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].IData - SySyncBuffer->CSymbolSyncBuff[0].IData / 3;
//			fSymbolSyncFI4 = SySyncBuffer->CSymbolSyncBuff[1].IData;

//			fSymbolSyncFQ1 = 0 * SySyncBuffer->CSymbolSyncBuff[3].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[2].QData + 0.5 * SySyncBuffer->CSymbolSyncBuff[1].QData - SySyncBuffer->CSymbolSyncBuff[0].QData / 6;
//			fSymbolSyncFQ2 = SySyncBuffer->CSymbolSyncBuff[3].QData / 6 + 0.5 * SySyncBuffer->CSymbolSyncBuff[2].QData - 1 * SySyncBuffer->CSymbolSyncBuff[1].QData + 0.5 * SySyncBuffer->CSymbolSyncBuff[0].QData;
//			fSymbolSyncFQ3 = -SySyncBuffer->CSymbolSyncBuff[3].QData / 6 + 1 * SySyncBuffer->CSymbolSyncBuff[2].QData - 0.5 * SySyncBuffer->CSymbolSyncBuff[1].QData - SySyncBuffer->CSymbolSyncBuff[0].QData / 3;
//			fSymbolSyncFQ4 = SySyncBuffer->CSymbolSyncBuff[1].QData;

//			memmove(SySyncBuffer->CSymbolSyncYBuff, SySyncBuffer->CSymbolSyncYBuff + 1, 2 * sizeof(Complex));
//			SySyncBuffer->CSymbolSyncYBuff[2].IData = ((fSymbolSyncFI1 * SySyncBuffer->fSymbolSyncU + fSymbolSyncFI2) * SySyncBuffer->fSymbolSyncU + fSymbolSyncFI3) * SySyncBuffer->fSymbolSyncU + fSymbolSyncFI4;
//			SySyncBuffer->CSymbolSyncYBuff[2].QData = ((fSymbolSyncFQ1 * SySyncBuffer->fSymbolSyncU + fSymbolSyncFQ2) * SySyncBuffer->fSymbolSyncU + fSymbolSyncFQ3) * SySyncBuffer->fSymbolSyncU + fSymbolSyncFQ4;
//			nSymbolSyncstrobe = (SySyncBuffer->nSymbolSyncKK + 1) % 2;
//			if (nSymbolSyncstrobe == 0)
//			{
//				BuffOut[*nBuffSizeOut].IData = SySyncBuffer->CSymbolSyncYBuff[2].IData;
//				BuffOut[*nBuffSizeOut].QData = SySyncBuffer->CSymbolSyncYBuff[2].QData;
//				fSymbolSyncIa = (SySyncBuffer->CSymbolSyncYBuff[2].IData + SySyncBuffer->CSymbolSyncYBuff[0].IData) / 2;
//				fSymbolSyncQa = (SySyncBuffer->CSymbolSyncYBuff[2].QData + SySyncBuffer->CSymbolSyncYBuff[0].QData) / 2;
//				SySyncBuffer->fSymbolSyncTimeError1 = SySyncBuffer->fSymbolSyncTimeError2;
//				SySyncBuffer->fSymbolSyncTimeError2 = (SySyncBuffer->CSymbolSyncYBuff[1].IData - fSymbolSyncIa) * (SySyncBuffer->CSymbolSyncYBuff[2].IData - SySyncBuffer->CSymbolSyncYBuff[0].IData) + (SySyncBuffer->CSymbolSyncYBuff[1].QData - fSymbolSyncQa) * (SySyncBuffer->CSymbolSyncYBuff[2].QData - SySyncBuffer->CSymbolSyncYBuff[0].QData);
//				SySyncBuffer->fSymbolSyncW = SySyncBuffer->fSymbolSyncPastW + fSymbolSyncFactor1 * (SySyncBuffer->fSymbolSyncTimeError2 - SySyncBuffer->fSymbolSyncTimeError1) + fSymbolSyncFactor2 * SySyncBuffer->fSymbolSyncTimeError1;
//				(*nBuffSizeOut)++;
//			}
//			SySyncBuffer->nSymbolSyncKK++;
//			SySyncBuffer->fSymbolSyncU = SySyncBuffer->fSymbolSyncPastN / SySyncBuffer->fSymbolSyncPastW;
//		}
//	}
//	//nSymbolSyncTimeSum += nSymbolSyncMS;  //nSymbolSyncMS每片清零
//}

//void Algorithm_Demodulation::Differential_Decoding(BYTE BuffIn[], BYTE Buffout[], int nSamleSize, int& BCache_Value)
//{
//	int BResultOut;
//	for (int i = 0; i < nSamleSize; i++)
//	{
//		BResultOut = BuffIn[i] - BCache_Value;
//		MOD(BResultOut, m_sAlgDemInit.nSignalDifferentialType, Buffout[i]);
//		BCache_Value = BuffIn[i];
//	}
//}

//void Algorithm_Demodulation::normolize(Complex BuffIn[], int nBuffInLen, Complex BuffOut[], int nSignalType, NormolizeBuffer* NormolizeBuffer)
//{
//	int nNormolizeByte = 0;
//	nNormolizeByte += (nBuffInLen >= m_sAlgDemInit.nBytenNormolizeLen) ? m_sAlgDemInit.nBytenNormolizeLen : nBuffInLen;
//	if (NormolizeBuffer->NormolizeFlag == 0)
//	{
//		for (int i = 10; i < nNormolizeByte; i++)
//		{
//			if (fabs(BuffIn[i].IData) > NormolizeBuffer->NormolizeMax)
//			{
//				NormolizeBuffer->NormolizeMax = fabs(BuffIn[i].IData);
//			}

//		}

//		NormolizeBuffer->NormolizeFlag = 1;
//	}
//	else
//	{
//		NormolizeBuffer->NormolizeMax = NormolizeBuffer->NormolizeMax;
//	}

//	switch (nSignalType)
//	{
//	case 6:
//		for (int i = 0; i < nBuffInLen; i++)
//		{
//			BuffOut[i].IData = BuffIn[i].IData / NormolizeBuffer->NormolizeMax;
//			BuffOut[i].QData = BuffIn[i].QData / NormolizeBuffer->NormolizeMax;
//		}
//		break;              //FSK
//	case 7:
//		for (int i = 0; i < nBuffInLen; i++)
//		{
//			BuffOut[i].IData = BuffIn[i].IData / NormolizeBuffer->NormolizeMax;
//			BuffOut[i].QData = BuffIn[i].QData / NormolizeBuffer->NormolizeMax;
//		}
//		break;              //CPM
//	case 8:
//		for (int i = 0; i < nBuffInLen; i++)
//		{
//			BuffOut[i].IData = 0.5 * BuffIn[i].IData / NormolizeBuffer->NormolizeMax;
//			BuffOut[i].QData = 0.5 * BuffIn[i].QData / NormolizeBuffer->NormolizeMax;
//		}
//		break;             //SBPSK
//	case 9:
//		for (int i = 0; i < nBuffInLen; i++)
//		{
//			BuffOut[i].IData = 0.5 * BuffIn[i].IData / NormolizeBuffer->NormolizeMax;
//			BuffOut[i].QData = 0.5 * BuffIn[i].QData / NormolizeBuffer->NormolizeMax;
//		}
//		break;             //SOQPSK
//	default:break;
//	}
//}

//void Algorithm_Demodulation::BlockFilter(Complex BuffIn[], int nBuffLen, Complex BuffOut[], Complex* cFilterBuff)    //已检测
//{
//	float* fFilterCoef = CFilter.fFilterCoef;
//	int	nFilterTaps = CFilter.nFilterTaps;
//	Complex cMulSum{ 0,0 };

//	for (int i = 0; i < nBuffLen; i++)
//	{
//		memmove(cFilterBuff, cFilterBuff + 1, (nFilterTaps - 1) * sizeof(Complex));
//		cFilterBuff[nFilterTaps - 1] = BuffIn[i];
//		if (CFilter.bCoefEvenSym)
//		{
//			for (int j = 0; j < floor(nFilterTaps / 2); j++)
//			{
//				cMulSum.IData += (cFilterBuff[j].IData + cFilterBuff[nFilterTaps - 1 - j].IData) * fFilterCoef[j];
//				cMulSum.QData += (cFilterBuff[j].QData + cFilterBuff[nFilterTaps - 1 - j].QData) * fFilterCoef[j];
//			}
//			if (nFilterTaps % 2 != 0)
//			{
//				cMulSum.IData += cFilterBuff[(nFilterTaps - 1) / 2].IData * fFilterCoef[(nFilterTaps - 1) / 2];
//				cMulSum.QData += cFilterBuff[(nFilterTaps - 1) / 2].QData * fFilterCoef[(nFilterTaps - 1) / 2];
//			}
//		}
//		else
//		{
//			for (int j = 0; j < nFilterTaps; j++)
//			{
//				cMulSum.IData += cFilterBuff[j].IData * CFilter.fFilterCoef[nFilterTaps - 1 - j];
//				cMulSum.QData += cFilterBuff[j].QData * CFilter.fFilterCoef[nFilterTaps - 1 - j];
//			}
//		}
//		BuffOut[i] = cMulSum;
//		SET_COMP_ZERO(cMulSum.IData, cMulSum.QData);
//	}
//	//是否需要清理空间
//	//DELETE_ARR(fFilterCoef);
//}

//void Algorithm_Demodulation::SingleBlockFilter(float BuffIn[], int nBuffLen, float BuffOut[], float* fFilterBuff)
//{
//	float* fFilterCoef = SFilter.fFilterCoef;
//	int	nFilterTaps = SFilter.nFilterTaps;
//	float cMulSum = 0;

//	for (int i = 0; i < nBuffLen; i++)
//	{
//		memmove(fFilterBuff, fFilterBuff + 1, (nFilterTaps - 1) * sizeof(Complex));
//		fFilterBuff[nFilterTaps - 1] = BuffIn[i];

//		if (SFilter.bCoefEvenSym)
//		{
//			for (int j = 0; j < floor(nFilterTaps / 2); j++)
//			{
//				cMulSum += (fFilterBuff[j] + fFilterBuff[nFilterTaps - 1 - j]) * fFilterCoef[j];
//			}
//			if (nFilterTaps % 2 != 0)
//			{
//				cMulSum += fFilterBuff[(nFilterTaps - 1) / 2] * fFilterCoef[(nFilterTaps - 1) / 2];
//			}
//		}
//		else
//		{
//			for (int j = 0; j < nFilterTaps; j++)
//			{
//				cMulSum += fFilterBuff[j] * fFilterCoef[nFilterTaps - 1 - j];
//			}
//		}
//		BuffOut[i] = cMulSum;
//		cMulSum = 0;
//	}
//	DELETE_ARR(fFilterCoef);
//}

//void Algorithm_Demodulation::FLLFilter(Complex BuffIn[], int nBuffLen, Complex BuffOut[], Complex* cFilterBuff)
//{
//	float* fFilterCoef = CFllFilter.fFilterCoef;
//	int	nFilterTaps = CFllFilter.nFilterTaps;
//	Complex cMulSum{ 0,0 };

//	for (int i = 0; i < nBuffLen; i++)
//	{
//		memmove(cFilterBuff, cFilterBuff + 1, (nFilterTaps - 1) * sizeof(Complex));
//		cFilterBuff[nFilterTaps - 1] = BuffIn[i];

//		if (CFllFilter.bCoefEvenSym)
//		{
//			for (int j = 0; j < floor(nFilterTaps / 2); j++)
//			{
//				cMulSum.IData += (cFilterBuff[j].IData + cFilterBuff[nFilterTaps - 1 - j].IData) * fFilterCoef[j];
//				cMulSum.QData += (cFilterBuff[j].QData + cFilterBuff[nFilterTaps - 1 - j].QData) * fFilterCoef[j];
//			}
//			if (nFilterTaps % 2 != 0)
//			{
//				cMulSum.IData += cFilterBuff[(nFilterTaps - 1) / 2].IData * fFilterCoef[(nFilterTaps - 1) / 2];
//				cMulSum.QData += cFilterBuff[(nFilterTaps - 1) / 2].QData * fFilterCoef[(nFilterTaps - 1) / 2];
//			}
//		}
//		else
//		{
//			for (int j = 0; j < nFilterTaps; j++)
//			{
//				cMulSum.IData += cFilterBuff[j].IData * fFilterCoef[nFilterTaps - 1 - j];
//				cMulSum.QData += cFilterBuff[j].QData * fFilterCoef[nFilterTaps - 1 - j];
//			}
//		}
//		BuffOut[i] = cMulSum;
//		SET_COMP_ZERO(cMulSum.IData, cMulSum.QData);
//	}
//	//是否需要清理空间
//	//DELETE_ARR(fFilterCoef);
//}
