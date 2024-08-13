#define _CRT_SECURE_NO_WARNINGS
#include"CDMA.h"

Demo_CDMA::Demo_CDMA(int index) :m_index(index)
{

}

Demo_CDMA::~Demo_CDMA()
{

}

bool Demo_CDMA::DemoduletionInit(const DemodulationInitParamater& demoParameter)
{
	m_demodulationInitParamater = demoParameter;
	subNum = 4096 - ((float)demoParameter.band / demoParameter.fs) * 4096 + 50;
	InitialDemodulation(demoParameter);
	return 1;
}

void Demo_CDMA::demoduletion(Complex* dataInputRealSlice, long long baseTime, DemodulationResult* demodulationResult)
{
	memcpy(((char*)Databuf) + sizeof(Complex) * flag * m_SamleSize, dataInputRealSlice, sizeof(Complex) * m_SamleSize);

	if (++flag < slice)
	{
		demodulationResult->burstData->softDistinguishDataLen = 0;

		return;
	}

	flag = 0;

	Demodulation(Databuf, baseTime, demodulationResult);
}

void Demo_CDMA::Demodulation(Complex* data_input_slice, long long llBaseTime, DemodulationResult* signal_16APSK_Out)
{
	//QFile file;
	//file.setFileName("D:/data/ceshi1.dat");
	//file.open(QIODevice::WriteOnly);
	//file.write((char*)data_input_slice, CDMA_Samplesize * 2 * 4);
	//file.close();

	//AGC
	//AGC(data_input_slice, fAGCPastVc);
	if (flag1)
	{
		float* inputI = new float[CDMA_Samplesize];
		float* inputQ = new float[CDMA_Samplesize];


		//memcpy(inputI, input_I, sizeof(int) * CDMA_Samplesize);
		//memcpy(inputQ, input_Q, sizeof(int) * CDMA_Samplesize);
		for (int i = 0; i < CDMA_Samplesize; i++)
		{
			inputI[i] = data_input_slice[i].IData;
			inputQ[i] = data_input_slice[i].QData;
		}


		// ��Ƶ���й���
		int* pnCode_temp = new int[pnLength];
		memset(pnCode_temp, 0x00, sizeof(int) * pnLength);

		estCDMASpreadCode(inputI, inputQ, fs, pnLength, pnCode_temp);

		memcpy(pnCode, pnCode_temp, sizeof(int) * pnLength);  // ��Ƶ����

		delete[] inputI;
		delete[] inputQ;
		delete[] pnCode_temp;

		flag1 = 0;
	}



	AGCCDMA(data_input_slice, fAGCPastVc, CDMA_Samplesize);


	float* input_I = new float[CDMA_Samplesize];
	memset(input_I, 0x00, sizeof(float) * CDMA_Samplesize);
	float* input_Q = new float[CDMA_Samplesize];
	memset(input_Q, 0x00, sizeof(float) * CDMA_Samplesize);

	for (int i = 0; i < CDMA_Samplesize; i++)
	{
		input_I[i] = data_input_slice[i].IData;
		input_Q[i] = data_input_slice[i].QData;
	}

	CDMADemo(input_I, input_Q, CDMA_Samplesize, signal_16APSK_Out);

	DELETE_ARR(input_I);
	DELETE_ARR(input_Q);
}

void Demo_CDMA::CDMADemo(float* input_I, float* input_Q, int signLen, DemodulationResult* signal_CDMA_Out)
{
	//FILE* fp5;
	//char file_path5[500];

	//sprintf(file_path5, "D:/data/input_I.txt");
	//fp5 = fopen(file_path5, "at");
	//for (int i = 0; i < signLen; i++) {
	//	fprintf(fp5, "%f\n", input_I[i]);
	//}
	//fclose(fp5);

	//sprintf(file_path5, "D:/data/input_Q.txt");
	//fp5 = fopen(file_path5, "at");
	//for (int i = 0; i < signLen; i++) {
	//	fprintf(fp5, "%f\n", input_Q[i]);
	//}
	//fclose(fp5);


	int capture_delay = 0;
	int signal_len = 0;
	signal_len = signLen;

	if (recvState != currentRecvState) {
		currentRecvState = recvState;
		if (!recvState)
		{
			flag_acqusiation = 1;    // ��һ����������-Ҫ���²���
			flag_trcking = 0;		 // ��һ����������-Ҫ���½���
			InitCostasPLL();
		}

	}
	if (flag_acqusiation)
	{
		//flag_acqusiation = 0;
		// ����
		Capture(input_I, input_Q);

		if (delayDot_buf[1] == 0)
		{
			//û���񵽣�������������һ�����ݽ������
			return;
		}

		capture_fc = delayDot_buf[0];
		capture_delay = delayDot_buf[1];
		signal_len = signLen - capture_delay;
		if (capture_delay < 0)
		{
			qDebug() << "333333333333333333333333";
		}
		flag_acqusiation = 0;
	}

	float* captureSign_I = new float[signal_len];
	memset(captureSign_I, 0, sizeof(float) * signal_len);
	float* captureSign_Q = new float[signal_len];
	memset(captureSign_Q, 0, sizeof(float) * signal_len);

	int signal_len_temp = signal_len + pncodelength_sps;  // �������鳤�ȿ��޸ģ�������pncodelength_sps��ֹ���
	float* data_track_I = new float[signal_len_temp];
	memset(data_track_I, 0, sizeof(float) * (signal_len_temp));
	float* data_track_Q = new float[signal_len_temp];
	memset(data_track_Q, 0, sizeof(float) * (signal_len_temp));


	// �±�Ƶ

	//for (int i = 0; i < signal_len; i++)
	//{
	//	captureSign_I[i] = input_I[i + capture_delay] * cos(2 * Pi * capture_fc * (i + 1) / fs) + input_Q[i + capture_delay] * sin(2 * Pi * capture_fc * (i + 1) / fs);
	//	captureSign_Q[i] = input_Q[i + capture_delay] * cos(2 * Pi * capture_fc * (i + 1) / fs) - input_I[i + capture_delay] * sin(2 * Pi * capture_fc * (i + 1) / fs);
	//}

	memcpy((char*)captureSign_I, input_I + capture_delay, sizeof(float) * signal_len);
	memcpy((char*)captureSign_Q, input_Q + capture_delay, sizeof(float) * signal_len);

	//sprintf(file_path1, "D:/data/capture_delay.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < 1; i++) {
	//	fprintf(fp, "%d\n", capture_delay);
	//}
	//fclose(fp);

	// ����
	//int symbol_num = floor(signal_len / pncodelength_sps); //֮ǰ�汾
	//Complex* PNsyndata = new Complex[symbol_num* pncodelength_sps];
	//memset(PNsyndata, 0, sizeof(Complex) * symbol_num* pncodelength_sps);
	//Tracking(captureSign_I, captureSign_Q, symbol_num, PNsyndata);

	// add on 2024/05/06 
	int symbol_num = 0;      // ������Ԫ��
	int track_len = 0;
	int track_len1 = 0;
	if (flag_trcking == 0)
	{
		//��һ�����ݸ���
		int a = floor((signal_len) / pncodelength_sps);
		quotient[flag_trcking] = a; // ��
		remainder[flag_trcking] = signal_len % pncodelength_sps;		 // ����
		symbol_num = quotient[flag_trcking];
		track_len = symbol_num * pncodelength_sps;

		if (remainder[flag_trcking])
		{
			// ������
			memcpy(data_track_temp_I, captureSign_I + track_len, sizeof(float) * remainder[flag_trcking]);
			memcpy(data_track_temp_Q, captureSign_Q + track_len, sizeof(float) * remainder[flag_trcking]);
		}
		else
		{
			// �����������²���
			//flag_acqusiation = 1;
		}

		// ���и���ģ������
		memcpy(data_track_I, captureSign_I, sizeof(float) * track_len);
		memcpy(data_track_Q, captureSign_Q, sizeof(float) * track_len);

	}
	else
	{
		// �ڶ��μ�֮�����ݸ���
		if (remainder[flag_trcking - 1] == 0)
		{
			quotient[flag_trcking] = floor((signal_len) / pncodelength_sps); // ��
			remainder[flag_trcking] = signal_len % pncodelength_sps;		  // ����
			symbol_num = quotient[flag_trcking];
			track_len = symbol_num * pncodelength_sps;

			if (remainder[flag_trcking])
			{
				// ������
				memcpy(data_track_temp_I, captureSign_I + track_len, sizeof(float) * remainder[flag_trcking]);
				memcpy(data_track_temp_Q, captureSign_Q + track_len, sizeof(float) * remainder[flag_trcking]);
			}
			else
			{
				// �����������²���
				//flag_acqusiation = 1;
			}

			// ���и���ģ������
			memcpy(data_track_I, captureSign_I, sizeof(float) * track_len);
			memcpy(data_track_Q, captureSign_Q, sizeof(float) * track_len);
		}
		else
		{
			//������
			int length_increas = pncodelength_sps - remainder[flag_trcking - 1];
			memcpy(data_track_temp_I + remainder[flag_trcking - 1], captureSign_I, sizeof(float) * length_increas); // ����������
			memcpy(data_track_temp_Q + remainder[flag_trcking - 1], captureSign_Q, sizeof(float) * length_increas);

			int length_temp = signal_len - length_increas;

			quotient[flag_trcking] = floor((length_temp) / pncodelength_sps); // ��
			remainder[flag_trcking] = length_temp % pncodelength_sps;		  // ����
			symbol_num = quotient[flag_trcking] + 1;
			track_len = symbol_num * pncodelength_sps;
			track_len1 = quotient[flag_trcking] * pncodelength_sps;

			// ���и���ģ������
			memcpy(data_track_I, data_track_temp_I, sizeof(float) * pncodelength_sps);
			memcpy(data_track_Q, data_track_temp_Q, sizeof(float) * pncodelength_sps);
			memcpy(data_track_I + pncodelength_sps, captureSign_I + length_increas, sizeof(float) * track_len1);
			memcpy(data_track_Q + pncodelength_sps, captureSign_Q + length_increas, sizeof(float) * track_len1);

			// ����������
			if (remainder[flag_trcking])
			{
				// ������
				memcpy(data_track_temp_I, captureSign_I + track_len1 + length_increas, sizeof(float) * remainder[flag_trcking]);
				memcpy(data_track_temp_Q, captureSign_Q + track_len1 + length_increas, sizeof(float) * remainder[flag_trcking]);
			}
			else
			{
				// �����������²���
				//flag_acqusiation = 1;
			}
		}
	}
	flag_trcking++; //��־λ�ۼ� 

	if (flag_trcking >= 1000)
	{
		// ����
		remainder[0] = remainder[999];
		// ��λ
		flag_trcking = 1;
	}

	//sprintf(file_path1, "D:/data/data_track_I.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < track_len; i++) {
	//	fprintf(fp, "%f\n", data_track_I[i]);
	//}
	//fclose(fp);

	//sprintf(file_path2, "D:/data/data_track_Q.txt");
	//fp2 = fopen(file_path2, "at");
	//for (int i = 0; i < track_len; i++) {
	//	fprintf(fp2, "%f\n", data_track_Q[i]);
	//}
	//fclose(fp2);

	Complex* PNsyndata = new Complex[track_len];
	memset(PNsyndata, 0, sizeof(Complex) * track_len);

	//Tracking(data_track_I, data_track_Q, symbol_num, PNsyndata);
	Tracking_new(data_track_I, data_track_Q, symbol_num, PNsyndata);


	//sprintf(file_path5, "D:/data/data_track_I.txt");
	//fp5 = fopen(file_path5, "at");
	//for (int i = 0; i < track_len; i++) {
	//	fprintf(fp5, "%f\n", data_track_I[i]);
	//}
	//fclose(fp5);

	//sprintf(file_path5, "D:/data/data_track_Q.txt");
	//fp5 = fopen(file_path5, "at");
	//for (int i = 0; i < track_len; i++) {
	//	fprintf(fp5, "%f\n", data_track_Q[i]);
	//}
	//fclose(fp5);

	//sprintf(file_path5, "D:/data/Signaldata_track.txt");
	//fp5 = fopen(file_path5, "at");
	//for (int i = 0; i < track_len; i++) {
	//	fprintf(fp5, "%f\n", PNsyndata[i].IData);
	//}
	//fclose(fp5);


	// д�ļ� add on 20240721
	//Complex* output = new Complex[track_len];
	//memset(output, 0x00, sizeof(Complex) * track_len);

	//for (int i = 0; i < track_len; i++)
	//{
	//	output[i].IData = data_track_I[i];
	//	output[i].QData = data_track_Q[i];
	//}

	//QFile file;
	//file.setFileName("D:/data/PNsyndata.dat");
	////file.open(QIODevice::WriteOnlyd);					// ����
	//file.open(QIODevice::WriteOnly| QIODevice::Append);  //׷��
	//file.write((char*)output, track_len * 2 * 4);
	//file.close();

	// ��һ�λ���
	//int integratenum = symbol_num * pncodelength_sps / sum_temp1;

	int integratenum = symbol_num;
	int sum_len_temp = pncodelength_sps;

	Complex* integrateSign = new Complex[integratenum];
	memset(integrateSign, 0, sizeof(Complex) * integratenum);

	for (int i = 0; i < integratenum; i++)
	{
		for (int j = 0; j < sum_len_temp; j++)
		{
			//integrateSign[i].IData += PNsyndata[sum_temp1 * i + j].IData;
			//integrateSign[i].QData += PNsyndata[sum_temp1 * i + j].QData;

			integrateSign[i].IData += PNsyndata[sum_len_temp * i + j].IData;
			integrateSign[i].QData += PNsyndata[sum_len_temp * i + j].QData;
		}
	}

	//��һ�� ��ʽһ
	float* Signal_sort = new float[integratenum];
	memset(Signal_sort, 0x00, sizeof(float) * integratenum);

	float* Signal_sortq = new float[integratenum];
	memset(Signal_sortq, 0x00, sizeof(float)* integratenum);

	for (int i = 0; i < integratenum; i++)
	{
		Signal_sort[i] = fabs(integrateSign[i].IData);  // ȡģֵ
		Signal_sortq[i] = fabs(integrateSign[i].QData);  // ȡģֵ
	}
	
	De_CDMA.Quick_Sort(Signal_sort, 0, integratenum - 1);
	De_CDMA.Quick_Sort(Signal_sortq, 0, integratenum - 1);
	float mean_power = mean_function(Signal_sort, integratenum * 9 / 10, integratenum * 9.5 / 10);
	float mean_powerq = mean_function(Signal_sortq, integratenum * 9 / 10, integratenum * 9.5 / 10);

	for (int i = 0; i < integratenum; i++) {
		integrateSign[i].IData = integrateSign[i].IData / mean_power;
		integrateSign[i].QData = integrateSign[i].QData / mean_powerq;
	}


	// ��һ����ʽ��
	//float* integrateSignI = new float[integratenum];
	//float* integrateSignQ = new float[integratenum];
	//Computer.LhtCopyRealAndimag((Ipp32fc*)integrateSign, integrateSignI, integrateSignQ, integratenum);
	//float reti = 0;
	//float retq = 0;
	//Computer.LhtSumFloatL1(integrateSignI, integratenum, &reti);
	//reti /= integratenum;
	//Computer.LhtMulCFloat(integrateSignI, 1 / reti, integratenum, integrateSignI);
	//Computer.LhtSumFloatL1(integrateSignQ, integratenum, &retq);
	//retq /= integratenum;
	//Computer.LhtMulCFloat(integrateSignQ, 1 / retq, integratenum, integrateSignQ);
	//Computer.LhtCopyRealAndimag2Complex((Ipp32fc*)integrateSign, integrateSignI, integrateSignQ, integratenum);



	//Complex* integrateSignFiler = new Complex[integratenum];
	//memset(integrateSignFiler, 0, sizeof(Complex)* integratenum);
	//De_CDMA.BlockFilter(integrateSign, m_SamleSize, integrateSignFiler, cFilterBuff);
	//AGCCDMA(integrateSignFiler, fAGCPastVc, integratenum);


	// �ز�ͬ��
	Complex* DataPLLBuff = new Complex[integratenum];
	memset(DataPLLBuff, 0x00, sizeof(Complex) * integratenum);
	Complex* DataFLLBuff = new Complex[integratenum];
	memset(DataFLLBuff, 0x00, sizeof(Complex) * integratenum);
	De_CDMA.CarrierSync(integrateSign, DataPLLBuff, DataFLLBuff, integratenum, 1, &fPLLNCO, &fPLLPastFreqPart, &FLLBuffer);

	//FILE* fp;
	//char file_path1[500];
	//sprintf(file_path1, "D:/data/DataPLLBuffI.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < integratenum; i++) {
	//	fprintf(fp, "%f\n", DataPLLBuff[i].IData);
	//}
	//fclose(fp);

	//sprintf(file_path1, "D:/data/DataPLLBuffQ.txt");
	//fp = fopen(file_path1, "at");
	//for (int i = 0; i < integratenum; i++) {
	//	fprintf(fp, "%f\n", DataPLLBuff[i].QData);
	//}
	//fclose(fp);


	// �о�
	float* demoResult = new float[symbol_num];		// ������demoResult
	memset(demoResult, 0, sizeof(float) * symbol_num);
	for (int i = 0; i < symbol_num; i++)
	{
		//for (int j = 0; j < sum_temp2; j++)
		//{
		//	demoResult[i] += DataPLLBuff[j + i * sum_temp2].IData;
		//}

		demoResult[i] = DataPLLBuff[i].IData;
		//demoResult[i] = demoResult[i] / 1;   // modify on 20240623
	}

	//��������ͼ���ƽӿ�
	int star_len = symbol_num;
	Complex* DataPLLBuff_star = new Complex[star_len];
	memset(DataPLLBuff_star, 0x00, sizeof(Complex) * star_len);

	for (int i = 0; i < star_len; i++)
	{
		DataPLLBuff_star[i].IData = demoResult[i] * 100;
		DataPLLBuff_star[i].QData = demoResult[i] * 100;
	}

	signal_CDMA_Out->burstData->softDistinguishData = new Complex[star_len];
	memcpy(signal_CDMA_Out->burstData->softDistinguishData, DataPLLBuff_star, (sizeof(Complex)) * star_len);
	signal_CDMA_Out->burstData->softDistinguishDataLen = star_len;


	// �о�
	char* c_symbol = new char[symbol_num];
	memset(c_symbol, 0x00, sizeof(char) * symbol_num);

	//for (int i = 0; i < symbol_num; i++)
	//{
	//	//c_symbol[i] = demoResult[i] > 0 ? 1 : -1;
	//	c_symbol[i] = demoResult[i] > 0 ? 1 : 0;  // modify on 20240622
	//}

	signal_CDMA_Out->burstData->nDemodulationByteI = new char[symbol_num];
	memcpy(signal_CDMA_Out->burstData->nDemodulationByteI, c_symbol, (sizeof(char)) * symbol_num);

	DELETE_ARR(captureSign_I);
	DELETE_ARR(captureSign_Q);
	DELETE_ARR(integrateSign);
	DELETE_ARR(DataPLLBuff_star);
	DELETE_ARR(DataPLLBuff);
	DELETE_ARR(DataFLLBuff);
	DELETE_ARR(demoResult);
	DELETE_ARR(data_track_I);
	DELETE_ARR(data_track_Q);
	// add on 20240626
	DELETE_ARR(PNsyndata);
	DELETE_ARR(c_symbol);
	DELETE_ARR(Signal_sort);
	DELETE_ARR(Signal_sortq);

}

void Demo_CDMA::max_index_array(float* data, int len, float* output)
{
	double max = data[0];
	int i;
	int temp = 0;
	for (i = 1; i < len; i++)
	{
		if ((data[i]) >= max)
		{
			max = (data[i]);
			temp = i;
		}
	}
	output[0] = max;
	output[1] = temp;
}

void Demo_CDMA::ippfft(float* data_real, float* data_imag, int FFT_LENGTH, int flag1, float* Ioutput, float* Qoutput, float* output)
{
	Ipp32fc* pDst = NULL;
	Ipp32fc* pSrc = NULL;
	int FFT_size = FFT_LENGTH;
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
		Ioutput[n] = (float)pDst[n].re;
		Qoutput[n] = (float)pDst[n].im;
		if (flag1 == -1)
		{
			Qoutput[n] = -Qoutput[n];
		}
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

void Demo_CDMA::ippifft(float* data_real, float* data_imag, int FFT_LENGTH, float* ifft_real, float* ifft_imag, float* Module)
{
	Ipp32fc* Src = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * FFT_LENGTH);
	memset(Src, 0, sizeof(Ipp32fc) * FFT_LENGTH);
	Ipp32fc* Dst = (Ipp32fc*)ippMalloc(sizeof(Ipp32fc) * FFT_LENGTH);
	memset(Dst, 0, sizeof(Ipp32fc) * FFT_LENGTH);

	for (int i = 0; i < FFT_LENGTH; i++) {
		Src[i].re = data_real[i];
		Src[i].im = -1 * data_imag[i];
	}
	int FFTOrder = (log(FFT_LENGTH) / log(2));
	IppsFFTSpec_C_32fc* pSpec = 0;

	Ipp8u* pMemSpec = 0;
	Ipp8u* pMemInit = 0;
	Ipp8u* pMemBuffer = 0;

	int sizeSpec = 0;
	int sizeInit = 0;
	int sizeBuffer = 0;

	int flag = IPP_FFT_NODIV_BY_ANY;

	// get sizes for required buffers
	ippsFFTGetSize_C_32fc(FFTOrder, flag, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuffer);

	// allocate memory for required buffers
	pMemSpec = (Ipp8u*)ippMalloc(sizeSpec);

	if (sizeInit > 0) {
		pMemInit = (Ipp8u*)ippMalloc(sizeInit);
	}

	if (sizeBuffer > 0) {
		pMemBuffer = (Ipp8u*)ippMalloc(sizeBuffer);
	}

	// initialize FFT specification structure
	ippsFFTInit_C_32fc(&pSpec, FFTOrder, flag, ippAlgHintNone, pMemSpec, pMemInit);

	// free initialization buffer
	if (sizeInit > 0) {
		ippFree(pMemInit);
	}

	// perform forward FFT
	ippsFFTFwd_CToC_32fc(Src, Dst, pSpec, pMemBuffer);


	// free buffers
	if (sizeBuffer > 0) {
		ippFree(pMemBuffer);
	}



	for (int n = 0; n < FFT_LENGTH; n++)
	{
		ifft_real[n] = (float)Dst[n].re / FFT_LENGTH;
		ifft_imag[n] = (float)-1 * Dst[n].im / FFT_LENGTH;
		Module[n] = sqrt(ifft_real[n] * ifft_real[n] + ifft_imag[n] * ifft_imag[n]);  //��ֵ����
	}

	ippFree(pMemSpec);
	ippFree(Dst);
	ippFree(Src);

}

void Demo_CDMA::Capture(float* inputI, float* inputQ)
{
	float threshold = 0;
	float flag_acq = 0;

	int* delay_index = new int[step_acqusiation];
	float* value_max = new float[step_acqusiation];
	memset(delay_index, 0x00, sizeof(int) * step_acqusiation);
	memset(value_max, 0x00, sizeof(float) * step_acqusiation);

	int pncodelength_sps1 = pncodelength_sps + chip_sps;	//////////////

	float* PNfftdataI = new float[pncodelength_sps1];
	memset(PNfftdataI, 0, sizeof(float) * pncodelength_sps1);
	float* PNfftdataQ = new float[pncodelength_sps1];
	memset(PNfftdataQ, 0, sizeof(float) * pncodelength_sps1);
	float* PNfftdata = new float[pncodelength_sps1];
	memset(PNfftdata, 0, sizeof(float) * pncodelength_sps1);

	float* CapturedataI = new float[pncodelength_sps1];
	memset(CapturedataI, 0, sizeof(float) * pncodelength_sps1);
	float* CapturedataQ = new float[pncodelength_sps1];
	memset(CapturedataQ, 0, sizeof(float) * pncodelength_sps1);

	float* SignalFFTI = new float[pncodelength_sps1];
	memset(SignalFFTI, 0, sizeof(float) * pncodelength_sps1);
	float* SignalFFTQ = new float[pncodelength_sps1];
	memset(SignalFFTQ, 0, sizeof(float) * pncodelength_sps1);
	float* Signaldata = new float[pncodelength_sps1];
	memset(Signaldata, 0, sizeof(float) * pncodelength_sps1);

	float* DotMultI = new float[pncodelength_sps1];
	memset(DotMultI, 0, sizeof(float) * pncodelength_sps1);
	float* DotMultQ = new float[pncodelength_sps1];
	memset(DotMultQ, 0, sizeof(float) * pncodelength_sps1);

	float* SignalifftI = new float[pncodelength_sps1];
	memset(SignalifftI, 0, sizeof(float) * pncodelength_sps1);
	float* SignalifftQ = new float[pncodelength_sps1];
	memset(SignalifftQ, 0, sizeof(float) * pncodelength_sps1);

	float* Signaldata_temp = new float[pncodelength_sps1];
	memset(Signaldata_temp, 0, sizeof(float) * pncodelength_sps1);


	// PN����Ƭ����
	for (int i = 0; i < pncodelength_sps; i++)
	{
		PNcode_est_sps[i] = pnCode[i / chip_sps];
	}

	//float * PNcode_est_sps_temp = new float[pncodelength_sps];
	//memcpy(PNcode_est_sps_temp, PNcode_est_sps, sizeof(float)*pncodelength_sps);

	float* PNcode_est_sps_temp = new float[pncodelength_sps1];
	memcpy(PNcode_est_sps_temp, PNcode_est_sps, sizeof(float) * pncodelength_sps);
	PNcode_est_sps_temp[pncodelength_sps] = -1;
	PNcode_est_sps_temp[pncodelength_sps + 1] = -1;
	PNcode_est_sps_temp[pncodelength_sps + 2] = -1;
	PNcode_est_sps_temp[pncodelength_sps + 3] = -1;

	// PN��fft  -1 �Ƕ�fft���ȡ����
	ippfft(PNcode_est_sps_temp, PNcode_est_sps_temp, pncodelength_sps1, -1, PNfftdataI, PNfftdataQ, PNfftdata);

	//FILE* fp;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/PNfftdataI.txt");
	//fp = fopen(file_path2, "at");
	//for (int i = 0; i < pncodelength_sps1; i++) {
	//	fprintf(fp, "%f\n", PNfftdataI[i]);
	//}
	//fclose(fp);

	//int step_temp = step_acqusiation;
	int step_temp = 1;

	for (int i = 0; i < step_temp; i++)
	{
		float datamax[2] = { 0 };
		for (int j = 0; j < CDMA_Samplesize; j = (j + pncodelength_sps))
		{


			// �±�Ƶ ...�е�����
			for (int k = 0; k < pncodelength_sps; k++)
			{
				CapturedataI[k] = inputI[j + k] * cos(2 * Pi * fc_acq_buff[i] * (k + 1) / fs) + inputQ[j + k] * sin(2 * Pi * fc_acq_buff[i] * (k + 1) / fs);
				CapturedataQ[k] = inputQ[j + k] * cos(2 * Pi * fc_acq_buff[i] * (k + 1) / fs) - inputI[j + k] * sin(2 * Pi * fc_acq_buff[i] * (k + 1) / fs);
			}

			for (int k = pncodelength_sps; k < pncodelength_sps1; k++)
			{
				CapturedataI[k] = 0;
				CapturedataQ[k] = 0;
			}

			// �ź�����iff
			ippfft(CapturedataI, CapturedataQ, pncodelength_sps1, 1, SignalFFTI, SignalFFTQ, Signaldata);

			// �ź���PN��fft�Ĺ��� ���г˻�
			for (int k = 0; k < pncodelength_sps1; k++)
			{
				DotMultI[k] = SignalFFTI[k] * PNfftdataI[k] - SignalFFTQ[k] * PNfftdataQ[k];
				DotMultQ[k] = SignalFFTI[k] * PNfftdataQ[k] + SignalFFTQ[k] * PNfftdataI[k];
			}

			// ifft
			ippifft(DotMultI, DotMultQ, pncodelength_sps1, SignalifftI, SignalifftQ, Signaldata);

			//FILE* fp;
			//char file_path[500];
			//sprintf(file_path, "D:/data/Signaldata_acqusiation.txt");
			//fp = fopen(file_path, "at");
			//for (int i = 0; i < pncodelength_sps1; i++) {
			//	fprintf(fp, "%f\n", Signaldata[i]);
			//}
			//fclose(fp);

			// ��ֵ����
			if (j == 0)
			{
				threshold = mean_function(Signaldata, 0, pncodelength_sps1) * 5; // ����10,����plot(Signaldata)Ч�������޸�
				//threshold = 100;
			}
			memcpy(Signaldata_temp, Signaldata, sizeof(float) * pncodelength_sps1);
			De_CDMA.Quick_Sort(Signaldata_temp, 0, pncodelength_sps1);
			if (Signaldata_temp[pncodelength_sps1 - 1] > threshold)
			{
				flag_acq = 1;
				max_index_array(Signaldata, pncodelength_sps1, datamax);
				value_max[i] = datamax[0];
				delay_index[i] = datamax[1] + j;
				break;
			}
			if (j > 10000)
			{
				// 10000�������㻹û�в��񵽣�ǿ���˳�����ֹ��ѭ��
				break;
			}
		}

		//if (flag_acq)
		//{
		//	cout << "�� " << i + 1 << " ��" << " �Ѿ�����" << endl;
		//}




	}
	float datamax[2] = { 0 };
	max_index_array(value_max, 3, datamax);
	int index = datamax[1];
	delayDot_buf[0] = fc_acq_buff[index]; //��Ƶ
	delayDot_buf[1] = delay_index[index]; //�ӳٵ�

	DELETE_ARR(PNcode_est_sps_temp);
	DELETE_ARR(PNfftdataI);
	DELETE_ARR(PNfftdataQ);
	DELETE_ARR(PNfftdata);
	DELETE_ARR(CapturedataI);
	DELETE_ARR(CapturedataQ);

	DELETE_ARR(SignalFFTI);
	DELETE_ARR(SignalFFTQ);
	DELETE_ARR(Signaldata);
	DELETE_ARR(DotMultI);
	DELETE_ARR(DotMultQ);
	DELETE_ARR(SignalifftI);
	DELETE_ARR(SignalifftQ);


}

void Demo_CDMA::Tracking(float* inputI, float* inputQ, int symbol_num, Complex* PNsyndata)
{
	int signallen = symbol_num * pncodelength_sps;
	int Symbol = pncodelength_sps;

	float* PNcodelocal = new float[Symbol];
	memset(PNcodelocal, 0, sizeof(float) * Symbol);
	memcpy(PNcodelocal, PNcode_est_sps, sizeof(float) * Symbol); // �������PN�븴��



	float* PNcodeCQ = new float[Symbol];
	memset(PNcodeCQ, 0, sizeof(float) * Symbol);
	float* PNcodeZH = new float[Symbol];
	memset(PNcodeZH, 0, sizeof(float) * Symbol);

	float* DLL_Phase_Part = new float[symbol_num + 1];
	memset(DLL_Phase_Part, 0, sizeof(float) * (symbol_num + 1));
	float* DLL_Freq_Part = new float[symbol_num + 1];
	memset(DLL_Freq_Part, 0, sizeof(float) * (symbol_num + 1));
	float* EPNsynch = new float[symbol_num + 1];
	memset(EPNsynch, 0, sizeof(float) * (symbol_num + 1));


	float cqreal = 0; float cqimag = 0;
	float zhreal = 0; float zhimag = 0;
	float cqsum = 0; float zhsum = 0;
	float temp = 0;

	for (int i = 0; i < symbol_num; i++)
	{
		// ���
		for (int j = 0; j < Symbol; j++)
		{

			PNsyndata[i * Symbol + j].IData = inputI[i * Symbol + j] * PNcodelocal[j];
			PNsyndata[i * Symbol + j].QData = inputQ[i * Symbol + j] * PNcodelocal[j];
		}

		// PN�볬ǰ
		for (int j = 0; j < Symbol - 2; j++)
		{
			PNcodeCQ[j] = PNcodelocal[j + 2];
		}
		PNcodeCQ[Symbol - 2] = PNcodelocal[0];
		PNcodeCQ[Symbol - 1] = PNcodelocal[1];

		//PN���ͺ�
		for (int j = 2; j < Symbol; j++)
		{
			PNcodeZH[j] = PNcodelocal[j - 2];
		}
		PNcodeZH[0] = PNcodelocal[Symbol - 2];
		PNcodeZH[1] = PNcodelocal[Symbol - 1];

		cqreal = 0; cqimag = 0;
		zhreal = 0; zhimag = 0;
		cqsum = 0; zhsum = 0;

		// ��ǰ�ͺ�����
		for (int j = 0; j < Symbol; j++)
		{
			cqreal = cqreal + inputI[i * Symbol + j] * PNcodeCQ[j];
			cqimag = cqimag + inputQ[i * Symbol + j] * PNcodeCQ[j];
			zhreal = zhreal + inputI[i * Symbol + j] * PNcodeZH[j];
			zhimag = zhimag + inputQ[i * Symbol + j] * PNcodeZH[j];
		}
		cqsum = cqreal * cqreal + cqimag * cqimag;
		zhsum = zhreal * zhreal + zhimag * zhimag;

		// ��·�˲�
		int jj = i + 1;
		DLL_Phase_Part[jj] = trackcoef.TrackingLoopCoef2 * (cqsum - zhsum) / (cqsum + zhsum);
		DLL_Freq_Part[jj] = 0.5 * (DLL_Freq_Part[jj - 1] + DLL_Phase_Part[jj]);
		EPNsynch[jj] = trackcoef.TrackingLoopCoef1 * (cqsum - zhsum) / (cqsum + zhsum) + DLL_Freq_Part[jj];

		if (EPNsynch[jj] > trackcoef.ThresholdPositive)
		{
			for (int k = 1; k < Symbol; k++)
			{
				PNcodelocal[k] = PNcodeCQ[k - 1];
			}
			PNcodelocal[0] = PNcodeCQ[Symbol - 1];
		}
		else if (EPNsynch[jj] < trackcoef.ThresholdNegative)
		{
			for (int k = 0; k < Symbol - 1; k++)
			{
				PNcodelocal[k] = PNcodeZH[k + 1];
			}
			PNcodelocal[Symbol - 1] = PNcodeZH[0];
		}
	}


	// ��ӡ����
	//FILE* fp;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/EPNsynch.txt");
	//fp = fopen(file_path2, "at");
	//for (int i = 0; i < symbol_num; i++) {
	//	fprintf(fp, "%f\n", EPNsynch[i]);
	//}
	//fclose(fp);

	// add on 20240624
	memcpy(PNcode_est_sps, PNcodelocal, sizeof(float) * Symbol); // ƫ�ƺ��PN����п���

	DELETE_ARR(PNcodelocal);
	DELETE_ARR(PNcodeCQ);
	DELETE_ARR(PNcodeZH);
	DELETE_ARR(DLL_Phase_Part);
	DELETE_ARR(DLL_Freq_Part);
	DELETE_ARR(EPNsynch);

}

void Demo_CDMA::Initial_Tracking()
{
	float ThresholdPositive, ThresholdNegative;
	float BL_PN, sigma, Ko_PN, Kd_PN, K_PN, T_nco, Wn_PN, C1, C2;
	if (pnLength == 127)
	{
		ThresholdPositive = 0.25;
		ThresholdNegative = -0.25;
	}
	else if (pnLength == 511)
	{
		//ThresholdPositive = 100;
		//ThresholdNegative = -100;
		// 20240720
		ThresholdPositive = 0.25;
		ThresholdNegative = -0.25;

	}
	else if (pnLength == 1023)
	{
		//ThresholdPositive = 110;
		//ThresholdNegative = -110;

		// 20240720
		//ThresholdPositive = 0.25;
		//ThresholdNegative = -0.25;

		ThresholdPositive = 10000;
		ThresholdNegative = -10000;

	}

	int Rb = round(fs / pncodelength_sps);
	BL_PN = Rb * 0.5;
	sigma = 0.707;
	Ko_PN = 0.01;
	Kd_PN = 1;
	K_PN = Ko_PN * Kd_PN;
	T_nco = 1.0 / fs;

	Wn_PN = (8 * sigma * BL_PN) / (1 + 4 * sigma * sigma);
	C1 = (8 * sigma * Wn_PN * T_nco) / (K_PN * (4 + 4 * sigma * Wn_PN * T_nco + (Wn_PN * T_nco) * (Wn_PN * T_nco)));
	C2 = (4 * (T_nco * Wn_PN) * (T_nco * Wn_PN)) / (K_PN * (4 + 4 * sigma * Wn_PN * T_nco + (Wn_PN * T_nco) * (Wn_PN * T_nco)));
	trackcoef.ThresholdPositive = ThresholdPositive;
	trackcoef.ThresholdNegative = ThresholdNegative;
	trackcoef.TrackingLoopCoef1 = C1;
	trackcoef.TrackingLoopCoef2 = C2;
}

void Demo_CDMA::InitCostasPLL()
{
	//int Rc = round(fs / sum_temp1);
	int Rc = round(fs / pncodelength_sps); // һ�λ��ֲ������������Ƶ��������
	float BL_f, Wn_f, C1_f, C2_f;
	//BL_f = Rc * 0.01;
	//BL_f = Rc * 0.003;         // add on 20240720
	BL_f = Rc * 0.001;         // add on 20240721

	//BL_f = Rc * 0.0005;  // modify on 20240623
	if (pnLength == 1023)
	{
		BL_f = Rc * 0.005;
	}
	float Ko_f = 1;
	float Kd_f = 1;
	float K_f = Ko_f * Kd_f;
	float sigma = 0.707;
	float T_nco_f = 1.0 / Rc;
	Wn_f = (8 * sigma * BL_f) / (1 + 4 * sigma * sigma);
	C1_f = (2 * sigma * Wn_f * T_nco_f) / (K_f);
	C2_f = ((T_nco_f * Wn_f) * (T_nco_f * Wn_f)) / (K_f);

	CostasPll.fPLLLoopFilterCoef1 = C1_f;
	CostasPll.fPLLLoopFilterCoef2 = C2_f;
	De_CDMA.CostasPll = CostasPll;

}

float Demo_CDMA::mean_function(float* data, int start, int end)
{
	float sum = 0;
	int i;
	for (i = start; i < end; i++)
	{
		sum += data[i];
	}
	return sum / (end - start);
}

void Demo_CDMA::InitialDemodulation(const DemodulationInitParamater& signalinit)
{
	fPLLNCO = 0;
	fPLLPastFreqPart = 0;

	delayDot_buf = new int[2];
	memset(delayDot_buf, 0x00, sizeof(int) * 2);
	capture_fc = 0;
	flag_acqusiation = 1;

	chip_sps = 4;
	sum_temp1 = chip_sps;
	sum_temp2 = signalinit.pnCodeLength;
	pncodelength_sps = chip_sps * signalinit.pnCodeLength;

	fc_acqusiation = 1000;
	step_acqusiation = 3;
	step_fc_acqusiation = fc_acqusiation / step_acqusiation;
	fc_acq_buff = new int[step_acqusiation];
	//fc_acq_buff[0] = -step_fc_acqusiation;
	//fc_acq_buff[1] = 10;
	//fc_acq_buff[2] = step_fc_acqusiation;

	fc_acq_buff[0] = 0;
	fc_acq_buff[1] = 0;
	fc_acq_buff[2] = 0;


	PNcode_est_sps = new float[pncodelength_sps];
	memset(PNcode_est_sps, 0, sizeof(float) * pncodelength_sps);

	data_track_temp_I = new float[pncodelength_sps];// add on 2024/05/06
	memset(data_track_temp_I, 0x00, sizeof(float) * pncodelength_sps);
	data_track_temp_Q = new float[pncodelength_sps];
	memset(data_track_temp_Q, 0x00, sizeof(float) * pncodelength_sps);

	quotient = new int[1000];					// 1000���޸�
	memset(quotient, 0x00, sizeof(int) * 1000);
	remainder = new int[1000];
	memset(remainder, 0x00, sizeof(int) * 1000);


	m_sAlgDemInit.llSliceLength = 8192;

	// �ⲿ������
	pnLength = signalinit.pnCodeLength;                    // ��Ƶ��������
	fs = signalinit.fs;									   // ������ģ�������-��Ӧ��Ƭ���ʵ��ı�	 		

	pnCode = new int[pnLength];

	//int pn_temp[1023] = { 1,1,-1,1,1,-1,1,1,1,1,1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,-1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,1,1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,-1,1,1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,1,-1,-1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,1 };
	//memcpy(pnCode, pn_temp, (sizeof(int)) * pnLength);

	memset(pnCode, 0x00, sizeof(int) * pnLength);
	//memcpy(pnCode, signalinit.pnCode, (sizeof(int)) * pnLength);
	//pnCode = signalinit.pnCode;

	InitCostasPLL();
	//Initial_Tracking();					//��Ӧ����Trackingһ��ʹ��
	InitBlockCFilter();
	Initial_Tracking_new();			//��Ӧ����Tracking_newһ��ʹ��


	// add on 20240720
	flag = 0;               // �������ݱ�־λ
	slice = 128;             // ��Ƶ�����Ƭ��
	//if (pnLength == 127)
	//{
	//	slice = 60;             // ��Ƶ�����Ƭ��
	//}

	//else if (pnLength == 1023)
	//{
	//	slice = 60;             // ��Ƶ�����Ƭ��
	//}
	m_SamleSize = 8192;     // �����źŽ��һ�ε���Ϊ8192
	CDMA_Samplesize = m_SamleSize * slice; // ��Ƶ���һ���ܵ���
	Databuf = new Complex[CDMA_Samplesize];

}

void Demo_CDMA::AGC(Complex* BuffIn, float target_power)
{
	float* data_pow2 = new float[CDMA_Samplesize];
	memset(data_pow2, 0, sizeof(float) * CDMA_Samplesize);
	for (int i = 0; i < CDMA_Samplesize; i++)
	{
		//data_pow2[i] = BuffIn[i].IData * BuffIn[i].IData                           //֮ǰ
		data_pow2[i] = (BuffIn[i].IData * BuffIn[i].IData + BuffIn[i].QData * BuffIn[i].QData) / 2;    //4.26�޸�
	}
	float input_power = meanArray(data_pow2, CDMA_Samplesize);
	float gain_factor = sqrt(target_power / input_power);
	for (int i = 0; i < CDMA_Samplesize; i++)
	{
		BuffIn[i].IData = BuffIn[i].IData * gain_factor;
		BuffIn[i].QData = BuffIn[i].QData * gain_factor;
	}
	DELETE_ARR(data_pow2);
}

float Demo_CDMA::meanArray(float* data, int len)
{
	float sum = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		sum += data[i];
	}
	return sum / len;
}

// add on 20240622
void Demo_CDMA::Tracking_new(float* inputI, float* inputQ, int symbol_num, Complex* PNsyndata)
{
	int signallen = symbol_num * pncodelength_sps;
	int Symbol = pncodelength_sps;

	float* PNcodelocal = new float[Symbol];
	memset(PNcodelocal, 0, sizeof(float) * Symbol);
	memcpy(PNcodelocal, PNcode_est_sps, sizeof(float) * Symbol); // �������PN�븴��

	float* PNcodelocal_temp = new float[Symbol];
	memset(PNcodelocal_temp, 0, sizeof(float) * Symbol);

	float* PNEarly = new float[Symbol];
	memset(PNEarly, 0, sizeof(float) * Symbol);
	float* PNPrompt = new float[Symbol];
	memset(PNPrompt, 0, sizeof(float) * Symbol);
	float* PNLate = new float[Symbol];
	memset(PNLate, 0, sizeof(float) * Symbol);

	int DLL_temp_len = symbol_num + 1;
	float* DLL_Phase_Part = new float[DLL_temp_len];
	memset(DLL_Phase_Part, 0, sizeof(float) * (DLL_temp_len));
	float* DLL_Freq_Part = new float[DLL_temp_len];
	memset(DLL_Freq_Part, 0, sizeof(float) * (DLL_temp_len));
	float* EPNsynch = new float[DLL_temp_len];
	memset(EPNsynch, 0, sizeof(float) * (DLL_temp_len));

	int jj = 0;
	float I_E = 0; float Q_E = 0;
	float I_P = 0; float Q_P = 0;
	float I_L = 0; float Q_L = 0;
	float E = 0;
	float P = 0;
	float L = 0;
	float loop_error_In = 0;
	for (int i = 0; i < symbol_num; i++)
	{

		// PN�볬ǰearly
		for (int j = 0; j < Symbol - 2; j++)
		{
			PNEarly[j] = PNcodelocal[j + 2];
		}
		PNEarly[Symbol - 2] = PNcodelocal[0];
		PNEarly[Symbol - 1] = PNcodelocal[1];

		// PN�뼴ʱprompt
		for (int j = 0; j < Symbol; j++)
		{
			PNPrompt[j] = PNcodelocal[j];
		}

		//PN���ͺ�late
		for (int j = 2; j < Symbol; j++)
		{
			PNLate[j] = PNcodelocal[j - 2];
		}
		PNLate[0] = PNcodelocal[Symbol - 2];
		PNLate[1] = PNcodelocal[Symbol - 1];


		// ���
		for (int j = 0; j < Symbol; j++)
		{
			PNsyndata[i * Symbol + j].IData = inputI[i * Symbol + j] * PNPrompt[j];
			PNsyndata[i * Symbol + j].QData = inputQ[i * Symbol + j] * PNPrompt[j];
		}

		I_E = 0; Q_E = 0;
		I_L = 0; Q_L = 0;

		// ��ǰ��ʱ�ͺ����
		for (int j = 0; j < Symbol; j++)
		{
			I_E = I_E + inputI[i * Symbol + j] * PNEarly[j];
			Q_E = Q_E + inputQ[i * Symbol + j] * PNEarly[j];

			I_P = I_P + inputI[i * Symbol + j] * PNPrompt[j];
			Q_P = Q_P + inputQ[i * Symbol + j] * PNPrompt[j];

			I_L = I_L + inputI[i * Symbol + j] * PNLate[j];
			Q_L = Q_L + inputQ[i * Symbol + j] * PNLate[j];
		}

		E = 0; P = 0; L = 0;

		E = sqrt(I_E * I_E + Q_E * Q_E);
		P = sqrt(I_P * I_P + Q_P * Q_P);
		L = sqrt(I_L * I_L + Q_L * Q_L);


		// ��·�˲�
		jj = i + 1;
		loop_error_In = 0.5 * (E * E - L * L) / (E * E + L * L);
		DLL_Phase_Part[jj] = trackcoef.TrackingLoopCoef_new1 * loop_error_In;
		DLL_Freq_Part[jj] = DLL_Freq_Part[jj - 1] + trackcoef.TrackingLoopCoef_new2 * loop_error_In;
		EPNsynch[jj] = DLL_Phase_Part[jj] + DLL_Freq_Part[jj];

		// �Ƚ�
		if (L > P && E < P)
		{
			// �ź��ͺ�PN��һ����
			memcpy(PNcodelocal_temp + 1, PNcodelocal, sizeof(float) * (Symbol - 1));
			PNcodelocal_temp[0] = PNcodelocal[Symbol - 1];
			memcpy(PNcodelocal, PNcodelocal_temp, sizeof(float) * Symbol);
		}
		else if (L < P && E > P)
		{
			// �źų�ǰPN��һ����
			memcpy(PNcodelocal_temp, PNcodelocal + 1, sizeof(float) * (Symbol - 1));
			PNcodelocal_temp[Symbol - 1] = PNcodelocal[0];
			memcpy(PNcodelocal, PNcodelocal_temp, sizeof(float) * Symbol);
		}
	}

	//FILE* fp;
	//char file_path2[500];
	//sprintf(file_path2, "D:/data/EPNsynch.txt");
	//fp = fopen(file_path2, "at");
	//for (int i = 0; i < symbol_num; i++) {
	//	fprintf(fp, "%f\n", EPNsynch[i]);
	//}
	//fclose(fp);


	// add on 20240720
	memcpy(PNcode_est_sps, PNcodelocal, sizeof(float) * Symbol); // ƫ�ƺ��PN����п���

	DELETE_ARR(PNcodelocal);
	DELETE_ARR(PNcodelocal_temp);
	DELETE_ARR(PNEarly);
	DELETE_ARR(PNPrompt);
	DELETE_ARR(PNLate);
	DELETE_ARR(DLL_Phase_Part);
	DELETE_ARR(DLL_Freq_Part);
	DELETE_ARR(EPNsynch);
}

void Demo_CDMA::AGCCDMA(Complex* BuffIn, float target_power, int len)
{

	float* data_pow2 = new float[len];
	memset(data_pow2, 0, sizeof(float) * len);
	for (int i = 0; i < len; i++)
	{
		//data_pow2[i] = BuffIn[i].IData * BuffIn[i].IData                           //֮ǰ
		data_pow2[i] = (BuffIn[i].IData * BuffIn[i].IData + BuffIn[i].QData * BuffIn[i].QData) / 2;    //4.26�޸�
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

void Demo_CDMA::InitBlockCFilter()
{
	if (fs <= 23e7)
	{
		// ����CFilter��һ���Ѷ���Ľṹ����࣬����nFilterTaps��fFilterCoef�����Ա

		CFilter.nFilterTaps = 9;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// ʹ�������ʼ����ֵ
		float tempCoef[] =
		{
			-0.042463605937143, 1.462582099439346e-17, 0.212318029685715, 0.500262571713977, 0.636954089057146,
			0.500262571713977, 0.212318029685715, 1.462582099439346e-17, -0.042463605937143
		};

		// ����ʼ��ֵ���Ƶ���̬�����������
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	else
	{

		CFilter.nFilterTaps = 23;

		CFilter.fFilterCoef = new float[CFilter.nFilterTaps];

		// ʹ�������ʼ����ֵ
		float tempCoef[] = { 1.29978452894083e-17,0.00643082701099247,-4.38566444872506e-18,-0.0101055853029882,5.68512058168064e-18,0.0181900535453787,-8.12160083097234e-18,-0.0424434582725503,1.46188814957502e-17,0.212217291362751,0.500025212632458,0.636651874088254,0.500025212632458,0.212217291362751,1.46188814957502e-17,-0.0424434582725503,-8.12160083097234e-18,0.0181900535453787,5.68512058168064e-18,-0.0101055853029882,-4.38566444872506e-18,0.00643082701099247,1.29978452894083e-17 };

		// ����ʼ��ֵ���Ƶ���̬�����������
		for (int i = 0; i < CFilter.nFilterTaps; ++i)
		{
			CFilter.fFilterCoef[i] = tempCoef[i];
		}
	}
	cFilterBuff = new Complex[CFilter.nFilterTaps];
	memset(cFilterBuff, 0x00, sizeof(Complex) * CFilter.nFilterTaps);

	CFilter.bCoefEvenSym = false;
	De_CDMA.CFilter = CFilter;
}

void Demo_CDMA::Initial_Tracking_new()
{
	float BL_PN, sigma, Ko_PN, Kd_PN, K_PN, T_nco, Wn_PN, C1, C2;

	float Rb = round(fs / pncodelength_sps);
	//float Rb = 127000 / 4;
	float coefficient_temp = 0.015;
	BL_PN = float(Rb * coefficient_temp);

	sigma = 0.707;
	Ko_PN = 1;
	Kd_PN = 1;
	K_PN = Ko_PN * Kd_PN;
	Wn_PN = (8 * sigma * BL_PN) / (1 + 4 * sigma * sigma);

	T_nco = (float)1.0 / Rb;

	C1 = (float)(2 * sigma * Wn_PN * T_nco) / (K_PN);
	C2 = (float)(T_nco * Wn_PN) * (T_nco * Wn_PN) / (K_PN);

	trackcoef.TrackingLoopCoef_new1 = C1;
	trackcoef.TrackingLoopCoef_new2 = C2;
}

// add on 20240803
void Demo_CDMA::estCDMASpreadCode(float* input_I, float* input_Q, float fs, int pnEstLen, int* pncode)
{

	/*FILE* fp;
	char file_path1[500];

	sprintf(file_path1, "D:/data/input_I11.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < CDMA_Samplesize; i++) {
		fprintf(fp, "%f\n", input_I[i]);
	}
	fclose(fp);

	sprintf(file_path1, "D:/data/input_Q11.txt");
	fp = fopen(file_path1, "at");
	for (int i = 0; i < CDMA_Samplesize; i++) {
		fprintf(fp, "%f\n", input_Q[i]);
	}
	fclose(fp);*/
	int frontLen = 20000;
	float* input_I1 = new float[CDMA_Samplesize - frontLen];
	memset(input_I1, 0, sizeof(float) * (CDMA_Samplesize - frontLen));
	for (int i = 0; i < CDMA_Samplesize - frontLen; i++)
	{
		input_I1[i] = input_I[frontLen - 1 + i];
	}
	/* ��Ƶ���й��� */
	int slideLen = 0;
	int Threshold1 = 0;
	int Threshold2 = 0;
	if (pnEstLen == 127)
	{
		slideLen = pnEstLen * 4;
		Threshold1 = 10;
		Threshold2 = 10;
	}
	else if (pnEstLen == 511)
	{
		slideLen = pnEstLen * 4;
		Threshold1 = 50;
		Threshold2 = 50;
	}
	else
	{
		slideLen = pnEstLen * 3;
		Threshold1 = 400;
		Threshold2 = 200;
	}
	int dataPart = 100;
	int sampleLen = 4 * pnEstLen;// ��Ƶ������󳤶�
	int samplePoint = 4;// ��Ƶ���������
	//int fluRange = 2 * 8 + 1;//������Χ
	/*int s_D[17] = { -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8 };
	float fluData[17] = { 0 };*/
	int fluRange = 2 * 2 + 1;//������Χ
	int s_D[5] = { -2,-1,0,1,2 };
	float fluData[5] = { 0 };
	float* signAndMove = new float[sampleLen];
	memset(signAndMove, 0, sizeof(float) * sampleLen);
	float* slideMax = new float[dataPart * slideLen];
	memset(slideMax, 0, sizeof(float) * dataPart * slideLen);
	float signAndMoveSum = 0;
	//��ֹ������ǰ��û�����ݣ��ӵڶ���sampleLen��ʼ����
	for (int i = 0; i < dataPart; i++)
	{
		for (int j = 0; j < slideLen; j++)
		{
			for (int k = 0; k < fluRange; k++)
			{
				for (int m = 0; m < sampleLen; m++)
				{
					signAndMove[m] = input_I1[(i + 1) * sampleLen + j + m] * input_I1[(i + 2) * sampleLen + j + s_D[k] + m];
					signAndMoveSum += signAndMove[m];
				}
				fluData[k] = abs(signAndMoveSum / sampleLen);
				signAndMoveSum = 0;
			}
			slideMax[j + i * slideLen] = max_function(fluData, 0, 3);
		}
	}
	DELETE_ARR(signAndMove);
	//������ʼ��
	int* initiPointMax = new int[dataPart];
	memset(initiPointMax, 0, sizeof(int) * dataPart);
	int* initiPointMin = new int[dataPart];
	memset(initiPointMin, 0, sizeof(int) * dataPart);
	int* initiPoint = new int[dataPart]; //��ʼ��ֵ��m����С1
	memset(initiPoint, 0, sizeof(int) * dataPart);
	float* dataCache = new float[slideLen];
	memset(dataCache, 0, sizeof(float) * slideLen);
	int transIndex = 0;
	for (int i = 0; i < dataPart; i++)
	{
		for (int j = 0; j < slideLen; j++)
		{
			dataCache[j] = slideMax[j + i * slideLen];
		}
		initiPointMax[i] = maxIndex(dataCache, 0, slideLen);
		initiPointMin[i] = minIndex(dataCache, 0, slideLen);
		if (initiPointMin[i] >= slideLen - 50 || initiPointMin[i] <= 50)
		{
			initiPointMin[i] = 0;
		}
		if (initiPointMax[i] >= slideLen - 50 || initiPointMax[i] <= 50)
		{
			initiPointMax[i] = 0;
		}
		if (initiPointMin[i] != 0 && initiPointMax[i] == 0)
		{
			initiPoint[transIndex] = initiPointMin[i] + i * sampleLen;
			transIndex += 1;
		}
		else if (initiPointMin[i] == 0 && initiPointMax[i] != 0)
		{
			initiPoint[transIndex] = initiPointMax[i] + i * sampleLen;
			transIndex += 1;
		}
		else if (initiPointMin[i] != 0 && initiPointMax[i] != 0)
		{
			initiPoint[transIndex] = initiPointMax[i] + i * sampleLen;
			transIndex += 1;
		}
	}
	DELETE_ARR(slideMax);
	DELETE_ARR(dataCache);
	DELETE_ARR(initiPointMax);
	DELETE_ARR(initiPointMin);
	//�������п��ܵ���Ƶ����
	int possiNum = 0;
	for (int i = 0; i < dataPart; i++)
	{
		if (initiPoint[i] != 0 && initiPoint[i + 1] == 0)
		{
			possiNum = i;
		}
	}
	int* possiSpreadCode = new int[(possiNum + 1) * pnEstLen];
	memset(possiSpreadCode, 0, sizeof(int) * (possiNum + 1) * pnEstLen);
	float* dataCache1 = new float[pnEstLen];
	memset(dataCache1, 0, sizeof(float) * pnEstLen);
	for (int i = 0; i < possiNum + 1; i++)
	{
		for (int j = 0; j < pnEstLen; j++)
		{
			for (int k = 0; k < samplePoint; k++)
			{
				dataCache1[j] += input_I1[initiPoint[i] + j * samplePoint + k + 1];
			}
			if (dataCache1[j] > 0)
			{
				possiSpreadCode[j + i * pnEstLen] = 1;
			}
			else if (dataCache1[j] == 0)
			{
				possiSpreadCode[j + i * pnEstLen] = 0;
			}
			else
			{
				possiSpreadCode[j + i * pnEstLen] = -1;
			}
			dataCache1[j] = 0.0;
		}
	}
	DELETE_ARR(initiPoint);
	//�ж��������ȷ����Ƶ����
	float* dataCache2 = new float[pnEstLen];
	memset(dataCache2, 0, sizeof(float) * pnEstLen);
	float* codeMultiSum = new float[(possiNum + 1)];
	memset(codeMultiSum, 0, sizeof(float) * (possiNum + 1));
	float codeSum = 0;
	for (int i = 0; i < possiNum + 1; i++)
	{
		for (int j = 0; j < pnEstLen; j++)
		{
			dataCache1[j] = possiSpreadCode[j + i * pnEstLen];
		}
		for (int j = 0; j < possiNum + 1; j++)
		{
			for (int k = 0; k < pnEstLen; k++)
			{
				dataCache2[k] = possiSpreadCode[k + j * pnEstLen];
			}
			for (int k = 0; k < pnEstLen; k++)
			{
				codeSum += dataCache1[k] * dataCache2[k];
			}
			if (abs(codeSum) >= pnEstLen - Threshold1)
			{
				codeSum = 1;
			}
			else
			{
				codeSum = 0;
			}
			codeMultiSum[i] += codeSum;
			codeSum = 0;
		}
	}
	int CodeMaxIndex[3000] = { 0 };
	CodeMaxIndex[0] = maxIndex(codeMultiSum, 0, possiNum + 1);
	transIndex = 1;
	for (int i = CodeMaxIndex[0] + 1; i < possiNum + 1; i++)
	{
		if (codeMultiSum[i] == codeMultiSum[CodeMaxIndex[0]])
		{
			CodeMaxIndex[transIndex] = i;
			transIndex += 1;
		}
	}
	DELETE_ARR(codeMultiSum);


	if (transIndex > 2)
	{
		for (int i = transIndex - 3; i < transIndex; i++)
		{
			for (int j = 0; j < pnEstLen; j++)
			{
				dataCache1[j] = possiSpreadCode[CodeMaxIndex[i] * pnEstLen + j];
			}
			for (int j = 0; j < transIndex; j++)
			{
				for (int k = 0; k < pnEstLen; k++)
				{
					dataCache2[k] = possiSpreadCode[CodeMaxIndex[j] * pnEstLen + k];
				}
				for (int k = 0; k < pnEstLen; k++)
				{
					codeSum += dataCache1[k] * dataCache2[k];
				}
				if (abs(codeSum) >= pnEstLen - Threshold2)
				{
					for (int k = 0; k < pnEstLen; k++)
					{
						pncode[k] = dataCache1[k];
					}
					goto endLoop;
				}
				codeSum = 0;
			}
		}
	}
endLoop:
	DELETE_ARR(dataCache1);
	DELETE_ARR(dataCache2);
	if (transIndex == 2)
	{
		for (int i = 0; i < pnEstLen; i++)
		{
			pncode[i] = possiSpreadCode[CodeMaxIndex[1] * pnEstLen + i];
		}
	}
	if (transIndex == 1)
	{
		for (int i = 0; i < pnEstLen; i++)
		{
			pncode[i] = possiSpreadCode[CodeMaxIndex[0] * pnEstLen + i];
		}
	}
	if (pncode[0] == 0)
	{
		for (int i = 0; i < pnEstLen; i++)
		{
			pncode[i] = possiSpreadCode[CodeMaxIndex[transIndex - 1] * pnEstLen + i];
		}
	}
	DELETE_ARR(possiSpreadCode);
	DELETE_ARR(input_I1);

}

int Demo_CDMA::maxIndex(float* data, int start, int end)
{
	int i, max = start;
	for (i = start; i < end; i++)
	{
		if (data[max] < data[i])
		{
			max = i;
		}
	}
	return max;
}

int Demo_CDMA::minIndex(float* data, int start, int end)
{
	int i, min = start;
	for (i = start; i < end; i++)
	{
		if (data[min] > data[i])
		{
			min = i;
		}
	}
	return min;
}

float Demo_CDMA::max_function(float* data, int start, int end)
{
	float max1 = data[start];
	int i;
	for (i = start + 1; i < end; i++)
	{
		if (max1 < data[i])
		{
			max1 = data[i];
		}
	}
	return max1;
}

//�±�Ƶ
//void Demo_CDMA::down_conversion(Complex BuffIn[], Complex BuffOut[], int nBuffSize, float fDownConversionFc, float* fDownConversionPhase)
//{
//	float phase0 = 2 * Pi * fDownConversionFc / m_sAlgDemInit.nSigFs;
//	for (int i = 0; i < nBuffSize; i++)
//	{
//		*fDownConversionPhase += phase0;
//		//���� 2022��2��22��
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