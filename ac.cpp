#include "ac.h"
#include "threadpool/threadpool.h"
#include <sstream>
threadpool::ThreadPool pool;
//
int ac_calibration::trace(const string &root,const string &rt)
{
	if (root=="-1")
	{
		return 0;
	}

	vector<string> &rfroms = allnode.at(root).from;
	vector<future<int>> vecfut(rfroms.size() - 1);
	for (std::size_t i = 1; i < rfroms.size(); i++)
	{
		// Run a function returning a result
		future<int> f = pool.run(std::bind(&ac_calibration::trace, this, rfroms.at(i),root));
		// ... do some other things in parallel to the thread pool.
		vecfut.at(i - 1)._Swap(f);
		// Wait for the function to finish and get the result
	}
	trace(*rfroms.begin(),root);
	for (vector<future<int>>::iterator it = vecfut.begin(); it != vecfut.end(); it++)
	{
		it->get();
	}
	//入流都算完了 可算河道和流域
	//展示一下
	ac_msg tmsg;
	tmsg.tid = std::this_thread::get_id();
	tmsg.cd = root;
	tmsg.rt = rt;
	resnotify.push(tmsg);
	return 1;
}

void ac_calibration::wait()
{
	pool.wait();
}

void ac_calibration::stringsplit(const string & r, vector<string>& v)
{
	stringstream ss(r);
	string s;
	while (getline(ss,s,','))
	{
		if (!s.empty())
		{
			v.push_back(s);
		}
	}
}

class XinanjiangModel
{
private:
	// FORCING
	double *m_pP;   // 降水数据
	long m_nSteps;  // 模型要运行的步长(一共m_nSteps步)
	long steps;
	// OUTPUT
	double *m_pR;   // 流域内每一步长的产流量(径流深度) 
	double *m_pRs;  // 每一步长的地表径流深(毫米) 
	double *m_pRi;  // 每一步长的壤中流深（毫米）
	double *m_pRg;  // 每一步长的地下径流深(毫米) 
	double *m_pE;   // 每一步长的蒸发(毫米)
	double *m_pQrs;  // 流域出口地表径流量
	double *m_pQri;   // 流域出口壤中流径流流量
	double *m_pQrg;  // 流域出口地下径流量
	double *m_pQ;   // 流域出口的总流量
	double m_U;     // for 24h. U=A(km^2)/3.6/delta_t
					// SOIL
	double *m_pW;     // 流域内土壤湿度
	double *m_pWu;	  // 流域内上层土壤湿度
	double *m_pWl;	  // 流域内下层土壤适度
	double *m_pWd;    // 流域内深层土壤湿度
	double m_Wum;	 // 流域内上层土壤蓄水容量
	double m_Wlm;   // 流域内下层土壤蓄水容量
	double m_Wdm;   // 流域内深层土壤蓄水容量，WDM=WM-WUM-WLM 
	double m_Wu, m_Wd, m_Wl;//流域
					// EVAPORATION
	double *m_pEu;  // 上层土壤蒸发量（毫米）
	double *m_pEl;  // 下层土壤蒸发量（毫米）
	double *m_pEd;  // 深层土壤蒸发量（毫米）
					//runoff
	double *RF;
	// PARAMETER
	double m_Kc;      // 流域蒸散发能力与实测蒸散发值的比
	double m_IM;     // 不透水面积占全流域面积之比
	double m_B;      // 蓄水容量曲线的方次，小流域（几平方公里）B0.1左右
					 // 中等面积（平方公里以内）.2~0.3，较大面积.3~0.4   
	double m_WM;     // 流域平均蓄水容量（毫米）(WM=WUM+WLM+WDM)
	double m_C;      // 流域内深层土壤蒸发系数，江南湿润地区：0.15-0.2，
					 //华北半湿润地区：.09-0.12
	double m_SM;    //自由水蓄水容量
	double m_EX;    //自由水蓄水容量～面积分布曲线指数
	double m_KG;    //地下水日出流系数
	double m_KS;    //地下水日出流系数
	double m_KI;    //壤中流日出流系数
	double m_CG;    //地下水消退系数
	double m_CI;    //壤中流消退系数
	double *m_UH;    // 单元流域上地面径流的单位线
	double m_WMM;     // 流域内最大蓄水容量
	double m_Area;    // 流域面积
	double m_Em;     //蒸发皿数据
	int  m_DeltaT;   // 每一步长的小时数
	int  m_PD;       // 给定数据，用以判断是否时行河道汇流计算
	double m_KKS ,m_KKR,m_KKG;
	double m_S0, m_Fr0;
	double m_Eu, m_El, m_Ed;
	double m_R;
public:
	XinanjiangModel(void);
	~XinanjiangModel(void);
	// 初始化模型
	void InitModel(long nSteps, double Area, int DeltaT, int PD, char *ForcingFile);
	// 设置模型参数
	void SetParameters(double *Params);
	// 运行新安江模型
	void RunModel(void);
	// 保存模拟结果到文件
	void SaveResults(char *FileName);
	// 记录出流数据，用以作图分析
	void Runoff(char *runoff);
private:
	// 进行汇流计算，将径流深度转换为流域出口的流量
	void Routing(void);

	void evaporation();
};


XinanjiangModel::XinanjiangModel(void)
{
	this->m_pP = NULL;
	this->m_pEm = NULL;
	this->m_pE = NULL;
	this->m_pEd = NULL;
	this->m_pEl = NULL;
	this->m_pEu = NULL;
	this->m_pW = NULL;
	this->m_pWd = NULL;
	this->m_pWl = NULL;
	this->m_pWu = NULL;
	this->m_pR = NULL;
	this->m_pRg = NULL;
	this->m_pRi = NULL;
	this->m_pRs = NULL;
	this->m_pQ = NULL;
	this->m_pQrg = NULL;
	this->m_pQri = NULL;
	this->m_pQrs = NULL;
}
XinanjiangModel::~XinanjiangModel(void)
{
	delete[] this->m_pP;
	delete[] this->m_pEm;
	delete[] this->m_pE;
	delete[] this->m_pEd;
	delete[] this->m_pEl;
	delete[] this->m_pEu;
	delete[] this->m_pW;
	delete[] this->m_pWd;
	delete[] this->m_pWl;
	delete[] this->m_pWu;
	delete[] this->m_pR;
	delete[] this->m_pRg;
	delete[] this->m_pRi;
	delete[] this->m_pRs;
	delete[] this->m_pQ;
	delete[] this->m_pQrg;
	delete[] this->m_pQrs;
	delete[] this->m_pQri;
}
// 初始化模型
void XinanjiangModel::InitModel(long nSteps, double Area, int DeltaT, int PD, char  * ForcingFile)
{
	FILE  * fp;
	int i;
	this->m_nSteps = nSteps;
	this->steps = this->m_nSteps + 18;
	// 驱动数据
	this->m_pP = new double[this->steps];
	this->m_pEm = new double[this->steps];
	// 模型输出，蒸散发项
	this->m_pE = new double[this->steps];
	this->m_pEd = new double[this->steps];
	this->m_pEl = new double[this->steps];
	this->m_pEu = new double[this->steps];
	// 模型输出，出流项，经过汇流的产流
	this->m_pQrg = new double[this->steps];
	this->m_pQrs = new double[this->steps];
	this->m_pQri = new double[this->steps];
	this->m_pQ = new double[this->steps];
	// 模型输出，产流项
	this->m_pR = new double[this->steps];
	this->m_pRg = new double[this->steps];
	this->m_pRi = new double[this->steps];
	this->m_pRs = new double[this->steps];
	// 模型状态量，土壤湿度
	this->m_pW = new double[this->steps];
	this->m_pWd = new double[this->steps];
	this->m_pWl = new double[this->steps];
	this->m_pWu = new double[this->steps];
	//runoff值
	this->RF = new double[this->steps];
	for (i = 0; i<this->steps; i++)
	{
		// 驱动数据
		this->m_pP[i] = 0.00;
		this->m_pEm[i] = 0.00;
		// 模型输出，蒸散发项
		this->m_pE[i] = 0.00;
		this->m_pEd[i] = 0.00;
		this->m_pEl[i] = 0.00;
		this->m_pEu[i] = 0.00;
		// 模型输出，出流项，经过汇流的产流
		this->m_pQrg[i] = 0.00;
		this->m_pQrs[i] = 0.00;
		this->m_pQri[i] = 0.00;
		this->m_pQ[i] = 0.00;
		// 模型输出，产流项
		this->m_pR[i] = 0.00;
		this->m_pRg[i] = 0.00;
		this->m_pRi[i] = 0.00;
		this->m_pRs[i] = 0.00;
		// 模型状态量，土壤湿度
		this->m_pW[i] = 0.00;
		this->m_pWd[i] = 0.00;
		this->m_pWl[i] = 0.00;
		this->m_pWu[i] = 0.00;
	}
	this->m_Area = Area;
	this->m_DeltaT = DeltaT;
	this->m_PD = PD;
	this->m_U = this->m_Area / (3.6 * this->m_DeltaT);
	// Forcing文件格式：第一列：降水（单位毫米）空格第二列水面蒸发（毫米）
	if ((fp = fopen(ForcingFile, "r")) == NULL)
	{
		printf("Can not open forcing file!\n"); return;
	}
	for (i = 0; i<this->m_nSteps; i++)
	{
		fscanf(fp, "%lf%lf", &(this->m_pP[i]), &(this->m_pEm[i]));
	}
	fclose(fp);
}
// 设置模型参数
void XinanjiangModel::SetParameters(double*  Params)
{
	this->m_Kc = Params[0];   // (1) 流域蒸散发能力与实测水面蒸发之比
	this->m_IM = Params[1];     // (2) 流域不透水面积占全流域面积之比
	this->m_B = Params[2];     // (3) 蓄水容量曲线的方次
	this->m_Wum = Params[3];     // (4) 上层蓄水容量
	this->m_Wlm = Params[4];     // (5) 下层蓄水容量
	this->m_Wdm = Params[5];     // (6) 深层蓄水容量
	this->m_C = Params[6];     // (7) 深层蒸散发系数
	this->m_SM = Params[7];     // (8)自由水蓄水容量
	this->m_EX = Params[8];   // (9)自由水蓄水容量～面积分布曲线指数
	this->m_KG = Params[9];     // (10)地下水日出流系数
	this->m_KI = Params[10];    // (11)壤中流日出流系数
	this->m_CG = Params[11];    // (12)地下水消退系数
	this->m_CI = Params[12];    // (13)壤中流消退系数
	this->m_WM = this->m_Wum + this->m_Wlm + this->m_Wdm;
	this->m_WMM = this->m_WM * (1.0 + this->m_B) / (1.0 - this->m_IM);
}
// 运行新安江模型
void XinanjiangModel::RunModel(void)
{
	long i;
	// 模型的状态变量
	double PE;  // > 0 时为净雨量;< 0 为蒸发不足量（mm）
	double Ep;  //m_Kc  *  m_pEm[i]
	double P;
	double R;   // 产流深度，包括地表径流、壤中流和地下径流（mm）
	double RB;  // 不透水面上产生的径流深度（mm）
	double RG;  // 地下径流深度（mm）
	double RI;  // 壤中流深度（mm）
	double RS;  // 地表径流深（mm）
	double A;   //土壤湿度为W时土壤含水量折算成的径流深度（mm）
	double E = 0.0;   // 蒸散发(mm)
	double EU = 0.0;   // 上层土壤蒸散发量（mm）
	double EL = 0.0;   // 下层土壤蒸散发量（mm）
	double ED = 0.0;   // 深层土壤蒸散发量（mm）
	double S;
	double FRo;
	double FR;
	double MS;
	double AU;
	double WU = 5.0;   // 流域内上层土壤湿度
	double WL = 55.0;  // 流域内下层土壤适度
	double WD = 40.0;  // 流域内深层土壤湿度
	double W = 100.0;
	double So = 5.0;
	MS = m_SM * (1 + m_EX);
	FRo = 1 - pow((1 - So / MS), m_EX);
	for (i = 0; i<this->m_nSteps; i++)
	{
		// ——————蒸散发计算————————————//
		RB = m_pP[i] * m_IM;     // RB是降在不透水面的降雨量
		P = m_pP[i] * (1 - m_IM);
	/*	Ep = m_Kc * m_pEm[i];
		if ((WU + P) >= Ep)
		{
			EU = Ep; EL = 0; ED = 0;
		}
		else if ((WU + P)<Ep)
		{
			EU = WU + P;
			if (WL >= (m_C * m_Wlm))
			{
				EL = (Ep - EU) * WL / m_Wlm;  ED = 0;
			}
			else if ((m_C * (Ep - EU)) <= WL&&WL<(m_C * m_Wlm))
			{
				EL = m_C * (Ep - EU);  ED = 0;
			}
			else if (WL<m_C * (Ep - EU))
			{
				EL = WL;  ED = m_C * (Ep - EU) - EL;
			}
		}*/
		E = EU + EL + ED;
		PE = P - E;

		W = WU + WL + WD;
		////三水源划分汇流计算
		if (PE>0)
		{
			FR = (R - RB) / PE;
			AU = MS * (1 - pow((1 - So * FRo / FR / m_SM), 1 / (1 + m_EX)));
			if (PE + AU<MS)
				RS = FR * (PE + So * FRo / FR - m_SM + m_SM * pow((1 - (PE
					+ AU) / MS), m_EX + 1));
			else if (PE + AU >= MS)
				RS = FR * (PE + So * Fro / FR - m_SM);
			S = So * Fro / FR + (R – RS) / FR;
			RI = m_KI * S * FR;
			RG = m_KG * S * FR;
			RS += RB;
			R = RS + RI + RG;
			So = S * (1 - m_KI - m_KG);
			FRo = FR;
		}
		else
		{
			S = So;
			FR = 1 - pow((1 – S / MS), m_EX);
			RI = 0.00;
			RG = 0.00;
			So = S * (1 - m_KI - m_KG);
			RS = RB;
			R = RS + RI + RG;
			FRo = FR;
		}
		////三水源划分计算结束
		/* 以下部分是状态量：总蒸发量、上、下和深层土壤的蒸发的保存*/
		/* 1 */	this->m_pE[i] = E;     // 当前步长的蒸发（模型重要输出）
		/* 2 */	this->m_pEu[i] = EU;   // 当前步长上层土壤蒸发
		/* 3 */	this->m_pEl[i] = EL;   // 当前步长下层土壤蒸发
		/* 4 */	this->m_pEd[i] = ED;   // 当前步长深层土壤蒸发
		/* 5 */	this->m_pW[i] = W;	   // 当前步长流域平均土壤含水量
		/* 6 */	this->m_pWu[i] = WU;   // 当前步长流域上层土壤含水量
		/* 7 */	this->m_pWl[i] = WL;   // 当前步长流域下层土壤含水量
		/* 8 */	this->m_pWd[i] = WD;   // 当前步长流域深层土壤含水量
		/* 9 */	this->m_pRg[i] = RG;   // 当前步长流域地下径流深度
		/* 10 */ this->m_pRi[i] = RI;   // 当前步长流域壤中流深度
		/* 11 */	this->m_pRs[i] = RS;   // 当前步长流域地表径流径流深度
		/* 12 */ this->m_pR[i] = R;     // 当前步长的总产流径流深度
	}
	this->Routing();
}
// 进行汇流计算，将径流深度转换为流域出口的流量
void XinanjiangModel::Routing(void)
{
	/////////////    地面径流汇流计算：单位线法       ///////////////////////
	int i, j;
	double B[10000] = { 0.00 };
	if (this->m_PD == 1)
	{
		double UH[] = { 3.71,12.99,38.96,94.63,131.74,154.00,166.99,176.27,178.12,
			172.55,146.58, 90.91,53.80, 31.54,18.55, 9.27, 3.71,0.00 };
		for (i = 0; i<this->m_nSteps; i++)
		{
			for (j = 0; j<18; j++)
			{
				B[i + j] += this->m_pRs[i] * UH[j] / 10.0;
			}
		}
	}
	else
	{
		double UH[] = { 7.18,23.38,63.20,143.10,221.75,365.18,447.40,491.29,
			506.93,504.82,468.46,388.56,309.91,166.49,84.26,40.37,17.56,3.46 };
		for (i = 0; i<this->m_nSteps; i++)
		{
			for (j = 0; j<18; j++)
			{
				B[i + j] += this->m_pRs[i] * UH[j] / 10.0;
			}
		}
	}
	for (i = 0; i<this->steps; i++)
		this->m_pQrs[i] = B[i];
	///// 壤中流汇流计算:线性水库 
	for (i = 1; i<this->steps; i++) {
		this->m_pQri[i] = this->m_CI * this->m_pQri[i - 1]
			+ (1.0 - this->m_CI) * this->m_pRi[i] * this->m_U;
	}
	///// 地下径流汇流计算：线性水库 
	for (i = 1; i<this->steps; i++) {
		this->m_pQrg[i] = this->m_pQrg[i - 1] * this->m_CG
			+ this->m_pRg[i] * (1.0 - this->m_CG) * this->m_U;
	}
	//////单元面积总入流计算
	for (i = 0; i<this->steps; i++)
	{
		this->m_pQ[i] = this->m_pQrs[i] + this->m_pQri[i] + this->m_pQrg[i];
	}
}
//
void XinanjiangModel::evaporation(double &prun)
{
	/* ———————蒸散发计算结束——————————— */
	double Ep = m_Kc * m_Em;//蒸散发能力
	if ((m_Wu + prun) >= Ep)
	{
		m_Eu = Ep; m_El = 0; m_Ed = 0;
	}
	else
	{
		m_Eu = m_Wu + prun;
		if (m_Wl >= (m_C * m_Wlm))
		{
			m_El = (Ep - m_Eu) * m_Wl / m_Wlm;  m_Ed = 0;
		}
		else if ((m_C * (Ep - m_Eu)) <= m_Wl)
		{
			m_El = m_C * (Ep - m_Eu);  m_Ed = 0;
		}
		else
		{
			m_El = m_Wl;  m_Ed = m_C * (Ep - m_Eu) - m_El;
		}
	}
	double PE = prun - m_El - m_Eu - m_Ed;
	//——————子流域产流量计算————————————//
	if (PE <= 0)
	{
		m_R = 0;
	}
	else
	{
		double W = m_Wu + m_Wl + m_Wd;
		double A = m_WMM * (1 - pow((1.0 - W / m_WM), 1.0 / (1 + m_B)));
		// 土壤湿度折算净雨量+降水后蒸发剩余雨量<流域内最大含水容量
		if ((A + PE) < m_WMM)
		{
			// 流域内的产流深度计算
			m_R = PE             /*  降水蒸发后的剩余量*/
				+ W          /* 流域内土壤初始蓄水容量*/
				+ m_WM * pow((1 - (PE + A) / m_WMM), (1 + m_B))
				- m_WM; /* 减去流域平均蓄水容量（m_WM:参数）  */
		}
		// 土壤湿度折算净雨量+降水后蒸发剩余雨量<流域内最大含水容量
		else
		{
			// 流域内的产流深度计算
			m_R = PE             /*  降水蒸发后的剩余量					              +  W   /*  流域内土壤湿度*/
				- m_WM  /*  减去流域平均蓄水容量  */
				+ W;
		}
	}
	//三层蓄水量的计算: WU, WL, WD
	if (prun - m_Eu - m_R > 0.0) {
		if (m_Wu + prun - m_Eu - m_R <= m_Wum) {
			m_Wu = m_Wu + prun - m_Eu - m_R;
			m_Wl = m_Wl - m_El;
			m_Wd = m_Wd - m_Ed;
		}
		else {
			m_Wu = m_Wum;
			if ((m_Wl - m_El) + (m_Wu + prun - m_Eu - m_R - m_Wum) <= m_Wlm) {
				m_Wl = (m_Wl - m_El) + (m_Wu + prun - m_Eu - m_R - m_Wum);
				m_Wd = m_Wd - m_Ed;
			}
			else {
				m_Wl = m_Wlm;
				if ((m_Wd - m_Ed) + (m_Wl - m_El) + (m_Wu + prun - m_Eu - m_R - m_Wum) - m_Wlm <= m_Wdm) {
					m_Wd = (m_Wd - m_Ed) + (m_Wl - m_El) + (m_Wu + prun - m_Eu - m_R - m_Wum) - m_Wlm;
				}
				else {
					m_Wd = m_Wdm;
				}
			}
		}
	}
	else {
		if (m_Wu > fabs(prun - m_Eu - m_R)) {
			m_Wu = m_Wu - fabs(prun - m_Eu - m_R);
			m_Wl = m_Wl-m_El;
			m_Wd = m_Wd - m_Ed;
		}
		else {
			m_Wu = 0;
			m_Wl = m_Wu - fabs(prun - m_Eu - m_R) + m_Wl - m_El;
			m_Wd = m_Wd - m_Ed;
			if (m_Wl < 0) {
				m_Wl = 0;
				m_Wd = m_Wu - fabs(prun - m_Eu - m_R) + (m_Wl - m_El) + (m_Wd - m_Ed);
				if (m_Wd < 0) {
					m_Wd = 0;
				}
			}
		}
	}
}

int _tmain(int argc, _TCHAR *  argv[])
{
	long  nSteps = 942;
	int  DeltaT = 24;
	double Area1 = 1603;
	XinanjiangModel Model1;
	Model1.InitModel(nSteps, Area1, DeltaT, 1, "LFForcingfile.txt");
	//模型参数/*Kc,IM,B,m_Wum,Wlm,Wdm,C,SM,EX,KG,KI,CG,CI 
	double Params1[] = { 0.50,0.01,0.30,10,60,40,0.18,32,1.2,0.075,0.072,0.94,0.7 };
	Model1.SetParameters(Params1);
	Model1.RunModel();
	Model1.SaveResults("来凤站日模型计算结果.txt");
	Model1.Runoff("LF_Q.txt");
	Model1.~XinanjiangModel();
	double  Area2 = 2991;
	XinanjiangModel Model2;
	Model2.InitModel(nSteps, Area2, DeltaT, 0, "YCForcingfile.txt");
	//模型参数/*Kc,IM,B,m_Wum,Wlm,Wdm,C,SM,EX,KG,KI,CG,CI  
	double Params2[] = { 0.75,0.01,0.32,10,60,40,0.18,27,1.2,0.065,0.067,0.96,0.8 };
	Model2.SetParameters(Params2);
	Model2.RunModel();
	Model2.SaveResults("file.txt");
	Model2.Runoff("YC_Q.txt");
	Model2.~XinanjiangModel();
	FILE *fp1, *fp2;
	double Q1[1000], Q2[1000], Q[1000] = { 0.00 };
	ofstream outfile;
	outfile.open("Q.txt");
	if ((fp1 = fopen("LF_Q.txt", "r")) == NULL)
	{
		printf("Can not open the file!\n"); return 0;
	}
	if ((fp2 = fopen("YC_Q.txt", "r")) == NULL)
	{
		printf("Can not open the file!\n"); return 0;
	}
	if (outfile.is_open())
	{
		for (int i = 0; i<960; i++)
		{
			fscanf(fp1, "%lf", &Q1[i]);
			fscanf(fp2, "%lf", &Q2[i]);
			Q[i] = Q1[i] + Q2[i];
			outfile << setprecision(3) << setiosflags(ios::fixed) << Q[i] << endl;
		}
		fclose(fp1);
		fclose(fp2);
	}
	outfile.close();
	return 0;
}
