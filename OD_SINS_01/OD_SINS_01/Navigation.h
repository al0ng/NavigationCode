//Navigation.h

#ifndef NAVIGATION_H
#define NAVIGATION_H


#include "math.h"

#define DDL long double
#define DBL long double
#define UI  unsigned int

#define TF				0.1					//时间更新周期(0.8s)

//宏定义
#define SAMPLE_TIME_S	0.01				//国基子惯导数据采样周期(s)

#define TM_COUNT 		100					//解算频率(Hz)


//物理参数
#define PI			3.14159265358979		//圆周率
#define Rad_Deg		57.29577951308232		//弧度化为度的刻度系数 180./PI
#define Rad_Sec		206264.8062470964		//弧度化为度的刻度系数 180./PI
#define Deg_Rad		0.0174532925199433		//度化为弧度的刻度系数 PI/180.
#define Min_Rad		2.908882086657216e-4	//角分化为弧度
#define Dph_Rps		4.8481368110953599e-6	//度/小时化为rad/s
#define dp05h		PI/180/60.0
//#define Re			6378160.				//地球长半轴, 单位：米
#define Re			6378137.				//地球长半轴, 单位：米
#define Rp			6356752.				//地球短半轴, 单位：米
#define e			0.00335281093			//地球椭圆度 1/298.2572 
#define Wie			7.2921151467e-5     	//地球自转角速率, 单位：弧度/s
#define KM			1000.					//公里
#define Knot_Ms 	0.5144444444			//节化为m/s的系数 1852./3600.0	
//#define g0			9.78049				//标定g
#define g0          9.7803267714

#define mg          9.7803267714/1000.0
//#define mg        	9.78049e-3			//mg


#define MAX_TIME   		30					//连续30秒未校正退出组合



#define STS15 			15					//状态维数
#define WAT  			6					//量测维数

//////////////////////////////////////////
typedef struct
{

//	double     Qn2b[4];			//姿态四元数，qn2b
//	double		Cb2n[3][3];			//姿态矩阵，Cb2n
//	double		Cn2b[3][3];			//姿态矩阵，Cn2b

    double 	Att[3];	
	double		Vg[3];				//速度初值
	double		Pos[3];				//解算位置信息
	
}NAVI_OUT;


/********** 导航解算结构体 **********/
typedef struct
{
	double     Tm;					//解算周期
	UI		T100MS_Flag;		//100MS标志
	double		ThetaA2[3];			//补偿后角增量
	double		DeltaV2[3];			//补偿后速度增量
	double		DV[3];

	NAVI_OUT OutData;			//对外输出数据



	double		Q[4];				//INS/OD组合四元数
	double		Cb_n[3][3];			//INS/OD组合转换矩阵表示Cb2n
	/***** Kalman滤波器使用变量 *****/
	double		Xk[STS15];			//状态估计值
	double		Pk[STS15][STS15];	//状态估计方差阵
	double		Qt[STS15];			//状态噪声方差阵
	double		Ft[STS15][STS15];	//状态转移矩阵

	UI      X2Over_Count;		//滤波溢出计数


	double		DFai[3];			//失准角修正量

	UI    	TransFlag;			//状态转移标志
	UI		TransCount;			//时间更新次数

	UI		FilterCount;		//滤波次数
	UI		MesFlag[1];			//滤波计算标志

	UI		NoFilterTime;	//未滤波时间(0-1: 对准; GPS;)

	//INS/GPS组合
	double		H1[WAT][STS15];	//量测矩阵             六维
	double		R1[WAT][WAT];		//量测噪声方差	
	double		Z1[WAT];			//量测值

}SYSTEM_DATA;
/********** 轨迹微分解算结构体 **********/
typedef struct
{

	double     dqn2b[4];			//姿态四元数微分

	double		dvn[3];				//速度微分
	double		dpos[3];			//位置微分

	double		wbib[3];			//XYZ轴角速率(度/s)
	double		fb[3];				//XYZ轴比力(m/s2)

	double		dvm[3];				//主系统速度微分

	double		dpsi;				//主系统速度航向微分

	double		dfw[3];				//主子系统弹性变形角速度微分
	double		dfa[3];				//主子系统弹性变形角微分
	
}Cacult_DATA;

/********** 地球参数结构体 **********/
typedef struct 				//地球相关参数计算
{
	double 	Rs;				//sinL
	double 	Rc;				//cosL

	double		Rmh;			//子午圈曲率半径
	double		Rnh;			//卯酉圈曲率半径
	double 	g;				//重力加速度

	double 	Wz;				//Wie天向分量
	double		Wn;				//Wie北向分量

	double		Wen_n[3];		//Wen在n系投影
	double		Win_n[3];		//Win在n系投影
	double 	Win_b[3];		//Win在b系投影
	double		A_Cori[3];		//哥氏加速度
}EARTH_DATA;

/********** 初始对准结构体 **********/
typedef struct
{

	double     Qib2b[4];			//载体系相对于载体惯性系姿态四元数
	double		Vib_f[3];			//载体惯性系比力积分
	double		Qin2n[4];			//导航系相对于导航惯性系姿态四元数
	double     Vin_g[3];			//导航惯性系比力积分

	double		mk;					//权重系数分母

	double		Vx[3];				//矢量定姿递推计算中间值
	double		Vy[3];				//矢量定姿递推计算中间值
	double		Vz[3];				//矢量定姿递推计算中间值

	double		Vb[3][3];			//矢量定姿递推计算中间值
	double		Vn[3][3];			//矢量定姿递推计算中间值

	double		Qin2ib[4];			//对准结果姿态四元数
	double		Qn2b[4];			//对准结果姿态四元数

	double		Euler_c[3];			//对准结果欧拉角
	double		Euler_t[3];			//对准结果欧拉角

	double		VPos_0[6];			//对准开始时刻载体速度、位置
	double		cL0;				//初始纬度对应余弦值
	double		sL0;				//初始纬度对应正弦值

	int		Align_Count;		//对准中参考信息调用次数
	double		Ref_T;				//参考位置、速度更新周期

	double		winie[3];			//地球自转角速度投影
	
}ALIGN_DATA;

void Navigation_Init(double MINS[9]);					//初始化函数,对准前调用1次
extern void Navigation(NAVI_OUT *pNaviData, double IMU1[6], double IMU2[6]);		

extern void Attitude_Amend(double Fai[3], double *pCb_n, double Att[3], double Q[4]);
extern double Limit_Angle(double Ang, int Flag);
extern void Refresh_Q(double DTheta_b[3], double Q[4]);
extern void Convert_Q_To_Cbn(double q[4], double *pCbn);
void Convert_Att_To_Q(double Att[3], double q[4]);
void Convert_Att_Cbn(double Att[3], double *pCbn);

void Matrix31_Mult(double *pa, double b[3], double c[3]);
double sdet_33(double *pM);

void Quat_Mul(double Q1[4], double Q2[4], double Q3[4]);

extern void Convert_Cbn_To_Att(double *pCbn, double Att[3]);
extern void Calculate_Earth(double Pos[3], double Vg[3], double *pCb_n, EARTH_DATA *pEarth);
extern void Cross(double A[3], double B[3], double C[3]);
void LCross(double V[3], double *pM);
extern double* VMulf(double *Vin, double f,double *Vout);
extern double* Rv2Quat(double *Rv, double *Qbk2bk_1);
extern void Matrix33_Tran(double *pM);
void Matrix33_Mult(double *pa, double *pb, double *pc);

void QQ2Phi(double Qb2nc[4], double Qb2n[4], double Phi[3]);
void Quat2Rv(double Qb2n[4], double Rv[3]);

#endif
