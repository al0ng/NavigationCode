//Navigation.h

#ifndef NAVIGATION_H
#define NAVIGATION_H


#include "math.h"

#define DDL long double
#define DBL long double
#define UI  unsigned int

#define TF				0.8					//时间更新周期(0.8s)

//宏定义
#define SAMPLE_TIME_S	0.0025				//国基子惯导数据采样周期(s)

#define TM_COUNT 		200					//解算频率(Hz)


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

//	DBL     Qn2b[4];			//姿态四元数，qn2b
//	DBL		Cb2n[3][3];			//姿态矩阵，Cb2n
//	DBL		Cn2b[3][3];			//姿态矩阵，Cn2b

    DBL 	Att[3];	
	DBL		Vg[3];				//速度初值
	DBL		Pos[3];				//解算位置信息
	
}NAVI_OUT;


/********** 导航解算结构体 **********/
typedef struct
{
	DBL     Tm;					//解算周期
	UI		T100MS_Flag;		//100MS标志
	DBL		ThetaA2[3];			//补偿后角增量
	DBL		DeltaV2[3];			//补偿后速度增量
	DBL		DV[3];

	NAVI_OUT OutData;			//对外输出数据



	DBL		Q[4];				//INS/OD组合四元数
	DBL		Cb_n[3][3];			//INS/OD组合转换矩阵表示Cb2n
	/***** Kalman滤波器使用变量 *****/
	DBL		Xk[STS15];			//状态估计值
	DBL		Pk[STS15][STS15];	//状态估计方差阵
	DBL		Qt[STS15];			//状态噪声方差阵
	DBL		Ft[STS15][STS15];	//状态转移矩阵

	UI      X2Over_Count;		//滤波溢出计数


	DBL		DFai[3];			//失准角修正量

	UI    	TransFlag;			//状态转移标志
	UI		TransCount;			//时间更新次数

	UI		FilterCount;		//滤波次数
	UI		MesFlag[1];			//滤波计算标志

	UI		NoFilterTime;	//未滤波时间(0-1: 对准; GPS;)

	//INS/GPS组合
	DBL		H1[WAT][STS15];	//量测矩阵             六维
	DBL		R1[WAT][WAT];		//量测噪声方差	
	DBL		Z1[WAT];			//量测值

}SYSTEM_DATA;
/********** 轨迹微分解算结构体 **********/
typedef struct
{

	DBL     dqn2b[4];			//姿态四元数微分

	DBL		dvn[3];				//速度微分
	DBL		dpos[3];			//位置微分

	DBL		wbib[3];			//XYZ轴角速率(度/s)
	DBL		fb[3];				//XYZ轴比力(m/s2)

	DBL		dvm[3];				//主系统速度微分

	DBL		dpsi;				//主系统速度航向微分

	DBL		dfw[3];				//主子系统弹性变形角速度微分
	DBL		dfa[3];				//主子系统弹性变形角微分
	
}Cacult_DATA;

/********** 地球参数结构体 **********/
typedef struct 				//地球相关参数计算
{
	DBL 	Rs;				//sinL
	DBL 	Rc;				//cosL

	DBL		Rmh;			//子午圈曲率半径
	DBL		Rnh;			//卯酉圈曲率半径
	DBL 	g;				//重力加速度

	DBL 	Wz;				//Wie天向分量
	DBL		Wn;				//Wie北向分量

	DBL		Wen_n[3];		//Wen在n系投影
	DBL		Win_n[3];		//Win在n系投影
	DBL 	Win_b[3];		//Win在b系投影
	DBL		A_Cori[3];		//哥氏加速度
}EARTH_DATA;

/********** 初始对准结构体 **********/
typedef struct
{

	DBL     Qib2b[4];			//载体系相对于载体惯性系姿态四元数
	DBL		Vib_f[3];			//载体惯性系比力积分
	DBL		Qin2n[4];			//导航系相对于导航惯性系姿态四元数
	DBL     Vin_g[3];			//导航惯性系比力积分

	DBL		mk;					//权重系数分母

	DBL		Vx[3];				//矢量定姿递推计算中间值
	DBL		Vy[3];				//矢量定姿递推计算中间值
	DBL		Vz[3];				//矢量定姿递推计算中间值

	DBL		Vb[3][3];			//矢量定姿递推计算中间值
	DBL		Vn[3][3];			//矢量定姿递推计算中间值


	DBL		Qin2ib[4];			//对准结果姿态四元数
	DBL		Qn2b[4];			//对准结果姿态四元数

	DBL		Euler_c[3];			//对准结果欧拉角
	DBL		Euler_t[3];			//对准结果欧拉角

	DBL		VPos_0[6];			//对准开始时刻载体速度、位置
	DBL		cL0;				//初始纬度对应余弦值
	DBL		sL0;				//初始纬度对应正弦值

	int		Align_Count;		//对准中参考信息调用次数
	DBL		Ref_T;				//参考位置、速度更新周期

	DBL		winie[3];			//地球自转角速度投影
	
}ALIGN_DATA;

void Navigation_Init(DBL MINS[9]);					//初始化函数,对准前调用1次
extern void Navigation(NAVI_OUT *pNaviData, DBL IMU1[6], DBL IMU2[6]);		

extern void Init_KF_Filter();
extern void Kalman_Filter();

void Cal_KF_FNav(DBL Lat, DBL Vg[3], DBL DeltaVb[3], DBL Cbn[3][3], EARTH_DATA earth); 


extern void Matrix31_Mult(DBL *pa, DBL b[3], DBL c[3]);
extern int Matrix_Inverse(DBL *pM);
extern int Matrix_Inversesix(DBL pM[6][6]);   //六阶求逆
extern void Attitude_Amend(DBL Fai[3], DBL *pCb_n, DBL Att[3], DDL Q[4]);
extern DBL Limit_Angle(DBL Ang, int Flag);
extern void Refresh_Q(DBL DTheta_b[3], DDL Q[4]);
extern void Convert_Q_To_Cbn(DDL q[4], DBL *pCbn);
void Convert_Att_To_Q(DBL Att[3], DDL q[4]);
void Convert_Att_Cbn(DBL Att[3], DBL *pCbn);

void Matrix31_Mult(DBL *pa, DBL b[3], DBL c[3]);
DBL sdet_33(DBL *pM);

void Quat_Mul(DBL Q1[4], DBL Q2[4], DBL Q3[4]);

extern void Convert_Cbn_To_Att(DBL *pCbn, DBL Att[3]);
extern void Calculate_Earth(DDL Pos[3], DBL Vg[3], DBL *pCb_n, EARTH_DATA *pEarth);
extern void Cross(DBL A[3], DBL B[3], DBL C[3]);
void LCross(DBL V[3], DBL *pM);
extern DBL* VMulf(DBL *Vin, DBL f,DBL *Vout);
extern DBL* Rv2Quat(DBL *Rv, DBL *Qbk2bk_1);
extern void Matrix33_Tran(DBL *pM);
void Matrix33_Mult(DBL *pa, DBL *pb, DBL *pc);

void QQ2Phi(DBL Qb2nc[4], DBL Qb2n[4], DBL Phi[3]);
void Quat2Rv(DBL Qb2n[4], DBL Rv[3]);

#endif
