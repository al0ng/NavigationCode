//Navigation.h

#ifndef NAVIGATION_H
#define NAVIGATION_H


#include "math.h"

#define DDL long double
#define DBL long double
#define UI  unsigned int

#define TF				0.1					//ʱ���������(0.8s)

//�궨��
#define SAMPLE_TIME_S	0.01				//�����ӹߵ����ݲ�������(s)

#define TM_COUNT 		100					//����Ƶ��(Hz)


//�������
#define PI			3.14159265358979		//Բ����
#define Rad_Deg		57.29577951308232		//���Ȼ�Ϊ�ȵĿ̶�ϵ�� 180./PI
#define Rad_Sec		206264.8062470964		//���Ȼ�Ϊ�ȵĿ̶�ϵ�� 180./PI
#define Deg_Rad		0.0174532925199433		//�Ȼ�Ϊ���ȵĿ̶�ϵ�� PI/180.
#define Min_Rad		2.908882086657216e-4	//�Ƿֻ�Ϊ����
#define Dph_Rps		4.8481368110953599e-6	//��/Сʱ��Ϊrad/s
#define dp05h		PI/180/60.0
//#define Re			6378160.				//���򳤰���, ��λ����
#define Re			6378137.				//���򳤰���, ��λ����
#define Rp			6356752.				//����̰���, ��λ����
#define e			0.00335281093			//������Բ�� 1/298.2572 
#define Wie			7.2921151467e-5     	//������ת������, ��λ������/s
#define KM			1000.					//����
#define Knot_Ms 	0.5144444444			//�ڻ�Ϊm/s��ϵ�� 1852./3600.0	
//#define g0			9.78049				//�궨g
#define g0          9.7803267714

#define mg          9.7803267714/1000.0
//#define mg        	9.78049e-3			//mg


#define MAX_TIME   		30					//����30��δУ���˳����



#define STS15 			15					//״̬ά��
#define WAT  			6					//����ά��

//////////////////////////////////////////
typedef struct
{

//	double     Qn2b[4];			//��̬��Ԫ����qn2b
//	double		Cb2n[3][3];			//��̬����Cb2n
//	double		Cn2b[3][3];			//��̬����Cn2b

    double 	Att[3];	
	double		Vg[3];				//�ٶȳ�ֵ
	double		Pos[3];				//����λ����Ϣ
	
}NAVI_OUT;


/********** ��������ṹ�� **********/
typedef struct
{
	double     Tm;					//��������
	UI		T100MS_Flag;		//100MS��־
	double		ThetaA2[3];			//�����������
	double		DeltaV2[3];			//�������ٶ�����
	double		DV[3];

	NAVI_OUT OutData;			//�����������



	double		Q[4];				//INS/OD�����Ԫ��
	double		Cb_n[3][3];			//INS/OD���ת�������ʾCb2n
	/***** Kalman�˲���ʹ�ñ��� *****/
	double		Xk[STS15];			//״̬����ֵ
	double		Pk[STS15][STS15];	//״̬���Ʒ�����
	double		Qt[STS15];			//״̬����������
	double		Ft[STS15][STS15];	//״̬ת�ƾ���

	UI      X2Over_Count;		//�˲��������


	double		DFai[3];			//ʧ׼��������

	UI    	TransFlag;			//״̬ת�Ʊ�־
	UI		TransCount;			//ʱ����´���

	UI		FilterCount;		//�˲�����
	UI		MesFlag[1];			//�˲������־

	UI		NoFilterTime;	//δ�˲�ʱ��(0-1: ��׼; GPS;)

	//INS/GPS���
	double		H1[WAT][STS15];	//�������             ��ά
	double		R1[WAT][WAT];		//������������	
	double		Z1[WAT];			//����ֵ

}SYSTEM_DATA;
/********** �켣΢�ֽ���ṹ�� **********/
typedef struct
{

	double     dqn2b[4];			//��̬��Ԫ��΢��

	double		dvn[3];				//�ٶ�΢��
	double		dpos[3];			//λ��΢��

	double		wbib[3];			//XYZ�������(��/s)
	double		fb[3];				//XYZ�����(m/s2)

	double		dvm[3];				//��ϵͳ�ٶ�΢��

	double		dpsi;				//��ϵͳ�ٶȺ���΢��

	double		dfw[3];				//����ϵͳ���Ա��ν��ٶ�΢��
	double		dfa[3];				//����ϵͳ���Ա��ν�΢��
	
}Cacult_DATA;

/********** ��������ṹ�� **********/
typedef struct 				//������ز�������
{
	double 	Rs;				//sinL
	double 	Rc;				//cosL

	double		Rmh;			//����Ȧ���ʰ뾶
	double		Rnh;			//î��Ȧ���ʰ뾶
	double 	g;				//�������ٶ�

	double 	Wz;				//Wie�������
	double		Wn;				//Wie�������

	double		Wen_n[3];		//Wen��nϵͶӰ
	double		Win_n[3];		//Win��nϵͶӰ
	double 	Win_b[3];		//Win��bϵͶӰ
	double		A_Cori[3];		//���ϼ��ٶ�
}EARTH_DATA;

/********** ��ʼ��׼�ṹ�� **********/
typedef struct
{

	double     Qib2b[4];			//����ϵ������������ϵ��̬��Ԫ��
	double		Vib_f[3];			//�������ϵ��������
	double		Qin2n[4];			//����ϵ����ڵ�������ϵ��̬��Ԫ��
	double     Vin_g[3];			//��������ϵ��������

	double		mk;					//Ȩ��ϵ����ĸ

	double		Vx[3];				//ʸ�����˵��Ƽ����м�ֵ
	double		Vy[3];				//ʸ�����˵��Ƽ����м�ֵ
	double		Vz[3];				//ʸ�����˵��Ƽ����м�ֵ

	double		Vb[3][3];			//ʸ�����˵��Ƽ����м�ֵ
	double		Vn[3][3];			//ʸ�����˵��Ƽ����м�ֵ

	double		Qin2ib[4];			//��׼�����̬��Ԫ��
	double		Qn2b[4];			//��׼�����̬��Ԫ��

	double		Euler_c[3];			//��׼���ŷ����
	double		Euler_t[3];			//��׼���ŷ����

	double		VPos_0[6];			//��׼��ʼʱ�������ٶȡ�λ��
	double		cL0;				//��ʼγ�ȶ�Ӧ����ֵ
	double		sL0;				//��ʼγ�ȶ�Ӧ����ֵ

	int		Align_Count;		//��׼�вο���Ϣ���ô���
	double		Ref_T;				//�ο�λ�á��ٶȸ�������

	double		winie[3];			//������ת���ٶ�ͶӰ
	
}ALIGN_DATA;

void Navigation_Init(double MINS[9]);					//��ʼ������,��׼ǰ����1��
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
