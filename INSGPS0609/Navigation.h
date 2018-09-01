//Navigation.h

#ifndef NAVIGATION_H
#define NAVIGATION_H


#include "math.h"

#define DDL long double
#define DBL long double
#define UI  unsigned int

#define TF				0.8					//ʱ���������(0.8s)

//�궨��
#define SAMPLE_TIME_S	0.0025				//�����ӹߵ����ݲ�������(s)

#define TM_COUNT 		200					//����Ƶ��(Hz)


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

//	DBL     Qn2b[4];			//��̬��Ԫ����qn2b
//	DBL		Cb2n[3][3];			//��̬����Cb2n
//	DBL		Cn2b[3][3];			//��̬����Cn2b

    DBL 	Att[3];	
	DBL		Vg[3];				//�ٶȳ�ֵ
	DBL		Pos[3];				//����λ����Ϣ
	
}NAVI_OUT;


/********** ��������ṹ�� **********/
typedef struct
{
	DBL     Tm;					//��������
	UI		T100MS_Flag;		//100MS��־
	DBL		ThetaA2[3];			//�����������
	DBL		DeltaV2[3];			//�������ٶ�����
	DBL		DV[3];

	NAVI_OUT OutData;			//�����������



	DBL		Q[4];				//INS/OD�����Ԫ��
	DBL		Cb_n[3][3];			//INS/OD���ת�������ʾCb2n
	/***** Kalman�˲���ʹ�ñ��� *****/
	DBL		Xk[STS15];			//״̬����ֵ
	DBL		Pk[STS15][STS15];	//״̬���Ʒ�����
	DBL		Qt[STS15];			//״̬����������
	DBL		Ft[STS15][STS15];	//״̬ת�ƾ���

	UI      X2Over_Count;		//�˲��������


	DBL		DFai[3];			//ʧ׼��������

	UI    	TransFlag;			//״̬ת�Ʊ�־
	UI		TransCount;			//ʱ����´���

	UI		FilterCount;		//�˲�����
	UI		MesFlag[1];			//�˲������־

	UI		NoFilterTime;	//δ�˲�ʱ��(0-1: ��׼; GPS;)

	//INS/GPS���
	DBL		H1[WAT][STS15];	//�������             ��ά
	DBL		R1[WAT][WAT];		//������������	
	DBL		Z1[WAT];			//����ֵ

}SYSTEM_DATA;
/********** �켣΢�ֽ���ṹ�� **********/
typedef struct
{

	DBL     dqn2b[4];			//��̬��Ԫ��΢��

	DBL		dvn[3];				//�ٶ�΢��
	DBL		dpos[3];			//λ��΢��

	DBL		wbib[3];			//XYZ�������(��/s)
	DBL		fb[3];				//XYZ�����(m/s2)

	DBL		dvm[3];				//��ϵͳ�ٶ�΢��

	DBL		dpsi;				//��ϵͳ�ٶȺ���΢��

	DBL		dfw[3];				//����ϵͳ���Ա��ν��ٶ�΢��
	DBL		dfa[3];				//����ϵͳ���Ա��ν�΢��
	
}Cacult_DATA;

/********** ��������ṹ�� **********/
typedef struct 				//������ز�������
{
	DBL 	Rs;				//sinL
	DBL 	Rc;				//cosL

	DBL		Rmh;			//����Ȧ���ʰ뾶
	DBL		Rnh;			//î��Ȧ���ʰ뾶
	DBL 	g;				//�������ٶ�

	DBL 	Wz;				//Wie�������
	DBL		Wn;				//Wie�������

	DBL		Wen_n[3];		//Wen��nϵͶӰ
	DBL		Win_n[3];		//Win��nϵͶӰ
	DBL 	Win_b[3];		//Win��bϵͶӰ
	DBL		A_Cori[3];		//���ϼ��ٶ�
}EARTH_DATA;

/********** ��ʼ��׼�ṹ�� **********/
typedef struct
{

	DBL     Qib2b[4];			//����ϵ������������ϵ��̬��Ԫ��
	DBL		Vib_f[3];			//�������ϵ��������
	DBL		Qin2n[4];			//����ϵ����ڵ�������ϵ��̬��Ԫ��
	DBL     Vin_g[3];			//��������ϵ��������

	DBL		mk;					//Ȩ��ϵ����ĸ

	DBL		Vx[3];				//ʸ�����˵��Ƽ����м�ֵ
	DBL		Vy[3];				//ʸ�����˵��Ƽ����м�ֵ
	DBL		Vz[3];				//ʸ�����˵��Ƽ����м�ֵ

	DBL		Vb[3][3];			//ʸ�����˵��Ƽ����м�ֵ
	DBL		Vn[3][3];			//ʸ�����˵��Ƽ����м�ֵ


	DBL		Qin2ib[4];			//��׼�����̬��Ԫ��
	DBL		Qn2b[4];			//��׼�����̬��Ԫ��

	DBL		Euler_c[3];			//��׼���ŷ����
	DBL		Euler_t[3];			//��׼���ŷ����

	DBL		VPos_0[6];			//��׼��ʼʱ�������ٶȡ�λ��
	DBL		cL0;				//��ʼγ�ȶ�Ӧ����ֵ
	DBL		sL0;				//��ʼγ�ȶ�Ӧ����ֵ

	int		Align_Count;		//��׼�вο���Ϣ���ô���
	DBL		Ref_T;				//�ο�λ�á��ٶȸ�������

	DBL		winie[3];			//������ת���ٶ�ͶӰ
	
}ALIGN_DATA;

void Navigation_Init(DBL MINS[9]);					//��ʼ������,��׼ǰ����1��
extern void Navigation(NAVI_OUT *pNaviData, DBL IMU1[6], DBL IMU2[6]);		

extern void Init_KF_Filter();
extern void Kalman_Filter();

void Cal_KF_FNav(DBL Lat, DBL Vg[3], DBL DeltaVb[3], DBL Cbn[3][3], EARTH_DATA earth); 


extern void Matrix31_Mult(DBL *pa, DBL b[3], DBL c[3]);
extern int Matrix_Inverse(DBL *pM);
extern int Matrix_Inversesix(DBL pM[6][6]);   //��������
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
