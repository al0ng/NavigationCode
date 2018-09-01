/************************************************************************/
/* 惯导解算函数库                                                        */
/*     进行不同子样数的惯导解算                                           */
/* V1.1  2018-08-23  修改部分函數錯誤     			                    */
/************************************************************************/
#include "myNavigation.h"
#include <math.h>
#include <stdio.h>

insdata ins;	//惯导解算数据结构体
avpdata avp;	//输出数据结构体
EARTH earth;	//地球参数结构体

//==============================================
//  姿态阵转四元数
//==============================================
int m2qua(matrix* m, matrix* q)
{
	double C11,C12,C13,C21,C22,C23,C31,C32,C33;
	C11 = matrix_get(m,0,0);
	C12 = matrix_get(m,0,1);
	C13 = matrix_get(m,0,2);
	C21 = matrix_get(m,1,0);
	C22 = matrix_get(m,1,1);
	C23 = matrix_get(m,1,2);
	C31 = matrix_get(m,2,0);
	C32 = matrix_get(m,2,1);
	C33 = matrix_get(m,2,2);

	matrix_set(q,0,0, 0.5*sqrt(abs(1+C11+C22+C33)));
	matrix_set(q,1,0, ((C32-C23<0)?-1:1)*0.5*sqrt(abs(1+C11-C22-C33)));
	matrix_set(q,2,0, ((C13-C31<0)?-1:1)*0.5*sqrt(abs(1-C11+C22-C33)));
	matrix_set(q,3,0, ((C21-C12<0)?-1:1)*0.5*sqrt(abs(1-C11-C22+C33)));

	return 0;
}

//==============================================
//  四元数转姿态角
//==============================================
int q2att(matrix* qnb, double att[])
{
	double q11,q12,q13,q14,q22,q23,q24,q33,q34,q44;
	double C12,C22,C31,C32,C33;
	q11 = matrix_get(qnb,0,0)*matrix_get(qnb,0,0);
	q12 = matrix_get(qnb,0,0)*matrix_get(qnb,1,0);
	q13 = matrix_get(qnb,0,0)*matrix_get(qnb,2,0);
	q14 = matrix_get(qnb,0,0)*matrix_get(qnb,3,0);
	q22 = matrix_get(qnb,1,0)*matrix_get(qnb,1,0);
	q23 = matrix_get(qnb,1,0)*matrix_get(qnb,2,0);
	q24 = matrix_get(qnb,1,0)*matrix_get(qnb,3,0);
	q33 = matrix_get(qnb,2,0)*matrix_get(qnb,2,0);
	q34 = matrix_get(qnb,2,0)*matrix_get(qnb,3,0);
	q44 = matrix_get(qnb,3,0)*matrix_get(qnb,3,0);
	C12=2*(q23-q14);
	C22=q11-q22+q33-q44;
	C31=2*(q24-q13); C32=2*(q34+q12); C33=q11-q22-q33+q44;

	att[0] = asin(C32); 
	att[1] = atan2(-C31,C33);
	att[2] = atan2(-C12,C22);
	
	return 0;
}

//==============================================
//  姿态角转四元数
//==============================================
int a2qua(double *att, matrix *qnb)
{
	double sp,sr,sy, cp,cr,cy;
	sp = sin(att[0]/2); sr = sin(att[1]/2); sy = sin(att[2]/2); 
	cp = cos(att[0]/2); cr = cos(att[1]/2); cy = cos(att[2]/2); 
	matrix_set(qnb,0,0, cp*cr*cy - sp*sr*sy);
	matrix_set(qnb,1,0, sp*cr*cy - cp*sr*sy);
	matrix_set(qnb,2,0, cp*sr*cy + sp*cr*sy);
	matrix_set(qnb,3,0, cp*cr*sy + sp*sr*cy);
	return 0;
}

//==============================================
//  旋转矢量转四元数
//==============================================
int rv2q(matrix *rvm, matrix* q)
{
	double rv[3] = {matrix_get(rvm,0,0),matrix_get(rvm,1,0),matrix_get(rvm,2,0)};
	double norm_rv;
	norm_rv = sqrt(rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2]);
	matrix_set(q,0,0, cos(norm_rv/2));
	matrix_set(q,1,0, sin(norm_rv/2)*rv[0]/norm_rv);
	matrix_set(q,2,0, sin(norm_rv/2)*rv[1]/norm_rv);
	matrix_set(q,3,0, sin(norm_rv/2)*rv[2]/norm_rv);

	return 0;
}

//==============================================
//  四元数转姿态阵
//==============================================
int q2mat(matrix* qnb, matrix* Cnb)
{
	double q1,q2,q3,q4;
	q1 = matrix_get(qnb,0,0);
	q2 = matrix_get(qnb,1,0);
	q3 = matrix_get(qnb,2,0);
	q4 = matrix_get(qnb,3,0);

	matrix_set(Cnb, 0,0, q1*q1+q2*q2-q3*q3-q4*q4);
	matrix_set(Cnb, 0,1, 2*(q2*q3-q1*q4));
	matrix_set(Cnb, 0,2, 2*(q2*q4+q1*q3));
	matrix_set(Cnb, 1,0, 2*(q2*q3+q1*q4));
	matrix_set(Cnb, 1,1, q1*q1-q2*q2+q3*q3-q4*q4);
	matrix_set(Cnb, 1,2, 2*(q3*q4-q1*q2));
	matrix_set(Cnb, 2,0, 2*(q2*q4-q1*q3));
	matrix_set(Cnb, 2,1, 2*(q3*q4+q1*q2));
	matrix_set(Cnb, 2,2, q1*q1-q2*q2-q3*q3+q4*q4);

	return 0;
}

//==============================================
//  四元数乘法
//==============================================
int qmul(matrix* p, matrix* q, matrix* qm)
{
	double p1,p2,p3,p4,q1,q2,q3,q4;

	p1 = matrix_get(p,0,0);
	p2 = matrix_get(p,1,0);
	p3 = matrix_get(p,2,0);
	p4 = matrix_get(p,3,0);
	q1 = matrix_get(q,0,0);
	q2 = matrix_get(q,1,0);
	q3 = matrix_get(q,2,0);
	q4 = matrix_get(q,3,0);

	matrix_set(qm,0,0, p1*q1-p2*q2-p3*q3-p4*q4);
	matrix_set(qm,1,0, p1*q2 + p2*q1 + p3*q4 - p4*q3);
	matrix_set(qm,2,0, p1*q3 + p3*q1 + p4*q2 - p2*q4);
	matrix_set(qm,3,0, p1*q4 + p4*q1 + p2*q3 - p3*q2);

	return 0;
}

//==============================================
//  四元数乘向量
//==============================================
int qmulv(matrix *qnb, matrix *vi, matrix *vo)
{
	matrix *Cnb = matrix_calloc(3,3);
	q2mat(qnb, Cnb);

	matrix_mul(Cnb, vi, vo);
	matrix_free(Cnb);
	return 0;
}

//==============================================
//  四元数规范化
//==============================================
int qnormlz(matrix *qnb)
{
	double norm_inv;
	norm_inv = 1/matrix_norm(qnb);
	matrix_mul(qnb, norm_inv, qnb);
	return 0;
}

//==============================================
//  惯导解算初始化
//==============================================
int insinitial(double *avp0, int nn)
{
	ins.ts = (float)avp0[9];
	ins.nn = nn;
	ins.nnts = ins.nn*ins.ts;
	ins.pos = matrix_calloc(3,1);
	ins.vn  = matrix_calloc(3,1);
	ins.vn1  = matrix_calloc(3,1);
	ins.qnb = matrix_calloc(4,1);
	matrix_set(ins.vn, 0,0, avp0[3]);
	matrix_set(ins.vn, 1,0, avp0[4]);
	matrix_set(ins.vn, 2,0, avp0[5]);
	matrix_set(ins.pos, 0,0, avp0[6]);
	matrix_set(ins.pos, 1,0, avp0[7]);
	matrix_set(ins.pos, 2,0, avp0[8]);
	a2qua(avp0, ins.qnb);

	for (int i=0; i<ins.nn; i++)
	{
		ins.dvm[i] = matrix_calloc(3,1);
		ins.dwm[i] = matrix_calloc(3,1);
	}
	ins.vmm = matrix_calloc(3,1);
	ins.wmm = matrix_calloc(3,1);

	earth.gn = matrix_calloc(3,1);
	earth.gcc = matrix_calloc(3,1);
	earth.wnie = matrix_calloc(3,1);
	earth.wnen = matrix_calloc(3,1);
	earth.wnin = matrix_calloc(3,1);
	earth.wnien = matrix_calloc(3,1);

	return 0;
}

//==============================================
//  地球参数更新
//==============================================
void earthupdate()
{
	matrix *vxw = matrix_calloc(3,1);
	double sq, sq2;
	earth.sl = sin(matrix_get(ins.pos,0,0));
	earth.cl = cos(matrix_get(ins.pos,0,0));
	earth.tl = earth.sl/earth.cl;
	earth.sl2 = earth.sl*earth.sl;
	earth.sl4 = earth.sl2*earth.sl2;
	sq = 1-glv_e2*earth.sl2;
	sq2 = sqrt(sq);
	earth.RMh = glv_Re*(1-glv_e2)/sq/sq2+matrix_get(ins.pos,2,0);
	earth.RNh = glv_Re/sq2+matrix_get(ins.pos,2,0);
	earth.clRNh = earth.cl*earth.RNh;
	matrix_set(earth.wnie,1,0, glv_wie*earth.cl);
	matrix_set(earth.wnie,2,0, glv_wie*earth.sl);
	matrix_set(earth.wnen,0,0, -matrix_get(ins.vn,1,0)/earth.RMh);
	matrix_set(earth.wnen,1,0, matrix_get(ins.vn,0,0)/earth.RNh);
	matrix_set(earth.wnen,2,0, matrix_get(ins.vn,0,0)/earth.RNh*earth.tl);
	matrix_add(earth.wnie, earth.wnen, earth.wnin);
	matrix_add(earth.wnie, earth.wnin, earth.wnien);
	earth.g = glv_g0*(1+5.27094e-3*earth.sl2+2.32718e-5*earth.sl4)-3.086e-6*matrix_get(ins.pos,2,0);
	matrix_set(earth.gn, 2,0, -earth.g);
	matrix_crossmul(ins.vn, earth.wnien, vxw);
	matrix_add(earth.gn, vxw, earth.gcc);
	matrix_free(vxw);
}

//==============================================
//  姿态更新
//==============================================
int attupdate()
{
	matrix *phi = matrix_calloc(3,1);
	matrix *dphi = matrix_calloc(3,1);
	matrix *wnints = matrix_calloc(3,1);
	matrix *qbb = matrix_calloc(4,1);
	matrix *qnn = matrix_calloc(4,1);

	switch (ins.nn)
	{
	case 1:
		matrix_memcpy(phi, ins.wmm);
		break;
	case 2:
		matrix_crossmul(ins.dwm[0], ins.dwm[1], dphi);
		matrix_mul(dphi, 2.0/3.0, dphi);
		matrix_add(phi, ins.wmm);
		matrix_add(phi, dphi);
		break;
	case 4:
		break;
	}
	
	matrix_mul(earth.wnin, -ins.nnts, wnints);
	rv2q(phi, qbb);
	rv2q(wnints, qnn);
	qmul(ins.qnb,qbb,ins.qnb);
	qmul(qnn,ins.qnb,ins.qnb);
	qnormlz(ins.qnb);

	matrix_free(qnn);
	matrix_free(qbb);
	matrix_free(wnints);
	matrix_free(phi);
	matrix_free(dphi);
	return 0;
}

//==============================================
//  速度更新
//==============================================
int vnupdate()
{
	matrix* vm1 = matrix_calloc(3,1);
	matrix* vm2 = matrix_calloc(3,1);
	matrix* cm1 = matrix_calloc(3,1);
	matrix* cm2 = matrix_calloc(3,1);
	matrix *rotm = matrix_calloc(3,1);
	matrix *rotm_half = matrix_calloc(3,1);
	matrix *scullm = matrix_calloc(3,1);
	matrix *gccts = matrix_calloc(3,1);
	matrix *wnints2 = matrix_calloc(3,1);
	matrix *qnn = matrix_calloc(4,1);

	matrix_crossmul(ins.wmm,ins.vmm,rotm);
	matrix_mul(rotm, 0.5, rotm_half);
	matrix_mul(earth.gcc, ins.nnts, gccts);
	matrix_mul(earth.wnin, -ins.nnts/2, wnints2);
	switch (ins.nn)
	{
	case 1:
		
		break;
	case 2:
		matrix_crossmul(ins.dvm[0], ins.dwm[1], cm1);
		matrix_crossmul(ins.dvm[1], ins.dwm[0], cm2);
		matrix_sub(cm1, cm2, scullm);
		matrix_mul(scullm, 2.0/3.0, scullm);
		break;
	case 4:
		break;
	}
	matrix_add(ins.vmm, rotm_half, vm1);
	matrix_add(vm1, scullm);
	rv2q(wnints2, qnn);
	qmulv(ins.qnb, vm1, vm2);
	qmulv(qnn, vm2, vm1);
	matrix_memcpy(ins.vn1, ins.vn);
	matrix_add(ins.vn, vm1);
	matrix_add(ins.vn, gccts);

	matrix_free(vm1);
	matrix_free(vm2);
	matrix_free(cm1);
	matrix_free(cm2);
	matrix_free(rotm);
	matrix_free(rotm_half);
	matrix_free(scullm);
	matrix_free(gccts);
	matrix_free(wnints2);
	matrix_free(qnn);
	return 0;
}

//==============================================
//  位置更新
//==============================================
int posupdate()
{
	matrix *dpos = matrix_calloc(3,1);
	matrix_add(ins.vn1, ins.vn);
	matrix_mul(ins.vn1, 0.5, ins.vn1);
	matrix_set(dpos,0,0, matrix_get(ins.vn1,1,0)/earth.RMh);
	matrix_set(dpos,1,0, matrix_get(ins.vn1,0,0)/earth.clRNh);
	matrix_set(dpos,2,0, matrix_get(ins.vn1,2,0));
	matrix_mul(dpos, ins.nnts, dpos);
	matrix_add(ins.pos, dpos);

	matrix_free(dpos);
	return 0;
}

//==============================================
//  惯导解算更新
//==============================================
void insupdate(matrix *wm, matrix *vm)
{
	if (ins.nn>1)
	{
		for (int k=0; k<ins.nn; k++)
		{
			matrix_getv(wm,k,ins.dwm[k]);
			matrix_getv(vm,k,ins.dvm[k]);
		}
	}
	matrix_sumrow(vm,ins.vmm);
	matrix_sumrow(wm,ins.wmm);
	earthupdate();
	vnupdate();
	posupdate();
	attupdate();
}

//==============================================
//  ins数据转成姿态速度位置信息
//==============================================
void ins2avp()
{
	q2att(ins.qnb, avp.att);
	avp.vn[0]  = matrix_get(ins.vn, 0,0);
	avp.vn[1]  = matrix_get(ins.vn, 1,0);
	avp.vn[2]  = matrix_get(ins.vn, 2,0);
	avp.pos[0] = matrix_get(ins.pos, 0,0);
	avp.pos[1] = matrix_get(ins.pos, 1,0);
	avp.pos[2] = matrix_get(ins.pos, 2,0);
}

