#include "myNavigation.h"
#define nSimple 2			//惯导解算子样数

double avp0[10] = {0,0,0, 0,0,0, 0.6108652539,1.8849555922,400,0.01};		//初始信息
double imu_eb = 1*glv_dph;
double imu_db = 1*glv_mg;
double imu_web = 0.01*glv_dpsh;
double imu_wdb = 0.01*glv_mg;
double mea_dvn[3] = {0.01, 0.01, 0.01};
double mea_dpos[3] = {1, 1, 1};
double phi0[3] = {1*glv_deg, 1*glv_deg, 1*glv_deg};
double dvn0 = 0.1;
double dpos0[3] = {1, 1, 1};



int navmain();

