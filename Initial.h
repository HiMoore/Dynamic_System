#ifndef __initial_h__
#define __initial_h__

#define ACD_MIUE 3.986005e14 /* Earth Force param */
#define ACD_PI 3.1415926535897932 /* Pi */
#define ACD_EPS 1.0e-14 /* a minimum number*/
#define ACD_Rearth 6378140.0 /* earth radus [m]*/
#define ACD_d2r 0.01745329251994
#define ACD_Wearth 7.27220521664e-5  /*地球转动角速度*/
#define ACD_J2 1.08264e-3  /*J2项摄动系数*/


/*======== ldb 2013-1-9 Start =========*/
extern double Wwheel_Pre[4][5];/* fei lun mai kuan */
extern unsigned long con1;
extern unsigned long con2;
extern unsigned long con3;
/*======== ldb 2013-1-9 End   =========*/

extern double gnm[7][7];  /* 高斯系数*/
extern double hnm[7][7] ;  /* 高斯系数*/

/*系统仿真参数*/
extern double H;    /*龙格库塔积分步长*/

/*卫星参数*/
extern double ACS_IBf[3][3] ;
/*{{  193.3, 0,0},
{ 0, 193.5, 0},
{0, 0, 48.08}};	卫星转动惯量*/
extern double ACS_Msat;   /*卫星质量*/
extern double ACS_MsatMin;
/*平台陀螺*/
extern double ACS_Pgyrob1[3][3];/*{{0,0,-1},{0,-1,0},{-1,0,0},{-0.57735,-0.57735,-0.57735}}; 平台陀螺方向余弦矩阵 */
extern double ACS_GyroR1[3];/*平台陀螺噪声系数*/
extern double ACS_GyroRC1[3];/*平台陀螺常值漂移*/
extern unsigned char ACS_NoGyro1[3];  /*平台陀螺测量有效标志*/

/* 光纤陀螺*/
extern double ACS_Pgyrob2[4][3];/*{{0,0,-1},{0,-1,0},{-1,0,0},{-0.57735,-0.57735,-0.57735}}; 光纤陀螺方向余弦矩阵 */
extern double ACS_GyroR2[4];/*光纤陀螺噪声系数*/
extern double ACS_GyroRC2[4];/*光纤陀螺常值漂移*/
extern unsigned char ACS_NoGyro2[3];  /*光纤陀螺测量有效标志*/

/*星敏感器*/
extern double ACS_StarsbA[4] ;/* {0.50161239751720732,  0.30224446310001007,  0.78737374568608542,  0.19255095893474417}; 星敏感器A到星体旋转四元数 */
extern double ACS_quatstarRA[3];/*星敏感器A噪声系数*/
extern double ACS_StarALimit[3];/*星敏感器A光轴指向分别与太阳\地球\月球矢量夹角的临界值,度*/
extern double ACS_StarsbB[4];/*{0.50161239751720732,  -0.30224446310001007,  0.78737374568608542,  -0.19255095893474417}; 星敏感器B到星体旋转四元数 */
extern double ACS_quatstarRB[3];/*星敏感器B噪声系数*/
extern double ACS_StarBLimit[3];/*星敏感器B光轴指向分别与太阳\地球\月球矢量夹角的临界值,度*/

/*数字太阳敏感器*/
extern double ACS_PFsunb1[3][3];/* 数字太阳敏感器1方向余弦矩阵 */
/*double ACS_PFsunb1[3][3]={{0,0,-1.0},{0.5,-0.8660254,0},{-0.8660254,-0.5,0}};*/
/*double ACS_PFsunb1[3][3]={{0,0,-1.0},{0.5,-0.8660254,0},{-0.8660254,-0.5,0}};*/
extern double ACS_S1[3]; /*数字太阳敏感器1安装面*/
extern double ACS_Fsun1Limit;/*数字太阳敏感器1视场测量极值*/
extern double ACS_Fsun1R;/*数字太阳敏感器1噪声系数*//*出于保护数字太阳模拟器的继电器考虑，此值设置较小*/
extern double ACS_PFsunb2[3][3];/* 数字太阳敏感器2方向余弦矩阵  */
/*double ACS_PFsunb2[3][3]={{-0.5,-0.8660254,0.0},{0.0,0.0,-1.0},{0.8660254,-0.5,0.0}};*{{0.5,0.8660254,0},{0,0,1},{0.8660254,-0.5,0}}; 数字太阳敏感器2方向余弦矩阵  */
/*double ACS_PFsunb2[3][3]={{0.8660254,0.5,0.0},{0.5,-0.8660254,0.0},{0.0,0.0,-1.0}};*/
extern double ACS_S2[3]; /*数字太阳敏感器2安装面*/
extern double ACS_Fsun2Limit;/*数字太阳敏感器2视场测量极值*/
extern double ACS_Fsun2R;/*数字太阳敏感器2噪声系数*//*出于保护数字太阳模拟器的继电器考虑，此值设置较小*/

/*01太阳敏感器*/
extern double ACS_Psun01b[3][3];/* 0-1太阳敏感器方向余弦矩阵*/

/*加速度计*/
extern double ACS_Psb[3][3];  /*加速度计安装矩阵*/
extern double ACS_AccelerR[3];/*加速度计噪声系数*/
extern double ACS_AccelerRC[3];/*加速度计常偏*/

/*磁强计*/
extern double ACS_PbBB[3][3];/* 磁强计方向余弦矩阵逆 */
extern double ACS_MagoutR[3];

/*GPS接收机*/
extern double ACS_GPSR[14];
extern unsigned char ACS_GPSSixFlag ;/*GPS平均瞬时根数标志,1平均,0瞬时*/
extern int GPSweek  ;
extern unsigned int GPSsecond  ;
extern double ACS_GPSDelay ;

extern double Distance;
extern double Theta_yx, Theta_yz;

/*气动力参数*/
extern double ACS_SaeroC ;   /*迎风面积常数项*/
extern double ACS_SaeroS ;   /*迎风面积周期项*/
extern double ACS_Cd ;       /*大气阻力系数*/
extern double ACS_Cp[3] ;    /*本体系下卫星质心到气动压心的矢径*/
extern double ACS_r0 ;       /*200km高度*/
extern double ACS_density0;
extern double ACS_H0 ;       /*密度标高  m*/
extern double ACS_miu0 ;     /*密度公式系数*/

/*磁参数*/
extern unsigned char ACS_Number_mag ;  /*磁场阶数*/
extern double ACS_Mdis[3];

/*太阳光压参数*/
extern double ACS_Psunb[3][3];
extern double ACS_niu ;   /*太阳帆板表面反射系数*/
extern double ACS_Ssun ;   /*太阳光照射面积*/
extern double ACS_Psun ;  /*太阳光压*/
extern double ACS_Cps[3];

/*推力器组参数*//*20091203推进系统初样设计报告*/
extern double ACS_DMsat[14];
extern double ACS_KFthurst[14];
extern double ACS_cosattis[14][3];/*double ACS_Cattis[14][3]={{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{-0.203,0.203,-1.354},{0.203,-0.203,-1.354},
{-0.203,-0.203,-1.354},{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{0.203,0.203,-1.354},
{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{-0.203,0.203,-1.354},{-0.203,0.203,-1.354},
{0.203,-0.203,-1.354},{0.203,-0.203,-1.354}};推力器组偏心矢量*/
extern double ACS_Cattis[14][3] ;

/*飞轮模型参数*/
extern double ACS_Pwheelb[4][3];
extern double ACS_Jwheel[4];
extern double ACS_Twheel_MAX  ;/* N.m*/
/*调偏流机构参数*/
extern double ACS_Rsa[3][3];
extern double ACS_DWa[3];
extern double ACS_Pab[3][3];

/*磁力矩器*/
extern double ACS_PMagb [3][3];


/*extern double 内部定义*/
/*extern double Time,Wwheel[4],Y[6],Wib[3],Quatib[4],Msat,ThTime[14];*/
extern double nine[9];    /*轨道根数，nine[9]:A,E,SI,DW,SW,DM,SU,theta,n*/
extern double Worbity;    /*轨道坐标系下轨道角速度y方向分量*/
extern double Tb[3];   /*卫星合力矩*/
extern double Jingdu_GreenWitch_inertia;    /*格林尼制视恒星时角*/
extern double Bci[3];    /*卫星所处位置的地磁场矢量在惯性系下的分量形式*/
extern double Td_geomagnetic[3];    /*地磁干扰力矩*/
extern double Yousun;  /*太阳黄经*/
extern double Aisun;  /*黄赤交角*/
extern double SunInertia[3];  /*太阳矢量在惯性坐标系下的分量*/
extern double Quatis[4];   /*太阳坐标系相对惯性坐标系的姿态四元素*/
extern double QuatsiA[4],QuatsiB[4];  /*星敏感器A、B安装坐标系到惯性坐标系的姿态四元素*/
extern double QuatsiA0[4],QuatsiB0[4];/*星敏感器理论输出*/
extern double Tc[3];       /*三轴控制力矩*/
extern double Hbody[3];        /*卫星角动量*/
extern double TorqureCross[3]; /*卫星陀螺力矩*/
extern double Hwheel[3];       /*飞轮系统角动量*/
extern double Sqe[4];   /*偏差四元素*/
extern double Swe[3];   /*偏差角速度*/
extern double Wgyro1[4];  /*光纤陀螺1测量角速度*/
extern double Wgyro10[4];/*光纤陀螺1理论角速度*/
extern double Wgyro2[4];  /*光纤陀螺2测量角速度*/
extern double Wgyro20[4];/*光纤陀螺2理论角速度*/
extern double Fsun1out;   /*数字太阳敏感器1测量输出*/
extern double Fsun1out0; /*数字太阳敏感器1理论输出*/
extern unsigned char Fsun1Flag;/*数字太阳敏感器1测量有效标志位*/
extern double Fsun2out;   /*数字太阳敏感器2测量输出*/
extern double Fsun2out0; /*数字太阳敏感器2理论输出*/
extern unsigned char Fsun2Flag;/*数字太阳敏感器2测量有效标志位*/
extern unsigned char Sun01out[5];   /*0-1式太阳敏感器测量输出*/
extern double Td_gravity[3];   /*重力梯度力矩*/
extern double Uwheel[3]; /*飞轮电源*/
extern double F_aero[3];   /*气动力*/
extern double Td_aero[3];  /*气动力矩*/
extern double F_Sun[3];   /*太阳光压力*/
extern double Td_Sun[3];   /*太阳光压力矩*/
extern double MomentumWheel[4];  /*飞轮角动量*/
extern double as[3];   /*加速度计测量输出*/
extern double as0[3];   /*加速度计理论输出*/
extern double Twheel[3];   /*飞轮输出力矩*/
extern double MagU[3];   /*磁力矩器控制电压*/
extern double Tmag[3];   /*磁力矩器输出力矩*/
extern double Magout[3];  /*磁强计测量输出*/
extern double Magout0[3];/*磁强计理论输出*/
extern unsigned char Umag;  /*磁力矩器有效标志*/
extern double Earthi[3];  /*惯性系下地心矢量*/
extern double Mooni[3];  /*惯性系下月球矢量*/
extern double GPSout[14]; /*GPS测量输出*/
extern double GPSout0[14];/*GPS理论输出*/
extern double Fb[3]; /*本体系下卫星所受合外力*/
extern double F_Th[3];/*本体系下推力器合推力*/
extern double T_Th[3];/*本体系下推力器合力矩*/
extern unsigned char EshadowFlag;/*地影标志位*/
extern unsigned char StarAFlag;/*星敏感器A测量有效标志位*/
extern unsigned char StarBFlag;/*星敏感器B测量有效标志位*/
extern double Twheelb[3]; 	/*本体系下飞轮产生力矩*/


extern struct OUTPUTstr
{
	float strY[6];
	float strNine[9];
	float strQuatib[4];
	float strWib[3];
	float strMsat;
	float strF_Th[3];
	float strT_Th[3];
	float strTwheelb[3];
	float strFsun1out0;
	float strFsun2out0;
	float strQuatsiA0[4];
	float strQuatsiB0[4];
	float strWgyro10[4];
	float strWgyro20[4];
	float stras0[3];
	float strMagout0[3];
	float strGPSout0[14];
	float strSunInertia[3];
	float strMooni[3];
	float strEarthi[3];
	float strBci[3];
	float strTd_geomagnetic[3];
	float strTd_gravity[3];
	float strTd_aero[3];
	float strTd_Sun[3];
	float strF_aero[3];
	float strF_Sun[3];
	float fEuler[3];
	float fHbody[3];
	unsigned char strEshadowFlag;
	unsigned char strFsun1Flag;
	unsigned char strFsun2Flag;
	unsigned char strSun01;
	unsigned char strStarAFlag;
	unsigned char strStarBFlag;
} outputLAN;



/*extern double  外部定义*/
extern double Time;
/*extern double Y[6]={6545298.2246124865,1154112.6737821733,0.00000000000000000,153.86531797261719,-872.61358048796740,7693.4020512707202};   */ /*卫星位置速度初始条件*/
extern double Y[6];
extern double Wib[3];
extern double Quatib[4];
extern double Wwheel[4];
extern double ThTime[14];
extern double Mags[3];


extern double AOCS_Y[32];    /*系统输出数据*/
extern double AOCS_Y0[100];  /*系统理论输出数据*/
extern unsigned char AOCS_Flag[9];/*系统输出标志*/
extern unsigned char AOCS_Flag0[6];/*系统理论输出标志*/
/*extern unsigned int OUT[35];*/
extern double AOCS_Xupdata[19];/*系统输入变量*/
extern double AOCS_Xinitial[14];/*系统输入初始条件*/
extern double Msat; /*卫星质量*/
extern unsigned char sFlag ;/*调偏流机构启动标志*/
extern int N;
extern double Euler[3];
extern double Quatio[4];

extern double para_WgyroX;
extern double para_WgyroY;
extern double para_WgyroZ;
extern double para_WgyroS;
extern double para_Kx ;
extern double para_Ky ;
extern double para_Kz ;
extern double para_Ks ;

extern long int cntN ,cntN_A ;


extern double Wheel_compute   ;
extern unsigned char THRUST_OPEN ;
extern unsigned char THRUST_CLOSE ;

extern unsigned char Thrust_Lv1Lv2  ;     /*发动机自锁阀LV1/LV2*/
extern unsigned char Thrust_Lv3Lv4  ;     /*发动机自锁阀LV3/LV4*/
extern unsigned char Thrust_Lv1 ,Thrust_Lv2 ,Thrust_Lv3 ,Thrust_Lv4 ;



#endif

