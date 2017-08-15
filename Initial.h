#ifndef __initial_h__
#define __initial_h__

#define ACD_MIUE 3.986005e14 /* Earth Force param */
#define ACD_PI 3.1415926535897932 /* Pi */
#define ACD_EPS 1.0e-14 /* a minimum number*/
#define ACD_Rearth 6378140.0 /* earth radus [m]*/
#define ACD_d2r 0.01745329251994
#define ACD_Wearth 7.27220521664e-5  /*����ת�����ٶ�*/
#define ACD_J2 1.08264e-3  /*J2���㶯ϵ��*/


/*======== ldb 2013-1-9 Start =========*/
extern double Wwheel_Pre[4][5];/* fei lun mai kuan */
extern unsigned long con1;
extern unsigned long con2;
extern unsigned long con3;
/*======== ldb 2013-1-9 End   =========*/

extern double gnm[7][7];  /* ��˹ϵ��*/
extern double hnm[7][7] ;  /* ��˹ϵ��*/

/*ϵͳ�������*/
extern double H;    /*����������ֲ���*/

/*���ǲ���*/
extern double ACS_IBf[3][3] ;
/*{{  193.3, 0,0},
{ 0, 193.5, 0},
{0, 0, 48.08}};	����ת������*/
extern double ACS_Msat;   /*��������*/
extern double ACS_MsatMin;
/*ƽ̨����*/
extern double ACS_Pgyrob1[3][3];/*{{0,0,-1},{0,-1,0},{-1,0,0},{-0.57735,-0.57735,-0.57735}}; ƽ̨���ݷ������Ҿ��� */
extern double ACS_GyroR1[3];/*ƽ̨��������ϵ��*/
extern double ACS_GyroRC1[3];/*ƽ̨���ݳ�ֵƯ��*/
extern unsigned char ACS_NoGyro1[3];  /*ƽ̨���ݲ�����Ч��־*/

/* ��������*/
extern double ACS_Pgyrob2[4][3];/*{{0,0,-1},{0,-1,0},{-1,0,0},{-0.57735,-0.57735,-0.57735}}; �������ݷ������Ҿ��� */
extern double ACS_GyroR2[4];/*������������ϵ��*/
extern double ACS_GyroRC2[4];/*�������ݳ�ֵƯ��*/
extern unsigned char ACS_NoGyro2[3];  /*�������ݲ�����Ч��־*/

/*��������*/
extern double ACS_StarsbA[4] ;/* {0.50161239751720732,  0.30224446310001007,  0.78737374568608542,  0.19255095893474417}; ��������A��������ת��Ԫ�� */
extern double ACS_quatstarRA[3];/*��������A����ϵ��*/
extern double ACS_StarALimit[3];/*��������A����ָ��ֱ���̫��\����\����ʸ���нǵ��ٽ�ֵ,��*/
extern double ACS_StarsbB[4];/*{0.50161239751720732,  -0.30224446310001007,  0.78737374568608542,  -0.19255095893474417}; ��������B��������ת��Ԫ�� */
extern double ACS_quatstarRB[3];/*��������B����ϵ��*/
extern double ACS_StarBLimit[3];/*��������B����ָ��ֱ���̫��\����\����ʸ���нǵ��ٽ�ֵ,��*/

/*����̫��������*/
extern double ACS_PFsunb1[3][3];/* ����̫��������1�������Ҿ��� */
/*double ACS_PFsunb1[3][3]={{0,0,-1.0},{0.5,-0.8660254,0},{-0.8660254,-0.5,0}};*/
/*double ACS_PFsunb1[3][3]={{0,0,-1.0},{0.5,-0.8660254,0},{-0.8660254,-0.5,0}};*/
extern double ACS_S1[3]; /*����̫��������1��װ��*/
extern double ACS_Fsun1Limit;/*����̫��������1�ӳ�������ֵ*/
extern double ACS_Fsun1R;/*����̫��������1����ϵ��*//*���ڱ�������̫��ģ�����ļ̵������ǣ���ֵ���ý�С*/
extern double ACS_PFsunb2[3][3];/* ����̫��������2�������Ҿ���  */
/*double ACS_PFsunb2[3][3]={{-0.5,-0.8660254,0.0},{0.0,0.0,-1.0},{0.8660254,-0.5,0.0}};*{{0.5,0.8660254,0},{0,0,1},{0.8660254,-0.5,0}}; ����̫��������2�������Ҿ���  */
/*double ACS_PFsunb2[3][3]={{0.8660254,0.5,0.0},{0.5,-0.8660254,0.0},{0.0,0.0,-1.0}};*/
extern double ACS_S2[3]; /*����̫��������2��װ��*/
extern double ACS_Fsun2Limit;/*����̫��������2�ӳ�������ֵ*/
extern double ACS_Fsun2R;/*����̫��������2����ϵ��*//*���ڱ�������̫��ģ�����ļ̵������ǣ���ֵ���ý�С*/

/*01̫��������*/
extern double ACS_Psun01b[3][3];/* 0-1̫���������������Ҿ���*/

/*���ٶȼ�*/
extern double ACS_Psb[3][3];  /*���ٶȼư�װ����*/
extern double ACS_AccelerR[3];/*���ٶȼ�����ϵ��*/
extern double ACS_AccelerRC[3];/*���ٶȼƳ�ƫ*/

/*��ǿ��*/
extern double ACS_PbBB[3][3];/* ��ǿ�Ʒ������Ҿ����� */
extern double ACS_MagoutR[3];

/*GPS���ջ�*/
extern double ACS_GPSR[14];
extern unsigned char ACS_GPSSixFlag ;/*GPSƽ��˲ʱ������־,1ƽ��,0˲ʱ*/
extern int GPSweek  ;
extern unsigned int GPSsecond  ;
extern double ACS_GPSDelay ;

extern double Distance;
extern double Theta_yx, Theta_yz;

/*����������*/
extern double ACS_SaeroC ;   /*ӭ�����������*/
extern double ACS_SaeroS ;   /*ӭ�����������*/
extern double ACS_Cd ;       /*��������ϵ��*/
extern double ACS_Cp[3] ;    /*����ϵ���������ĵ�����ѹ�ĵ�ʸ��*/
extern double ACS_r0 ;       /*200km�߶�*/
extern double ACS_density0;
extern double ACS_H0 ;       /*�ܶȱ��  m*/
extern double ACS_miu0 ;     /*�ܶȹ�ʽϵ��*/

/*�Ų���*/
extern unsigned char ACS_Number_mag ;  /*�ų�����*/
extern double ACS_Mdis[3];

/*̫����ѹ����*/
extern double ACS_Psunb[3][3];
extern double ACS_niu ;   /*̫��������淴��ϵ��*/
extern double ACS_Ssun ;   /*̫�����������*/
extern double ACS_Psun ;  /*̫����ѹ*/
extern double ACS_Cps[3];

/*�����������*//*20091203�ƽ�ϵͳ������Ʊ���*/
extern double ACS_DMsat[14];
extern double ACS_KFthurst[14];
extern double ACS_cosattis[14][3];/*double ACS_Cattis[14][3]={{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{-0.203,0.203,-1.354},{0.203,-0.203,-1.354},
{-0.203,-0.203,-1.354},{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{0.203,0.203,-1.354},
{-0.203,-0.203,-1.354},{0.203,0.203,-1.354},{-0.203,0.203,-1.354},{-0.203,0.203,-1.354},
{0.203,-0.203,-1.354},{0.203,-0.203,-1.354}};��������ƫ��ʸ��*/
extern double ACS_Cattis[14][3] ;

/*����ģ�Ͳ���*/
extern double ACS_Pwheelb[4][3];
extern double ACS_Jwheel[4];
extern double ACS_Twheel_MAX  ;/* N.m*/
/*��ƫ����������*/
extern double ACS_Rsa[3][3];
extern double ACS_DWa[3];
extern double ACS_Pab[3][3];

/*��������*/
extern double ACS_PMagb [3][3];


/*extern double �ڲ�����*/
/*extern double Time,Wwheel[4],Y[6],Wib[3],Quatib[4],Msat,ThTime[14];*/
extern double nine[9];    /*���������nine[9]:A,E,SI,DW,SW,DM,SU,theta,n*/
extern double Worbity;    /*�������ϵ�¹�����ٶ�y�������*/
extern double Tb[3];   /*���Ǻ�����*/
extern double Jingdu_GreenWitch_inertia;    /*���������Ӻ���ʱ��*/
extern double Bci[3];    /*��������λ�õĵشų�ʸ���ڹ���ϵ�µķ�����ʽ*/
extern double Td_geomagnetic[3];    /*�شŸ�������*/
extern double Yousun;  /*̫���ƾ�*/
extern double Aisun;  /*�Ƴཻ��*/
extern double SunInertia[3];  /*̫��ʸ���ڹ�������ϵ�µķ���*/
extern double Quatis[4];   /*̫������ϵ��Թ�������ϵ����̬��Ԫ��*/
extern double QuatsiA[4],QuatsiB[4];  /*��������A��B��װ����ϵ����������ϵ����̬��Ԫ��*/
extern double QuatsiA0[4],QuatsiB0[4];/*���������������*/
extern double Tc[3];       /*�����������*/
extern double Hbody[3];        /*���ǽǶ���*/
extern double TorqureCross[3]; /*������������*/
extern double Hwheel[3];       /*����ϵͳ�Ƕ���*/
extern double Sqe[4];   /*ƫ����Ԫ��*/
extern double Swe[3];   /*ƫ����ٶ�*/
extern double Wgyro1[4];  /*��������1�������ٶ�*/
extern double Wgyro10[4];/*��������1���۽��ٶ�*/
extern double Wgyro2[4];  /*��������2�������ٶ�*/
extern double Wgyro20[4];/*��������2���۽��ٶ�*/
extern double Fsun1out;   /*����̫��������1�������*/
extern double Fsun1out0; /*����̫��������1�������*/
extern unsigned char Fsun1Flag;/*����̫��������1������Ч��־λ*/
extern double Fsun2out;   /*����̫��������2�������*/
extern double Fsun2out0; /*����̫��������2�������*/
extern unsigned char Fsun2Flag;/*����̫��������2������Ч��־λ*/
extern unsigned char Sun01out[5];   /*0-1ʽ̫���������������*/
extern double Td_gravity[3];   /*�����ݶ�����*/
extern double Uwheel[3]; /*���ֵ�Դ*/
extern double F_aero[3];   /*������*/
extern double Td_aero[3];  /*��������*/
extern double F_Sun[3];   /*̫����ѹ��*/
extern double Td_Sun[3];   /*̫����ѹ����*/
extern double MomentumWheel[4];  /*���ֽǶ���*/
extern double as[3];   /*���ٶȼƲ������*/
extern double as0[3];   /*���ٶȼ��������*/
extern double Twheel[3];   /*�����������*/
extern double MagU[3];   /*�����������Ƶ�ѹ*/
extern double Tmag[3];   /*���������������*/
extern double Magout[3];  /*��ǿ�Ʋ������*/
extern double Magout0[3];/*��ǿ���������*/
extern unsigned char Umag;  /*����������Ч��־*/
extern double Earthi[3];  /*����ϵ�µ���ʸ��*/
extern double Mooni[3];  /*����ϵ������ʸ��*/
extern double GPSout[14]; /*GPS�������*/
extern double GPSout0[14];/*GPS�������*/
extern double Fb[3]; /*����ϵ���������ܺ�����*/
extern double F_Th[3];/*����ϵ��������������*/
extern double T_Th[3];/*����ϵ��������������*/
extern unsigned char EshadowFlag;/*��Ӱ��־λ*/
extern unsigned char StarAFlag;/*��������A������Ч��־λ*/
extern unsigned char StarBFlag;/*��������B������Ч��־λ*/
extern double Twheelb[3]; 	/*����ϵ�·��ֲ�������*/


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



/*extern double  �ⲿ����*/
extern double Time;
/*extern double Y[6]={6545298.2246124865,1154112.6737821733,0.00000000000000000,153.86531797261719,-872.61358048796740,7693.4020512707202};   */ /*����λ���ٶȳ�ʼ����*/
extern double Y[6];
extern double Wib[3];
extern double Quatib[4];
extern double Wwheel[4];
extern double ThTime[14];
extern double Mags[3];


extern double AOCS_Y[32];    /*ϵͳ�������*/
extern double AOCS_Y0[100];  /*ϵͳ�����������*/
extern unsigned char AOCS_Flag[9];/*ϵͳ�����־*/
extern unsigned char AOCS_Flag0[6];/*ϵͳ���������־*/
/*extern unsigned int OUT[35];*/
extern double AOCS_Xupdata[19];/*ϵͳ�������*/
extern double AOCS_Xinitial[14];/*ϵͳ�����ʼ����*/
extern double Msat; /*��������*/
extern unsigned char sFlag ;/*��ƫ������������־*/
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

extern unsigned char Thrust_Lv1Lv2  ;     /*������������LV1/LV2*/
extern unsigned char Thrust_Lv3Lv4  ;     /*������������LV3/LV4*/
extern unsigned char Thrust_Lv1 ,Thrust_Lv2 ,Thrust_Lv3 ,Thrust_Lv4 ;



#endif

