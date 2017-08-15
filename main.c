
#include"math.h"
#include"stdio.h"

#include"Initial.h"
#include"Functions.h"



/*动力学计算机全系统程序*/
void AOCS(double AOCS_Xupdata[19],double AOCS_Xinitial[14],int N,double AOCS_Y0[],unsigned char AOCS_Flag0[6],double AOCS_Y[],unsigned char AOCS_Flag[9]);
void mymain(void);


void mymain(void)
{


	N=N+1;
	/*输入数据*/
	//输出AOCS_Xinitial[14],分别为惯性系下卫星的位置速度矢量[6]，本体系下本体系相对惯性系的姿态角速度[3]
	//惯性系到本体系的姿态四元数[4]，以及卫星质量[1]。   整个AOCS_Xinitial[14]作为动力学递推计算的输入量
	Xinitial2all(Y,Wib,Quatib,Msat,AOCS_Xinitial);
	//输出AOCS_Xupdata[19],分别为卫星时间[1]，飞轮转速[4], 推力器开启时间[14]
	//整个AOCS_Xupdata[19]作为动力学的递推计算的输入量
	Xupdate2all(Time,Wwheel,ThTime,AOCS_Xupdata);
    //输入量：N未知，可能是递归参数；Y0为系统理论输出数据，Flag0为系统理论输出标志
	//Y为系统输出数据，包括太敏，陀螺，星敏，GPS测量输出， Flag为系统输出标志
	AOCS(AOCS_Xupdata,AOCS_Xinitial,N,AOCS_Y0,AOCS_Flag0,AOCS_Y,AOCS_Flag);
	/*  change(AOCS_Y,AOCS_Flag,OUT);*/

	Time=Time+H;

}


void AOCS(double AOCS_Xupdata[19],double AOCS_Xinitial[14],int N,double AOCS_Y0[100],unsigned char AOCS_Flag0[6],double AOCS_Y[],unsigned char AOCS_Flag[9])
{
	int i;
	double Worbity_CHK,nine_CHK[9];
	all2Xupdate(AOCS_Xupdata,&Time,Wwheel,ThTime); //输出飞轮转速[4]和推力器开启时间[14]
	all2Xinitial(AOCS_Xinitial,Y,Wib,Quatib,&Msat); // 输出卫星位置速度[6]，本体相对惯性系在本体系下的姿态角速度[3]，本体系相对惯性系的四元数[4]，卫星质量[1]

	/*质量更新模型*/
	MsatUPDATE(&Msat,ACS_DMsat,ThTime); //返回质量

	/*环境力矩*/
	//输入为卫星时间，太阳黄经，黄赤交角； 输出为格林尼治视恒星时角，太阳矢量[3]，测量系相对惯性系的姿态四元数[4]
	SunVectori(Time,&Yousun,&Aisun,&Jingdu_GreenWitch_inertia, SunInertia, Quatis);/*太阳矢量*/
	//输入为太阳矢量，本体系相对惯性系的四元数[3]，太阳帆板安装矩阵[3][3]，太阳帆板表面反射系数=0.5，太阳光照射面积=5，
	//太阳光压=9e-6，太阳光压压心到卫星质心的矢径[3]，输出本体系下太阳光压力[3]和光压力矩[3]
	SunTorque(SunInertia,Quatib,ACS_Psunb,ACS_niu,ACS_Ssun,ACS_Psun,ACS_Cps,F_Sun,Td_Sun);   /*太阳光压力矩模型*/

	/*执行机构*/
	//输入为推力器开启时间[14]，推力器推力系数[14]，推力器组喷嘴中心的方向余弦角[14][3]，推力器组偏心矢量[14][3]
	//输出本体系下推力器推力[3]和力矩[3]
	Fthrust(ThTime,ACS_KFthurst,ACS_cosattis,ACS_Cattis, F_Th, T_Th);
	//输入为飞轮转速[4]，飞轮转动惯量[4]为什么是一维矩阵? 飞轮安装矩阵[4][3], 输出本体系下飞轮组产生的力矩[3]
	Wheel(Wwheel,ACS_Jwheel,ACS_Pwheelb, Twheelb);
	Magnetic(Bci,Quatib,Mags,ACS_PMagb, Tmag);

	/*轨道递推*/
	for(i=0;i<3;i++)
        Fb[i]=F_Th[i]+F_aero[i]+F_Sun[i];/*计算本体系下卫星合外力*/
	orbJ2_rk4(Y,H,Fb,Quatib,Msat); /*轨道动力学模型，利用龙格库塔法递推，得到下一时刻的位置和速度，需要给出初始位置和速度*/
	PV2Six(Y,nine,&Worbity);  /*位置速度转换六根数模型，由位置和速度计算轨道六要素*/
	PV2Six_OLD(Y,nine_CHK,&Worbity_CHK); /*nine[9]:0A,1E,2SI,3DW,4SW,5DM,6SU,7theta,8n*/
	Iner2orbitquat(nine[3],nine[2],nine[6],Quatio);

	/*姿态递推*/
	QuatTatti_rk4(H,Wib,Quatib);/*姿态运动学，递推下一时刻的姿态四元素Qib  */
	for(i=0;i<3;i++)
		Tb[i]=Twheelb[i]+T_Th[i]+Td_geomagnetic[i]+Td_gravity[i]+Td_aero[i]+Td_Sun[i]+Tmag[i];   /*计算卫星合力矩*/
	//输入飞轮转动惯量[4]，飞轮转速[4]； 输出飞轮角动量[4]
	WheelMomentum(ACS_Jwheel,Wwheel,MomentumWheel);    /*飞轮角动量计算*/
	//输入飞轮角动量NMS[4](不知道是啥？)，飞轮安装矩阵[4][3]； 输出飞轮系统角动量[4]
	WheelSysToThree(MomentumWheel,ACS_Pwheelb,Hwheel);  /*飞轮系统角动量计算 	*/
	//输入飞轮系统角动量，本体相对惯性系的角速度[3]，卫星转动惯量[3][3]； 输出卫星角动量[3]
	formHB(Hwheel,Wib,ACS_IBf,Hbody);                /*卫星角动量计算*/
	//输入积分步长[1]，卫星角动量[3]，卫星合力矩[3]，卫星转动惯量[3][3]; 输出本体系相对惯性系的角速度[3]
	BiasattiDynamic_rk4(H,Hbody,Tb,ACS_IBf,Wib);   /*姿态动力学,带轮控系统的刚体动力学，递推下一时刻的角速度Wib*/

	/*天体矢量*/
	//输入为卫星时间+积分步长，太阳黄经，黄赤交角, 输出惯性系下的太阳矢量[3]，太阳系相对惯性系的四元数[4]，格林尼治视恒星时角[1]*/
	SunVectori(Time+H,&Yousun,&Aisun,&Jingdu_GreenWitch_inertia,SunInertia,Quatis);  /*太阳矢量模型*/
	//输入卫星时间+积分步长， 卫星位置速度[6]； 输出惯性系下月球矢量[3]
	MoonVectori(Time+H,Y,Mooni); /*月球矢量模型*/
	//输入卫星位置速度[6]； 输出惯性系下地球矢量[3]
	EarthVectori(Y,Earthi); /*地心矢量模型*/
    //输入卫星位置速度[6]，惯性系下太阳矢量[3]， 输出地影标志位：EshadowFlag=1时进入地影；
	EarthShadow(Y,SunInertia,&EshadowFlag);/*地影计算*/

	/*敏感器*/
    //输入太阳矢量[3]，本体系相对惯性系的四元数[4]，数字太敏1的安装矩阵[3][3]，数字太敏1的视场测量极值[1],数字太敏1的噪声系数[1]，地影标志位，数字太敏1测量有效标志位
	//输出理论值和带噪声值：太阳矢量在测量系yz和xz平面上的投影分别与z轴的夹角
	Fsun(SunInertia,Quatib,ACS_PFsunb1,ACS_Fsun1Limit,ACS_Fsun1R,EshadowFlag,&Fsun1Flag,&Fsun1out0,&Fsun1out);  /*数字太阳敏感器1测量模型*/
	Fsun(SunInertia,Quatib,ACS_PFsunb2,ACS_Fsun2Limit,ACS_Fsun2R,EshadowFlag,&Fsun2Flag,&Fsun2out0,&Fsun2out);  /*数字太阳敏感器2测量模型*/
	//输入太阳矢量，地球矢量，月球矢量，星敏A光轴指向分别与太阳、地球、月球矢量夹角的临界值[3]（单位是度°）,星敏A到星体的旋转四元数[4]
	//本体系相对惯性系的四元数[4]，星敏A噪声系数[3],星敏A测量有效标志位[1]； 输出星敏A/B安装
	quatstar(SunInertia,Earthi,Mooni,ACS_StarALimit,ACS_StarsbA,Quatib,ACS_quatstarRA,&StarAFlag,QuatsiA0,QuatsiA);    /*星敏感器A测量模型，输出测量四元素*/
	quatstar(SunInertia,Earthi,Mooni,ACS_StarBLimit,ACS_StarsbB,Quatib,ACS_quatstarRB,&StarBFlag,QuatsiB0,QuatsiB);    /*星敏感器B测量模型，输出测量四元素*/
	Gyro(Wib,ACS_Pgyrob1,ACS_NoGyro1,ACS_GyroR1,ACS_GyroRC1,Wgyro10,Wgyro1);        /*平台陀螺测量模型，输出测量角速度*/
	Gyro(Wib,ACS_Pgyrob2,ACS_NoGyro2,ACS_GyroR2,ACS_GyroRC2,Wgyro20,Wgyro2);        /*光纤陀螺测量模型，输出测量角速度*/
	Acceler(F_Th,F_aero,F_Sun,Msat,ACS_Psb,ACS_AccelerR,ACS_AccelerRC,as0,as);   /*加速度计测量模型，输出测量角加速度*/
	Magnet(Bci,Quatib,ACS_PbBB,ACS_MagoutR,Magout0,Magout);  /*磁强计模型，输出测量磁场强度*/

    Distance_in_Radar(r1,r2, &Distance);
    Angel_in_Camera(r1,r2,ACS_StarsbA,Quatib, &Theta_yz, &Theta_yx);

	/*系统输出*/
	Output(AOCS_Y0,AOCS_Flag0,AOCS_Y,AOCS_Flag);

}

int main()
{
    int i;
    for(i=0; i<10; i++)
    {
        mymain();
    }
    return 0;
}
