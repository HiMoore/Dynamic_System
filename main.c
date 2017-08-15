
#include"math.h"
#include"stdio.h"

#include"Initial.h"
#include"Functions.h"



/*����ѧ�����ȫϵͳ����*/
void AOCS(double AOCS_Xupdata[19],double AOCS_Xinitial[14],int N,double AOCS_Y0[],unsigned char AOCS_Flag0[6],double AOCS_Y[],unsigned char AOCS_Flag[9]);
void mymain(void);


void mymain(void)
{


	N=N+1;
	/*��������*/
	//���AOCS_Xinitial[14],�ֱ�Ϊ����ϵ�����ǵ�λ���ٶ�ʸ��[6]������ϵ�±���ϵ��Թ���ϵ����̬���ٶ�[3]
	//����ϵ������ϵ����̬��Ԫ��[4]���Լ���������[1]��   ����AOCS_Xinitial[14]��Ϊ����ѧ���Ƽ����������
	Xinitial2all(Y,Wib,Quatib,Msat,AOCS_Xinitial);
	//���AOCS_Xupdata[19],�ֱ�Ϊ����ʱ��[1]������ת��[4], ����������ʱ��[14]
	//����AOCS_Xupdata[19]��Ϊ����ѧ�ĵ��Ƽ����������
	Xupdate2all(Time,Wwheel,ThTime,AOCS_Xupdata);
    //��������Nδ֪�������ǵݹ������Y0Ϊϵͳ����������ݣ�Flag0Ϊϵͳ���������־
	//YΪϵͳ������ݣ�����̫�������ݣ�������GPS��������� FlagΪϵͳ�����־
	AOCS(AOCS_Xupdata,AOCS_Xinitial,N,AOCS_Y0,AOCS_Flag0,AOCS_Y,AOCS_Flag);
	/*  change(AOCS_Y,AOCS_Flag,OUT);*/

	Time=Time+H;

}


void AOCS(double AOCS_Xupdata[19],double AOCS_Xinitial[14],int N,double AOCS_Y0[100],unsigned char AOCS_Flag0[6],double AOCS_Y[],unsigned char AOCS_Flag[9])
{
	int i;
	double Worbity_CHK,nine_CHK[9];
	all2Xupdate(AOCS_Xupdata,&Time,Wwheel,ThTime); //�������ת��[4]������������ʱ��[14]
	all2Xinitial(AOCS_Xinitial,Y,Wib,Quatib,&Msat); // �������λ���ٶ�[6]��������Թ���ϵ�ڱ���ϵ�µ���̬���ٶ�[3]������ϵ��Թ���ϵ����Ԫ��[4]����������[1]

	/*��������ģ��*/
	MsatUPDATE(&Msat,ACS_DMsat,ThTime); //��������

	/*��������*/
	//����Ϊ����ʱ�䣬̫���ƾ����Ƴཻ�ǣ� ���Ϊ���������Ӻ���ʱ�ǣ�̫��ʸ��[3]������ϵ��Թ���ϵ����̬��Ԫ��[4]
	SunVectori(Time,&Yousun,&Aisun,&Jingdu_GreenWitch_inertia, SunInertia, Quatis);/*̫��ʸ��*/
	//����Ϊ̫��ʸ��������ϵ��Թ���ϵ����Ԫ��[3]��̫�����尲װ����[3][3]��̫��������淴��ϵ��=0.5��̫�����������=5��
	//̫����ѹ=9e-6��̫����ѹѹ�ĵ��������ĵ�ʸ��[3]���������ϵ��̫����ѹ��[3]�͹�ѹ����[3]
	SunTorque(SunInertia,Quatib,ACS_Psunb,ACS_niu,ACS_Ssun,ACS_Psun,ACS_Cps,F_Sun,Td_Sun);   /*̫����ѹ����ģ��*/

	/*ִ�л���*/
	//����Ϊ����������ʱ��[14]������������ϵ��[14]�����������������ĵķ������ҽ�[14][3]����������ƫ��ʸ��[14][3]
	//�������ϵ������������[3]������[3]
	Fthrust(ThTime,ACS_KFthurst,ACS_cosattis,ACS_Cattis, F_Th, T_Th);
	//����Ϊ����ת��[4]������ת������[4]Ϊʲô��һά����? ���ְ�װ����[4][3], �������ϵ�·��������������[3]
	Wheel(Wwheel,ACS_Jwheel,ACS_Pwheelb, Twheelb);
	Magnetic(Bci,Quatib,Mags,ACS_PMagb, Tmag);

	/*�������*/
	for(i=0;i<3;i++)
        Fb[i]=F_Th[i]+F_aero[i]+F_Sun[i];/*���㱾��ϵ�����Ǻ�����*/
	orbJ2_rk4(Y,H,Fb,Quatib,Msat); /*�������ѧģ�ͣ�����������������ƣ��õ���һʱ�̵�λ�ú��ٶȣ���Ҫ������ʼλ�ú��ٶ�*/
	PV2Six(Y,nine,&Worbity);  /*λ���ٶ�ת��������ģ�ͣ���λ�ú��ٶȼ�������Ҫ��*/
	PV2Six_OLD(Y,nine_CHK,&Worbity_CHK); /*nine[9]:0A,1E,2SI,3DW,4SW,5DM,6SU,7theta,8n*/
	Iner2orbitquat(nine[3],nine[2],nine[6],Quatio);

	/*��̬����*/
	QuatTatti_rk4(H,Wib,Quatib);/*��̬�˶�ѧ��������һʱ�̵���̬��Ԫ��Qib  */
	for(i=0;i<3;i++)
		Tb[i]=Twheelb[i]+T_Th[i]+Td_geomagnetic[i]+Td_gravity[i]+Td_aero[i]+Td_Sun[i]+Tmag[i];   /*�������Ǻ�����*/
	//�������ת������[4]������ת��[4]�� ������ֽǶ���[4]
	WheelMomentum(ACS_Jwheel,Wwheel,MomentumWheel);    /*���ֽǶ�������*/
	//������ֽǶ���NMS[4](��֪����ɶ��)�����ְ�װ����[4][3]�� �������ϵͳ�Ƕ���[4]
	WheelSysToThree(MomentumWheel,ACS_Pwheelb,Hwheel);  /*����ϵͳ�Ƕ������� 	*/
	//�������ϵͳ�Ƕ�����������Թ���ϵ�Ľ��ٶ�[3]������ת������[3][3]�� ������ǽǶ���[3]
	formHB(Hwheel,Wib,ACS_IBf,Hbody);                /*���ǽǶ�������*/
	//������ֲ���[1]�����ǽǶ���[3]�����Ǻ�����[3]������ת������[3][3]; �������ϵ��Թ���ϵ�Ľ��ٶ�[3]
	BiasattiDynamic_rk4(H,Hbody,Tb,ACS_IBf,Wib);   /*��̬����ѧ,���ֿ�ϵͳ�ĸ��嶯��ѧ��������һʱ�̵Ľ��ٶ�Wib*/

	/*����ʸ��*/
	//����Ϊ����ʱ��+���ֲ�����̫���ƾ����Ƴཻ��, �������ϵ�µ�̫��ʸ��[3]��̫��ϵ��Թ���ϵ����Ԫ��[4]�����������Ӻ���ʱ��[1]*/
	SunVectori(Time+H,&Yousun,&Aisun,&Jingdu_GreenWitch_inertia,SunInertia,Quatis);  /*̫��ʸ��ģ��*/
	//��������ʱ��+���ֲ����� ����λ���ٶ�[6]�� �������ϵ������ʸ��[3]
	MoonVectori(Time+H,Y,Mooni); /*����ʸ��ģ��*/
	//��������λ���ٶ�[6]�� �������ϵ�µ���ʸ��[3]
	EarthVectori(Y,Earthi); /*����ʸ��ģ��*/
    //��������λ���ٶ�[6]������ϵ��̫��ʸ��[3]�� �����Ӱ��־λ��EshadowFlag=1ʱ�����Ӱ��
	EarthShadow(Y,SunInertia,&EshadowFlag);/*��Ӱ����*/

	/*������*/
    //����̫��ʸ��[3]������ϵ��Թ���ϵ����Ԫ��[4]������̫��1�İ�װ����[3][3]������̫��1���ӳ�������ֵ[1],����̫��1������ϵ��[1]����Ӱ��־λ������̫��1������Ч��־λ
	//�������ֵ�ʹ�����ֵ��̫��ʸ���ڲ���ϵyz��xzƽ���ϵ�ͶӰ�ֱ���z��ļн�
	Fsun(SunInertia,Quatib,ACS_PFsunb1,ACS_Fsun1Limit,ACS_Fsun1R,EshadowFlag,&Fsun1Flag,&Fsun1out0,&Fsun1out);  /*����̫��������1����ģ��*/
	Fsun(SunInertia,Quatib,ACS_PFsunb2,ACS_Fsun2Limit,ACS_Fsun2R,EshadowFlag,&Fsun2Flag,&Fsun2out0,&Fsun2out);  /*����̫��������2����ģ��*/
	//����̫��ʸ��������ʸ��������ʸ��������A����ָ��ֱ���̫������������ʸ���нǵ��ٽ�ֵ[3]����λ�Ƕȡ㣩,����A���������ת��Ԫ��[4]
	//����ϵ��Թ���ϵ����Ԫ��[4]������A����ϵ��[3],����A������Ч��־λ[1]�� �������A/B��װ
	quatstar(SunInertia,Earthi,Mooni,ACS_StarALimit,ACS_StarsbA,Quatib,ACS_quatstarRA,&StarAFlag,QuatsiA0,QuatsiA);    /*��������A����ģ�ͣ����������Ԫ��*/
	quatstar(SunInertia,Earthi,Mooni,ACS_StarBLimit,ACS_StarsbB,Quatib,ACS_quatstarRB,&StarBFlag,QuatsiB0,QuatsiB);    /*��������B����ģ�ͣ����������Ԫ��*/
	Gyro(Wib,ACS_Pgyrob1,ACS_NoGyro1,ACS_GyroR1,ACS_GyroRC1,Wgyro10,Wgyro1);        /*ƽ̨���ݲ���ģ�ͣ�����������ٶ�*/
	Gyro(Wib,ACS_Pgyrob2,ACS_NoGyro2,ACS_GyroR2,ACS_GyroRC2,Wgyro20,Wgyro2);        /*�������ݲ���ģ�ͣ�����������ٶ�*/
	Acceler(F_Th,F_aero,F_Sun,Msat,ACS_Psb,ACS_AccelerR,ACS_AccelerRC,as0,as);   /*���ٶȼƲ���ģ�ͣ���������Ǽ��ٶ�*/
	Magnet(Bci,Quatib,ACS_PbBB,ACS_MagoutR,Magout0,Magout);  /*��ǿ��ģ�ͣ���������ų�ǿ��*/

    Distance_in_Radar(r1,r2, &Distance);
    Angel_in_Camera(r1,r2,ACS_StarsbA,Quatib, &Theta_yz, &Theta_yx);

	/*ϵͳ���*/
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
