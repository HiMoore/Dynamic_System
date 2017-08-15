
#include"math.h"
#include"stdio.h"
#include"stdlib.h"
#include "Initial.h"
#include"Functions.h"


void orbJ2_rk4(double Y[6],double H,double Fb[3],double Quatib[4],double Msat)
/*���������������λ�ú��ٶ�*/
{
	/*�������ϵ������λ���ٶ� Y
	���ֲ��� H
	����ϵ�����Ǻ����� Fb
	����ϵ��Թ���ϵ��̬��Ԫ�� Quatib
	�������� Msat
	�����һʱ�̹���ϵ������λ���ٶ� Y*/
	double DY[6],RK1[6],RK2[6],RK3[6],YSs[6];
	double DT;
	double Quatbi[4],Fi[3],ai[3];
	int I;

	QuatConjuction(Quatib,Quatbi);
	VectorQuatFrameTrans(Fb,Quatbi,Fi);  /*�������ϵ�����Ǻ�����*/
	for (I=0;I<3;I++)
	{
		ai[I]=Fi[I]/Msat;  /*�������ϵ�����Ǽ��ٶ�*/
	}

	DT=0.5*H;
	orbJ2_diff(Y,ai,DY);
	for (I=0;I<=5;I++){
		RK1[I]=DY[I]*DT;
		YSs[I]=Y[I];
		Y[I]=YSs[I]+RK1[I];
	}
	orbJ2_diff(Y,ai,DY);
	for (I=0;I<=5;I++){
		RK2[I]=DY[I]*DT;
		Y[I]=YSs[I]+RK2[I];
	}
	orbJ2_diff(Y,ai,DY);
	for (I=0;I<=5;I++){
		RK3[I]=DY[I]*H;
		Y[I]=YSs[I]+RK3[I];
	}
	orbJ2_diff(Y,ai,DY);
	for (I=0;I<=5;I++){
		Y[I]=YSs[I]+(2.0*(RK1[I]+RK3[I])+4.0*RK2[I]+DY[I]*H)/6.0;
	}
}
void orbJ2_diff(double Y[6],double ai[3],double DY[6])
/*�������ѧ��������*/
{
	/*�������ϵ������λ���ٶ� Y
	����ϵ�����Ǽ��ٶ� ai
	�������ϵ������λ���ٶ����� DY*/
	double R;
	double Cj,Cjz;

	R=sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]); /*���ľ�*/
	Cj=1+1.5*ACD_J2*ACD_Rearth*ACD_Rearth*(1-5*Y[2]*Y[2]/R/R)/R/R;
	Cjz=1+1.5*ACD_J2*ACD_Rearth*ACD_Rearth*(3-5*Y[2]*Y[2]/R/R)/R/R; /*�㶯ϵ��*/
	DY[0]=Y[3];
	DY[1]=Y[4];
	DY[2]=Y[5];
	DY[3]=-ACD_MIUE*Y[0]*Cj/R/R/R+ai[0];
	DY[4]=-ACD_MIUE*Y[1]*Cj/R/R/R+ai[1];
	DY[5]=-ACD_MIUE*Y[2]*Cjz/R/R/R+ai[2];
}

void PV2Six(double Y[6],double nine[9],double *Worbity)
/*��λ�ú��ٶȵ��ƹ����Ҫ�غ͹�����ٶ�*/
/*nine[9]:A,E,SI,DW,SW,DM,SU,theta,n*/
/*�볤��,ƫ����,������,������ྭ,���ĵ�Ǿ�,ƽ�����,γ�ȷ���,������,ƽ��������ٶ�*/
{
	double R[3], V[3];

	double H[3],LengthR,a,esinE,ecosE,e,E,f,i,Omiga,ksinu,kcosu,u,omiga,tmp,M,b,Worbit;

	R[0] = Y[0];R[1] = Y[1];R[2] = Y[2];V[0] = Y[3];V[1] = Y[4];V[2] = Y[5];
	LengthR=LengthVector(R);
	CrossProduct(R,V,H);if (LengthR<ACD_Rearth)LengthR=ACD_Rearth;
	a=1.0/(2.0/LengthR-DotProduct(V,V)/ACD_MIUE);if (a<ACD_Rearth)a=ACD_Rearth;
	ecosE=1.0-LengthR/a;
	esinE=DotProduct(R,V)/sqrt(a*ACD_MIUE);
	e=sqrt(esinE*esinE+ecosE*ecosE);	/*e=LimitUpLow(e,1.0-ACS_EPS,0.0);*/

	if (e>=1)
	{
	    e=1-ACD_EPS;
	    printf("Error : e >= 1");
	}
    else if (e<0)
	{
		e=0;
		printf("Error : e < 0");
	}

	E=atan2(esinE,ecosE);
	f=2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*E));
	tmp=LengthVector(H);	if (tmp<My_abs(H[2]))tmp=My_abs(H[2]);
	i=acos(H[2]/tmp);if (My_abs(i)<ACD_EPS) i = ACD_EPS;
	Omiga=atan2(H[0],-H[1]);
	ksinu=R[2]/sin(i);
	kcosu=R[0]*cos(Omiga)+R[1]*sin(Omiga);
	u=atan2(ksinu,kcosu);
	omiga=u-f;
	M=E-esinE;
	Worbit=sqrt(ACD_MIUE/a/(1-e*e))*(1+e*cos(f))/a/(1-e*e);

	b=a*sqrt(1-e*e);

	nine[0]=a;    /*m*/
	nine[1]=e;
	nine[2]=i;
	nine[3]=Omiga;
	nine[4]=omiga;
	nine[5]=M;
	nine[6]=u;
	nine[7]=f;
	nine[8]=LengthVector(H)/(a*b);   /*ƽ��������ٶ�*/
	*Worbity=-Worbit;

}

void PV2Six_OLD(double Y[6],double nine[9],double *Worbity)
/*��λ�ú��ٶȵ��ƹ����Ҫ�غ͹�����ٶ�*/
/*nine[9]:A,E,SI,DW,SW,DM,SU,theta,n*/
/*�볤��,ƫ����,������,������ྭ,���ĵ�Ǿ�,ƽ�����,γ�ȷ���,������,ƽ��������ٶ�*/
{
	/*�������ϵ������λ���ٶ� Y
	���������� nine
	˲ʱ������ٶ� Worbity*/
	double amu=ACD_MIUE/1.0e9;
	double SF,SU;
	double A,E,SI,DW,SW,DM,b;
	double SR,RV,V,EE1,EE2,EE3;
	double H,HX,HY,HZ,SN,CDE,DE,SDE,Worbit;
	double YY[6];
	double n[3],R[3],c[3];
	unsigned char j;
	for (j=0;j<6;j++)
		YY[j]=Y[j]/1000.0;
	for (j=0;j<3;j++) R[j]=YY[j];
	SR=sqrt(YY[1]*YY[1]+YY[2]*YY[2]+YY[0]*YY[0]);	/*r*/
	V=sqrt(YY[4]*YY[4]+YY[5]*YY[5]+YY[3]*YY[3]);  /*v*/
	RV=YY[1]*YY[4]+YY[2]*YY[5]+YY[0]*YY[3];       /*RV=r.v*/
	/*a---A*/
	A=1/(2/SR-V*V/amu);        				/*got by the energy equation*/
	/*e---E*/
	EE1=((V*V-amu/SR)*YY[0]-RV*YY[3])/amu;
	EE2=((V*V-amu/SR)*YY[1]-RV*YY[4])/amu;
	EE3=((V*V-amu/SR)*YY[2]-RV*YY[5])/amu;
	E=sqrt(EE1*EE1+EE2*EE2+EE3*EE3);

	/*i---SI  (rad)*/
	HX=YY[1]*YY[5]-YY[2]*YY[4];
	HY=YY[2]*YY[3]-YY[0]*YY[5];
	HZ=YY[0]*YY[4]-YY[1]*YY[3];
	H=sqrt(HX*HX+HY*HY+HZ*HZ);        /*H=r*v   the momentum of the satellite*/
	SI=acos(HZ/H);

	SN=sqrt(HX*HX+HY*HY);
	if(SN==0)
	{
		n[0]=1;
		n[1]=0;
		n[2]=0; /*��������ڳ����,�����Ϊ����ϵX��*/
	}
	else
	{
		n[0]=-HY;
		n[1]=HX;
		n[2]=0;
	}
	/*w---SW  (rad)*/  		 /*nine[9]:0A,1E,2SI,3DW,4SW,5DM,6SU,7theta,8n*/

	if(E<1e-7)
		SW=0;
	else
		SW=acos((EE2*n[1]+EE1*n[0])/(SN*E));
	if(EE3<1.0e-10) SW=-SW;

	/*W---DW  ( rad )*/


	/*	if(HX>1.0e-8)  DW=acos(n[0]/SN); *//*+-90 deviation with angle between HY and axis x*/

	/*	else	  	DW=-acos(n[0]/SN);*/
	if(HX>1.0e-8)  DW=acos(n[0]/sqrt(n[0]*n[0]+n[1]*n[1]));

	else	  	DW=-acos(n[0]/sqrt(n[0]*n[0]+n[1]*n[1]));

	/*E---DE*/
	SN=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (E<1e-7)
	{
		CDE=(YY[0]*n[0]+YY[1]*n[1]+YY[2]*n[2])/SR/SN;
		CrossProduct(n,R,c);
		SDE=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])/SR/SN;
	}
	else
	{
		CDE=(1.0-SR/A)/E;
		SDE=(RV/sqrt(amu*A))/E;
	}
	DE=dftan(CDE,SDE);


	/*M---DM  ( rad )*/
	DM=DE-E*SDE;

	/*	while (fabs(DE-E0)>=1.e-15)
	{
	E0=DE;
	DE=DM+E*sin(E0);
	}
	*/
	EE1=sqrt((1+E)/(1-E));
	if (DM<0) DM=DM+ACD_PI*2;

	/*f---SF*/
	SF=2*atan(EE1*tan(DE/2));
	if (DE>ACD_PI) SF=SF+ACD_PI*2;
	/*u---SU*/
	SU=SW+SF;
	/*Worbit*/
	Worbit=sqrt(amu/A/(1-E*E))*(1+E*cos(SF))/A/(1-E*E);

	b=A*sqrt(1-E*E);

	nine[0]=A*1000;    /*m*/
	nine[1]=E;
	nine[2]=SI;
	nine[3]=DW;
	nine[4]=SW;
	nine[5]=DM;
	nine[6]=SU;
	nine[7]=SF;
	nine[8]=H/(A*b);   /*ƽ��������ٶ�*/
	*Worbity=-Worbit;
}

double dftan(double X, double Y)
/*�������н�*/
{
	/*��������ֵ X
	����ֵ Y
	����Ƕ� dftmp*/
	double dftmp;

	dftmp=atan(Y/X);
	if (X<0) dftmp=dftmp+ACD_PI;
	if (X>0) {
		if (Y<=0) dftmp=dftmp+ACD_PI*2;
	}
	return(dftmp);
}
void Iner2orbitquat(double Omega,double i,double u, double Quatio[4])
/* �������ϵ��Թ�������ϵ����̬��Ԫ��*/
{
	/*����������ྭ Omega
	������ i
	γ�ȷ��� u
	������ϵ��Թ���ϵ��̬��Ԫ�� Quatio*/
	double c[4], c1[4], c2[4], q[4], Temp;
	unsigned char j;
	/*Quatio=Qz(Omega)Qx(i-pi/2)Qy(-u-pi/2)*/
	ZQuat(Omega, c1);
	Temp=i-ACD_PI/2.0;
	XQuat(Temp, c2);
	QuatProduct(c1, c2, c);
	Temp=-u-ACD_PI/2.0;
	YQuat(Temp, c1);
	QuatProduct(c, c1, q);
	snormal(q);
	for (j=0;j<4;j++)
		Quatio[j]=q[j];
}

void XQuat(double a, double q[])
{
	q[0] = cos(a / 2.0);
	q[1] = sin(a / 2.0);
	q[2] = q[3] = 0.0;
}


void YQuat(double a, double q[])
{
	q[0] = cos(a / 2.0);
	q[2] = sin(a / 2.0);
	q[1] = q[3] = 0.0;
}


void ZQuat(double a,double q[4])
/*z������̬��Ԫ��*/
{
	q[0] = cos(a / 2.0);
	q[3] = sin(a / 2.0);
	q[1] = q[2] = 0.0;
}

void QuatProduct(double l[], double p[], double q[])
/*��Ԫ���������  q=l*p */
{
	double out[4];
	int j;

	out[0]=l[0]*p[0]-l[1]*p[1]-l[2]*p[2]-l[3]*p[3];
	out[1]=l[0]*p[1]+l[1]*p[0]+l[2]*p[3]-l[3]*p[2];
	out[2]=l[0]*p[2]+l[2]*p[0]+l[3]*p[1]-l[1]*p[3];
	out[3]=l[0]*p[3]+l[3]*p[0]+l[1]*p[2]-l[2]*p[1];
	for(j=0;j<4;j++)q[j]=out[j];
}

void snormal(double Quat[])
/*��Ԫ����һ�� */
{	double d, Nd;

d=Quat[0]*Quat[0]+Quat[1]*Quat[1]+Quat[2]*Quat[2]+Quat[3]*Quat[3];
if (d<ACD_EPS)
Quat[0]=sqrt(fabs(1.0-(Quat[1]*Quat[1]+Quat[2]*Quat[2]+Quat[3]*Quat[3])));
else{
	Nd = sqrt(d);
	Quat[0]=Quat[0]/Nd;
	Quat[1]=Quat[1]/Nd;
	Quat[2]=Quat[2]/Nd;
	Quat[3]=Quat[3]/Nd;
}

if (Quat[0]<0.0){
	Quat[0]=-Quat[0];
	Quat[1]=-Quat[1];
	Quat[2]=-Quat[2];
	Quat[3]=-Quat[3];
}
}
//sfalg--ƫ��������������־��Rsa, Dwa, PabΪ��ƫ������������TbΪ�����أ�IΪת��������HbodyΪ�Ƕ�����WibΪ������ϵ��Թ�������ϵ����̬���ٶ�����ϵ�µķ�����ʽ����
void BiasattiDynamic_rk4(double H,double Hbody[3],double Tb[3],double I[3][3],double Wib[3])
/*��������������Ʊ�������ϵ��Թ�������ϵ�Ľ��ٶȣ��ڱ���ϵ�б�ʾ��*/
{

	/*ATTITUDE DYNAMICS INTEGRATION STEP SIZE = H*/

	double	RK1[3],RK2[3],RK3[3];
	double  YSw[3],DWib[3];
	double  DT;
	int j;

	DT=0.5*H;
	BiasattiDynamic_diff(Tb,I,Hbody,Wib,DWib);

	for (j=0;j<=2;j++){
		RK1[j]=DWib[j]*DT;
		YSw[j]=Wib[j];
		Wib[j]=YSw[j]+RK1[j];
	}

	BiasattiDynamic_diff(Tb,I,Hbody,Wib,DWib);
	for (j=0;j<=2;j++){
		RK2[j]=DWib[j]*DT;
		Wib[j]=YSw[j]+RK2[j];
	}

	BiasattiDynamic_diff(Tb,I,Hbody,Wib,DWib);
	for (j=0;j<=2;j++){
		RK3[j]=DWib[j]*H;
		Wib[j]=YSw[j]+RK3[j];
	}

	BiasattiDynamic_diff(Tb,I,Hbody,Wib,DWib);
	for (j=0;j<=2;j++){
		Wib[j]=YSw[j]+(2.0*(RK1[j]+RK3[j])+4.0*RK2[j]+DWib[j]*H)/6.0;
	}
}

//sfalg--ƫ��������������־��Rsa, Dwa, PabΪ��ƫ������������TΪ�����أ�IΪת��������HbodyΪ�Ƕ�����W��DWΪ���ٶȺͽǼ��ٶȣ�
void BiasattiDynamic_diff(double T[3],double I[3][3],double Hbody[3],double W[3],double DW[3])
/*��̬����ѧ����*/
{
	double I_INV[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	double W_MAT[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	double V1[3]={0,0,0};
	double V2[3]={0,0,0};
	double MAT1[3][3]={{0,0,0},{0,0,0},{0,0,0}};

	InvMatrix(I,I_INV);/*I_INV��I(-1),����I�������*/
	antisymMatrix(W,W_MAT);/*��W��б����W_MAT*/ //б����Ӧ��Ϊ���ԳƷ���
	MatrixProduct(I_INV,W_MAT,MAT1);//������ˣ�I�������W�ķ���/*MAT1=I(-1)*W_MAT*/
	MatrixProductVector(I_INV,T,V1);//�����������I�������T����  /*V1��I(-1)*Tc*/
	MatrixProductVector(MAT1,Hbody,V2);//�����������MAT1�������Hbody�Ƕ��� /* V2��I(-1)*W*Hbody */


	DW[0]=V1[0]-V2[0];/* DW=V1-V2-V3 */
	DW[1]=V1[1]-V2[1];
	DW[2]=V1[2]-V2[2];
}

void MatrixProductVector(double b[3][3],double c[3],double a[3])
/*������������*/
{
	int i,j;

	for(i=0;i<3;i++)
	{
		a[i]=0;
		for(j=0;j<3;j++)
		{
			a[i]=a[i]+b[i][j]*c[j];
		}
	}
}

void InvMatrix(double Mat[3][3],double InvMat[3][3])
{
	/*
	�������棺
	�������double Mat[3][3]
	�������double InvMat[3][3]��InvMat=inv��Mat��
	*/
	unsigned char i,j,k,n,m[3];
	double Mattmp[3][3],InputMat[3][3],y;

	n=3;

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			if (i==j)
			{
				Mattmp[i][j]=1.0;
			}
			else
			{
				Mattmp[i][j]=0.0;
			}
		}
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			InputMat[i][j]=Mat[i][j];
		}
	}

	for (i=0;i<n;i++)
	{
		y=ACD_EPS;
		for (k=0;k<n;k++)
		{
			if (My_abs(InputMat[i][k])>My_abs(y))
			{
				y=InputMat[i][k];
				m[i]=k;
			}
			else{;}
		}
		for (k=0;k<n;k++)
		{
			InputMat[i][k]=InputMat[i][k]/y;
			Mattmp[i][k]=Mattmp[i][k]/y;
		}
		for (j=0;j<n;j++)
		{
			if (j!=i)
			{
				y=InputMat[j][m[i]];
				for (k=0;k<n;k++)
				{
					InputMat[j][k]=InputMat[j][k]-InputMat[i][k]*y;
					Mattmp[j][k]=Mattmp[j][k]-Mattmp[i][k]*y;
				}
			}
			else{;}
		}
	}
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			InvMat[m[i]][j]=Mattmp[i][j];
		}
	}
}

void MatrixProduct(double b[3][3],double c[3][3],double a[3][3])
/*�������*/
{
	int i,j,k;

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			a[i][j]=0;
			for(k=0;k<3;k++)
			{
				a[i][j]=b[i][k]*c[k][j]+a[i][j];
			}
		}
	}
}

void antisymMatrix(double a[3],double a_mat[3][3])
/*��б����*/
{
	a_mat[0][0]=0.0;
	a_mat[0][1]=-a[2];
	a_mat[0][2]=a[1];

	a_mat[1][0]=a[2];
	a_mat[1][1]=0.0;
	a_mat[1][2]=-a[0];

	a_mat[2][0]=-a[1];
	a_mat[2][1]=a[0];
	a_mat[2][2]=0.0;
}

void QuatTatti_rk4(double H,double Wib[3],double Q[4])
/*��������������Ʊ�������ϵ��Թ������ϵ����̬��Ԫ��*/
{

	/*ATTITUDE DYNAMICS INTEGRATION STEP SIZE = H*/

	double	RK1[4],RK2[4],RK3[4];
	double  YS[4],DQ[4];
	double  DT;
	int i;


	DT=0.5*H;
	QuatTatti_diff(Q,Wib,DQ);

	for (i=0;i<=3;i++){
		RK1[i]=DQ[i]*DT;
		YS[i]=Q[i];
		Q[i]=YS[i]+RK1[i];
	}

	QuatTatti_diff(Q,Wib,DQ);
	for (i=0;i<=3;i++){
		RK2[i]=DQ[i]*DT;
		Q[i]=YS[i]+RK2[i];
	}

	QuatTatti_diff(Q,Wib,DQ);
	for (i=0;i<=3;i++){
		RK3[i]=DQ[i]*H;
		Q[i]=YS[i]+RK3[i];
	}

	QuatTatti_diff(Q,Wib,DQ);
	for (i=0;i<=3;i++){
		Q[i]=YS[i]+(2.0*(RK1[i]+RK3[i])+4.0*RK2[i]+DQ[i]*H)/6.0;
	}
	snormal(Q);
}

void QuatTatti_diff(double Q[4],double Wib[3],double DQ[4])
/*��̬�˶�ѧ����*/
{

	double W_Q[4],Q_temp[4];
	int I;

	W_Q[0]=0;
	W_Q[1]=Wib[0];
	W_Q[2]=Wib[1];
	W_Q[3]=Wib[2];

	QuatProduct(Q,W_Q,Q_temp);

	for (I=0;I<=3;I++)
	{
		DQ[I]=0.5*Q_temp[I];
	}
}

void Geomagnetic(double Y[6],double Jingdu_GreenWitch_inertia,double gnm[7][7],double hnm[7][7],unsigned char Number_mag,double Bci[3])
{
	/*�ڱ���������ϵ�м���ų���Ȼ��ת����J2000����ϵ�� */
	/*
	Geomagnetic intensity compute
	Input:  S1_Position,Satellite Orbit Parameter in inertia frame J2000,[x,y,z],m
	gnm:Guass Parameter;
	hnm:Guass Parameter;
	Numbermag,Compute order

	Ridus -- position distance to earth center [m]
	Weidu -- position earth center Latitude [rad]
	Jingdu -- position Longitude [rad]
	Number -- computer order
	gnm, hnm -- compile Guass parameter
	Bneg -- out result in north-east-gravity coordenate frame
	Output: Bci -- out result in initial coordenate frame. T */
	int i;
	double S1_Position[3];
	unsigned char m, n, k,Numbermag;
	double adr, adrn2, Knm, Pnm[7][7], DPnm[7][7], BBr, BBe, BBn ;
	static double Bneg[3], Ridus=1.0;
	double cosml[7], sinml[7], cosJingdu, sinJingdu, costheta, sintheta, cosChijing, sinChijing;
	double sinJingdu_GreenWitch_inertia,cosJingdu_GreenWitch_inertia;
	/* s00=1;sn0=s(n-1)0*(2n-1)/n,snm=sn(m-1)sqrt((n-m+1)(km+1)/(n+m)) */
	static	double Snm[7][7] = {
		{ 1.0,          0.0,          0.0,         0.0,         0.0,         0.0,         0.0},
		{ 1.0000000000, 1.0000000000, 0.0,         0.0,         0.0,         0.0,         0.0},
		{ 1.5000000000, 1.7320508076, 0.8660254038,0.0,         0.0,         0.0,         0.0},
		{ 2.5000000000, 3.0618621785, 1.9364916731,0.7905694150,0.0,         0.0,         0.0},
		{ 4.3750000000, 5.5339859053, 3.9131189606,2.0916500663,0.7395099729,0.0,         0.0},
		{ 7.8750000000,10.1665812838, 7.6852130745,4.7062126493,2.2185299187,0.7015607600,0.0},
		{14.4375000000,18.9031247417,14.9442322695,9.9628215130,5.4568620791,2.3268138086,0.6716932894}};
		/*static double Jingdu0=-0.52694514152154; 1999082200 */

		/*
		�μ� ����ֹ����ֲᣬ page 247
		2000-01-01-00h Greenwitch inertia star time angle 99.968deg,
		GPS time is 11404813s=132day+13s ���	1999.08.22.00h
		Jingdu0=99.968*ACD_d2r-0.729211585e-4*11404813+132.0*2*pi;*/
		for(i=0;i<3;i++)
		{
			S1_Position[i]=Y[i];
		}

		/* �������ǵ��ľ� */
		Ridus=LengthVector(S1_Position);
		if (Ridus<ACD_Rearth)
		{
			Ridus=ACD_Rearth;
		}
		else{;}

		/* �������µ����γ�Ⱥ;��� */
		/*
		Weidu=asin(S1_Position[2]/Ridus);
		if((My_abs(S1_Position[0])<ACD_EPS)&&(My_abs(S1_Position[1])<ACD_EPS))S1_Position[0]=1.0;
		Jingdu=atan2(S1_Position[1],S1_Position[0])-Jingdu_GreenWitch_inertia;
		cosJingdu = cos(Jingdu);
		sinJingdu = sin(Jingdu);

		Chijing=Jingdu+Jingdu_GreenWitch_inertia;
		cosChijing=cos(Chijing);
		sinChijing=sin(Chijing);

		theta=ACD_PI/2.0-Weidu;
		costheta=cos(theta);
		sintheta=sin(theta);*/
		/* �������ǳྭ����γ */
		costheta=S1_Position[2]/Ridus;
		sintheta=sqrt(1.0-costheta*costheta);

		if (sintheta<0.001)
		{
			sintheta=0.001;
		}
		else{;}
		cosChijing=S1_Position[0]/Ridus/sintheta;
		sinChijing=S1_Position[1]/Ridus/sintheta;

		sinJingdu_GreenWitch_inertia=sin(Jingdu_GreenWitch_inertia);
		cosJingdu_GreenWitch_inertia=cos(Jingdu_GreenWitch_inertia);
		cosJingdu = cosChijing*cosJingdu_GreenWitch_inertia+sinChijing*sinJingdu_GreenWitch_inertia;
		sinJingdu = sinChijing*cosJingdu_GreenWitch_inertia-cosChijing*sinJingdu_GreenWitch_inertia;

		/* Lengend����ʽ��ʼ�� */
		Pnm[0][0]=1.0;
		DPnm[0][0]=0.0;
		adr=ACD_Rearth/Ridus;
		adrn2=adr*adr;

		Numbermag=Number_mag;
		if (Numbermag>6)
		{
			Numbermag=6;
		}
		else{;}

		/* ���㾭�ȵ������� */
		cosml[0]=1.0;
		sinml[0]=0.0;
		for (k=1;k<=Numbermag;k++)
		{
			cosml[k]=cosml[k-1]*cosJingdu-sinml[k-1]*sinJingdu;
			sinml[k]=sinml[k-1]*cosJingdu+cosml[k-1]*sinJingdu;
		}

		/* ��ʼ�������شų����� */
		Bneg[0]=Bneg[1]=Bneg[2]=0.0;

		for(n=1;n<=Numbermag;n++)
		{
			/* Lengend����ʽ��ʼ������ */
			DPnm[n][n]=sintheta*DPnm[n-1][n-1]+costheta*Pnm[n-1][n-1];
			Pnm[n][n]=sintheta*Pnm[n-1][n-1];

			BBr=BBe=BBn=0.0;
			for(m=0;m<=n;m++)
			{
				if(n==1)
				{
					Knm=0.0;
				}
				else
				{
					Knm=((n-1.0)*(n-1.0)-m*m)/(2.0*n-1.0)/(2.0*n-3.0);
				}

				if(m==(n-1))
				{
					Pnm[n][m]=costheta* Pnm[n-1][m];
					DPnm[n][m]=costheta*DPnm[n-1][m]-sintheta*Pnm[n-1][m];
				}
				else if (m<(n-1))
				{
					Pnm[n][m]=costheta* Pnm[n-1][m]-Knm*Pnm[n-2][m];
					DPnm[n][m]=costheta*DPnm[n-1][m]-sintheta*Pnm[n-1][m]-Knm*DPnm[n-2][m];
				}
				else{;}
				BBn=BBn+  ( gnm[m][n]*cosml[m]+hnm[m][n]*sinml[m])*Snm[n][m]*DPnm[n][m];
				BBe=BBe+m*(-gnm[m][n]*sinml[m]+hnm[m][n]*cosml[m])*Snm[n][m]* Pnm[n][m];
				BBr=BBr+  ( gnm[m][n]*cosml[m]+hnm[m][n]*sinml[m])*Snm[n][m]* Pnm[n][m];
			}
			adrn2=adrn2*adr;
			Bneg[0]=Bneg[0]+adrn2*BBn;
			Bneg[1]=Bneg[1]-adrn2*BBe;
			Bneg[2]=Bneg[2]-adrn2*(n+1)*BBr;
		}
		Bneg[1]=Bneg[1]/sintheta;

		/* ��λ�任: nT->T */
		for (k=0;k<3;k++)
		{
			Bneg[k]=Bneg[k]/1.0e9;
		}
		/* ����ϵ�任 Bci=Az(-Chijing)Ay(180-costheta)Bneg */
		Bci[0]=-cosChijing*costheta*Bneg[0]-sinChijing*Bneg[1]-cosChijing*sintheta*Bneg[2];
		Bci[1]=-sinChijing*costheta*Bneg[0]+cosChijing*Bneg[1]-sinChijing*sintheta*Bneg[2];
		Bci[2]=sintheta*Bneg[0]-costheta*Bneg[2];
}

double LengthVector(double Vector[])
{
	/*
	������ģnorm or length of a vector  ���룺
	������Ԫ��double Vector[]
	*/
	double in,out;

	in=Vector[0]*Vector[0]+Vector[1]*Vector[1]+Vector[2]*Vector[2];
	out=sqrt(in);
	return(out);
	/*
	return(sqrt(Vector[0]*Vector[0]+Vector[1]*Vector[1]+Vector[2]*Vector[2]));*/
}

void GeomagneticTorque(double Bci[3],double Quatib[4],double Mdis[3],double Td_geomagnetic[3])
{
	double Bcb[3];

	VectorQuatFrameTrans(Bci,Quatib,Bcb);
	CrossProduct(Mdis,Bcb,Td_geomagnetic);
}

void CrossProduct(double Vin1[],double Vin2[],double Vout[])
/* ������� */
{
	unsigned char j;
	double out[3];

	out[0]=Vin1[1]*Vin2[2]-Vin1[2]*Vin2[1];
	out[1]=Vin1[2]*Vin2[0]-Vin1[0]*Vin2[2];
	out[2]=Vin1[0]*Vin2[1]-Vin1[1]*Vin2[0];
	for (j=0;j<3;j++)
		Vout[j]=out[j];
}

void VectorQuatFrameTrans(double Rin[],double QuatRin2Rout[],double Rout[])
{
	/*
	��Ԫ������ϵ�任  Rout=QuatRin2Rout'*Rin*QuatRin2Rout
	����ʸ���ķ���,double Rin[3]
	������̬��Ԫ�� double QuatRin2Rout[4]
	�����Ӧʸ������ double Rout[3]
	*/
	unsigned char j;
	double QuatRout2Rin[4],QTEMP1[4],QTEMP2[4];

	QuatRout2Rin[0]=QuatRin2Rout[0];
	for (j=1;j<4; j ++)
	{
		QuatRout2Rin[j]=-QuatRin2Rout[j];
	}

	QTEMP1[0]=0.0;
	for (j=1;j<4; j ++)
	{
		QTEMP1[j]=Rin[j-1];
	}
	QuatProduct(QuatRout2Rin,QTEMP1,QTEMP2);
	QuatProduct(QTEMP2,QuatRin2Rout,QTEMP1);
	for (j=0;j<3; j ++)
	{
		Rout[j]=QTEMP1[j+1];
	}
}

void SunVectori(double Time,double *Yousun,double *Aisun,double *Jingdu_GreenWitch_inertia,double SunInertia[],double Quatis[])
{
	/*
	̫�����Ĳ�������
	Input:  time,time from 1980010600,s;GPSTime=0,Initial Time
	Input:  time,time from 1999082200,s;GPSTime+1024*7+13=0
	Input:  time,time from 2004010100,s;GPSTime+1024*7+365.25*4+132+13=0
	J2000:  2000-01-01-12h;2000010112,s,JD2000=GPSTime+1024*7+132.5+13=0
	Output: Sun Ecliptic, Rad*/
	double Esun,Omigasun,ais,cosai,sinai,sinEcsuna;
	double Msun, Ecsuna,JD2000,c1[4],c2[4],Jingdu_GreenWitch_inertia0;
	static double Time0=0.0;/* Beijing Time 2004010100.00->UTC(s): (365.25*4-0.5)*86400.0-8*3600.0 */

	/*	JD2000=(Time-13.0-(132.5+1024.0*7)*86400.0)/86400.0/36525.0;	��ʼʱ��s need update !*/
	/*  JD2000=(Time+(365.25*4-0.5)*86400.0)/86400.0/36525.0;	 ��ʼʱ��s need update !*/
	JD2000=(Time+Time0)/86400.0/36525.0;	/* ��ʼʱ��s need update !*/

	/* Page161 Book ���������ȷ��*/
	/* ����ƽ����ƽ����� */
	Msun=357+31.0/60+44.76/3600+129596581.04/3600*JD2000-0.562/3600*JD2000*JD2000;
	Msun=Msun*ACD_d2r;

	/* ����ƽ����ƫ���� */
	Esun=0.01670862-0.00004204*JD2000-0.00000124*JD2000*JD2000;

	/* ����ƽ�����Ƶ���ǣ��Ƴཻ�ǣ� */
	ais=23.0+26.0/60+21.448/3600-46.8150/3600*JD2000-0.00059/3600*JD2000*JD2000;
	ais=ais*ACD_d2r;

	/* ����ƽ�������ص���ǣ�ƽ�ƾ��� */
	Omigasun=282.0+56.0/60+14.45/3600+6190.32/3600*JD2000+1.655/3600*JD2000*JD2000;
	Omigasun=Omigasun*ACD_d2r;

	/* ����ƽ����̫��γ�ȷ��ǣ�̫���ƾ��� */
	/*	Ecsuna = Msun+2.0*Esun*sin(Msun)+1.25*Esun*Esun*sin(2.0*Msun)+Omigasun ;*/
	Ecsuna=M2f(Msun,Esun)+Omigasun;

	Ecsuna=AngleWrap(Ecsuna);

	/* ����̫����J2000����ϵ�е�λ������ */
	cosai=cos(ais);
	sinai=sin(ais);
	sinEcsuna=sin(Ecsuna);
	SunInertia[0]=cos(Ecsuna);
	SunInertia[1]=sinEcsuna*cosai;
	SunInertia[2]=sinEcsuna*sinai;

	/* ����̫���������ϵ��J2000����ϵ�е���̬��Ԫ�� */
	XQuat(ais-ACD_PI/2.0, c1);
	YQuat(-Ecsuna-ACD_PI/2.0, c2);
	QuatProduct(c1, c2, Quatis);
	snormal(Quatis);

	*Yousun=Ecsuna;
	*Aisun=ais;
	/* ������������Ӻ���ʱ�� Page106 Book ����̽����������
	*Jingdu_GreenWitch_inertia=2.0*ACD_PI*(67310.54841/86400.0+(876600.0/24.0+8640184.812866/86400.0)*JD2000
	+0.093104/86400.0*JD2000*JD2000-6.2e-6/86400.0*JD2000*JD2000*JD2000);
	������������Ӻ���ʱ�� Page11 Book ������������� */
	Jingdu_GreenWitch_inertia0=(18.6973746+879000.0513367*JD2000
		+0.093104/3600.0*JD2000*JD2000
		-6.2e-6/3600.0*JD2000*JD2000*JD2000)*15.0*ACD_d2r;

	*Jingdu_GreenWitch_inertia=AngleWrap(Jingdu_GreenWitch_inertia0);
}

double M2f(double M, double e)
/* ƽ����ǵ�������ת�� */
{
	double f;
	f=M+(2.0*e-e*e*e/4.0)*sin(M)
		+(1.25*e*e-11.0/24.0*e*e*e*e)*sin(2.0*M)
		+(13.0/12.0*e*e*e)*sin(3.0*M)
		+(103.0/96.0*e*e*e*e)*sin(4.0*M);
	return(f);
}

double AngleWrap(double u)
/*��������pi����pi֮��*/
{
	double out;

	out = u-2.0*ACD_PI*(int)(u/(2.0*ACD_PI));
	if (out>=ACD_PI)
		out-=2.0*ACD_PI;
	else if (u<-ACD_PI)
		out+=2.0*ACD_PI;
	else{}

	return(out);
}

void quatstar(double SunInertia[3],double Earthi[3],double Mooni[3],double StarLimit[3],double Starsb[4],double Quatib[4],double quatstarR[3],unsigned char *StarFlag,double Quatsi0[4],double Quatis[4])
{
	/* ���������������������С */
	double Quatbi[4],Starbs[4];
	double rx,ry,rz;
	double Qx[4],Qy[4],Qz[4],Qr[4];
	double Qtmp[4];
	double Sunb[3],Earthb[3],Moonb[3];
	double Suns[3],Earths[3],Moons[3];
	double aSun,aEarth,aMoon;
	double StarZ[3]={0,0,1};/*������������*/
	double Quatsi[4];
	static int tS=0;
	QuatConjuction(Quatib,Quatbi);
	QuatProduct(Starsb,Quatbi,Quatsi0);
	snormal(Quatsi0);/*����ֵ*/
	srand((unsigned int)time(0)+tS++);
	rx=((double)(rand()*2.0/RAND_MAX-1))*quatstarR[0];
	ry=((double)(rand()*2.0/RAND_MAX-1))*quatstarR[1];
	rz=((double)(rand()*2.0/RAND_MAX-1))*quatstarR[2];
	XQuat(rx,Qx);
	YQuat(ry,Qy);
	ZQuat(rz,Qz);
	QuatProduct(Qx,Qy,Qtmp);
	QuatProduct(Qtmp,Qz,Qr);/*��������*/
	QuatProduct(Quatsi0,Qr,Quatsi);
	snormal(Quatsi);/*�������Ĳ������*/
	/*������������Ч��*/
	QuatConjuction(Starsb,Starbs);
	VectorQuatFrameTrans(SunInertia,Quatib,Sunb);
	VectorQuatFrameTrans(Sunb,Starbs,Suns);
	aSun=acos(DotProduct(Suns,StarZ)/LengthVector(Suns));
	*StarFlag=0;
	if(aSun<StarLimit[0]*ACD_PI/180)
	{
		*StarFlag=*StarFlag+2;/*������̫��ʸ���н�С���ٽ�ֵ,�����Ч*/
	}
	VectorQuatFrameTrans(Earthi,Quatib,Earthb);
	VectorQuatFrameTrans(Earthb,Starbs,Earths);
	aEarth=acos(DotProduct(Earths,StarZ)/LengthVector(Earths));
	if(aEarth<StarLimit[1]*ACD_PI/180)
	{
		*StarFlag=*StarFlag+4;/*���������ʸ���н�С���ٽ�ֵ,�����Ч*/
	}
	VectorQuatFrameTrans(Mooni,Quatib,Moonb);
	VectorQuatFrameTrans(Moonb,Starbs,Moons);
	aMoon=acos(DotProduct(Moons,StarZ)/LengthVector(Moons));
	if(aMoon<StarLimit[2]*ACD_PI/180)
	{
		*StarFlag=*StarFlag+8;/*����������ʸ���н�С���ٽ�ֵ,�����Ч*/
	}
	if(*StarFlag==0)
	{
		*StarFlag=1;/*�����Ч*/
	}
	QuatConjuction(Quatsi,Quatis);
}

void QuatConjuction(double Quat[],double QuatConj[])
/* ������Ԫ�� */
{
	QuatConj[0]=Quat[0];
	QuatConj[1]=-Quat[1];
	QuatConj[2]=-Quat[2];
	QuatConj[3]=-Quat[3];
}

void Quat2Tatti(double Quat[], double Tatti[3][3])
{
	/*
	��Ԫ�����������Ҿ���ת��  Tatti=(Quat0-Quat*)^2+Quat*Quat'
	������Ԫ��,double Quat[4]
	�����Ӧ��̬���� double Tatti[3][3],
	*/
	/*
	attitude Quaternion to attitude matrix
	Input:  Quat,Input Quaternion
	Output: Tatti,Attitude Matrix*/
	Tatti[0][0] = 1.0 - 2.0 * ( Quat[2]*Quat[2] + Quat[3]*Quat[3]);
	Tatti[0][1] =       2.0 * ( Quat[1]*Quat[2] + Quat[0]*Quat[3]);
	Tatti[0][2] =       2.0 * (-Quat[0]*Quat[2] + Quat[1]*Quat[3]);
	Tatti[1][0] =       2.0 * ( Quat[1]*Quat[2] - Quat[0]*Quat[3]);
	Tatti[1][1] = 1.0 - 2.0 * ( Quat[1]*Quat[1] + Quat[3]*Quat[3]);
	Tatti[1][2] =       2.0 * ( Quat[0]*Quat[1] + Quat[2]*Quat[3]);
	Tatti[2][0] =       2.0 * ( Quat[0]*Quat[2] + Quat[1]*Quat[3]);
	Tatti[2][1] =       2.0 * (-Quat[0]*Quat[1] + Quat[2]*Quat[3]);
	Tatti[2][2] = 1.0 - 2.0 * ( Quat[1]*Quat[1] + Quat[2]*Quat[2]);
}


void Tatti2Quat(double Tatti[3][3], double Quat[4])
{
	/*
	�������Ҿ�����Ԫ��ת��
	������̬����,double Tatti[3][3]
	�����Ӧ��Ԫ�� double Quat[4],
	*/
	/* attitude matrix to attitude Quaternion */
	double Trace;

	Trace=Tatti[0][0]+Tatti[1][1]+Tatti[2][2];
	if (Trace>=0.0)
	{
		Quat[0]=sqrt(1.0+Trace)/2.0;
		Quat[1]=(Tatti[1][2]-Tatti[2][1])/Quat[0]/4.0;
		Quat[2]=(Tatti[2][0]-Tatti[0][2])/Quat[0]/4.0;
		Quat[3]=(Tatti[0][1]-Tatti[1][0])/Quat[0]/4.0;
	}
	else
	{
		if ((Tatti[1][1]>Tatti[0][0])&&(Tatti[1][1]>Tatti[2][2]))
		{
			Quat[2]=sqrt(My_abs(1.0-Tatti[0][0]+Tatti[1][1]-Tatti[2][2]))/2.0;

			if (Quat[2]<ACD_EPS)
			{
				Quat[2]=ACD_EPS;
			}
			else{;}
			Quat[0]=(Tatti[2][0]-Tatti[0][2])/Quat[2]/4.0;
			Quat[1]=(Tatti[1][0]+Tatti[0][1])/Quat[2]/4.0;
			Quat[3]=(Tatti[2][1]+Tatti[1][2])/Quat[2]/4.0;
		}
		else if (Tatti[2][2]>Tatti[0][0])
		{
			Quat[3]=sqrt(My_abs(1.0-Tatti[0][0]-Tatti[1][1]+Tatti[2][2]))/2.0;

			if (Quat[3]<ACD_EPS)
			{
				Quat[3]=ACD_EPS;
			}
			else{;}
			Quat[0]=(Tatti[0][1]-Tatti[1][0])/Quat[3]/4.0;
			Quat[1]=(Tatti[2][0]+Tatti[0][2])/Quat[3]/4.0;
			Quat[2]=(Tatti[2][1]+Tatti[1][2])/Quat[3]/4.0;
		}
		else
		{
			Quat[1]=sqrt(My_abs(1.0+Tatti[0][0]-Tatti[1][1]-Tatti[2][2]))/2.0;

			if (Quat[1]<ACD_EPS)
			{
				Quat[1]=ACD_EPS;
			}
			else{;}
			Quat[0]=(Tatti[1][2]-Tatti[2][1])/Quat[1]/4.0;
			Quat[2]=(Tatti[1][0]+Tatti[0][1])/Quat[1]/4.0;
			Quat[3]=(Tatti[2][0]+Tatti[0][2])/Quat[1]/4.0;
		}
	}
	snormal(Quat);
}

double My_abs(double x)
{
	/*
	����ֵ���  ���룺
	�������double x
	*/
	if (x<0.0)
	{
		return(-x);
	}
	else
	{
		return(x);
	}
}

void WheelSysToThree(double WheelSys[4],double Pwheelb[4][3],double WheelThree[3])
{
	/*����ϵͳ�Ƕ������� */
	/*
	wheel Control torque Limitup
	Input:  WheelSys[4],���ֽǶ���NMS
	Pwheelb[4][3],���ְ�װ����
	Output: WheelThree[3],����ϵͳ�Ƕ���*/
	int j,k;
	for (j=0;j<3;j++)
	{
		WheelThree[j]=0.0;
		for (k=0;k<4;k++)
		{
			WheelThree[j]=WheelThree[j]+Pwheelb[k][j]*WheelSys[k];
		}
	}
}

void formHB(double Hwheel[3],double Wbody[3],double IBf[3][3],double Hbody[3])
{
	/*
	���ǽǶ�������
	Input:  Hwheel[],Wheel monmentum, Nms;
	Wbody[]=[Wbx Wby Wbz],Satellite error anglar Velocity,rad/s;
	IBf[][], Inertial matrix,kgmm
	Output: Hbody,Angelar monmentun of system,Nms*/
	unsigned char j,k;
	double out[3];

	for(j=0;j<3;j++)
	{
		out[j]=Hwheel[j];
		for(k=0;k<3;k++)
		{
			out[j]=out[j]+IBf[j][k]*Wbody[k];
		}
	}
	for(j=0;j<3;j++)
	{
		Hbody[j]=out[j];
	}
}

int My_Sign(double x)
{
	if (x < 0)	return(-1);
	else 		return(1);
}

void WheelMomentum(double Jwheel[4],double Wwheel[4],double MomentumWheel[4])
{
	/*������ֽǶ���*/
	int i;
	for (i=0;i<4;i++)
	{
		MomentumWheel[i]=Jwheel[i]*Wwheel[i];
	}
}

void Gyro(double w[3],double Pgyrob[3][3],unsigned char NoGyro[3],double GyroR[3],double GyroRC[3],double Wgyro0[3],double Wgyro[3])
{
	/* ���ݲ������ */
	double we[3];
	static int tG=0;
	int k,j;

	for (j=0;j<3;j++)
	{
		we[j] = 0;
		for (k=0;k<3;k++)
		{
			/* ��ȡ�µ����ݷ������� */
			we[j] = we[j]+Pgyrob[j][k]*w[k];
		}
	}

	for (j=0;j<3;j++)
	{
		Wgyro0[j]=0;
		Wgyro[j]=0;
	}
	srand((unsigned int)time(0)+tG++);
	for (j=0;j<3;j++)
	{
		Wgyro0[j]=we[j];
		Wgyro[j]=we[j]+GyroRC[j]+((double)(rand()*2.0/RAND_MAX-1))*GyroR[j];
	}
}

//��ʸ��1��aϵ�任��2ϵ
void VectorCosMatrixFrameTrans(double Rin[],double CosMatrixOutIn[3][3],double Rout[])
/*������������ϵ�任  Rout=CosMatrixOutIn*Rin */
{
	unsigned char j,k;
	double out[3];

	for (k=0;k<3; k++) {
		out[k]=0.0;
		for (j=0;j<3;j++)
			out[k]=out[k]+CosMatrixOutIn[k][j]*Rin[j];
	}
	for (j=0;j<3;j++)
		Rout[j]=out[j];
}

void Fsun(double SunInertia[3],double Quatib[4],double PFsunb[3][3],double FsunLimit,double FsunR,unsigned char EshadowFlag,unsigned char *FsunFlag,double *Fsunout0,double *Fsunout)
{
	/*̫������������ģ��*/
	/*����̫��ʸ��
	���̫��ʸ���ڲ�������ϵyz��xzƽ���ϵ�ͶӰ�ֱ���z��ļн�*/
	double Sunb[3],Suns[3];
	double S[3];/*̫����������װ����,����װ����ϵZ��*/
	/*	double PFsunb[3][3];*/
	double a;
	static int tF=0;
	/*	InvMatrix(PbFsun,PFsunb);      */              /*��������̫���������������Ҿ���*/
	if(EshadowFlag==1)
	{
		*FsunFlag=0;/*�����Ӱ,�����������Ч*/
		*Fsunout0=0;/* ����ֵ*/
		*Fsunout=0;/*���������*/
	}
	else
	{
		VectorQuatFrameTrans(SunInertia,Quatib,Sunb);  /*���㱾������ϵ�µ�̫��ʸ������*/
		S[0]=PFsunb[2][0];
		S[1]=PFsunb[2][1];
		S[2]=PFsunb[2][2];
		a=DotProduct(Sunb,S)/LengthVector(Sunb);
		if(a<=0)
		{
			*FsunFlag=0;/*̫��ʸ����������������,�����������Ч*/
			*Fsunout0=0;/* ����ֵ*/
			*Fsunout=0;/*���������*/
		}
		else
		{
			VectorCosMatrixFrameTrans(Sunb,PFsunb,Suns);  /*�����������ϵ�µ�̫��ʸ������*/
			*Fsunout0=atan(Suns[1]/Suns[2]);/* ����ֵ*/
				srand((unsigned)time(0)+tF++);
			*Fsunout=*Fsunout0+((double)(rand()*2.0/RAND_MAX-1))*FsunR;/*���������*/
			if(My_abs(*Fsunout)>FsunLimit)
			{
				*FsunFlag=0;/*�����ӳ���Χ,�����Ч*/
			}
			else
			{
				*FsunFlag=1;
			}
		}
	}
}

/*void Fsun(double SunInertia[3],double Quatib[4],double PbFsun[3][3],double FsunLimit[2],double FsunR[2],unsigned char EshadowFlag,unsigned char *FsunFlag,double Fsunout0[2],double Fsunout[2])
{*/
/*̫������������ģ��*/
/*����̫��ʸ��
���̫��ʸ���ڲ�������ϵyz��xzƽ���ϵ�ͶӰ�ֱ���z��ļн�*/
/*double Sunb[3],Suns[3];
double PFsunb[3][3];
static int tF=0;
if(EshadowFlag==1)
{
*FsunFlag=0;*//*�����Ӱ,�����������Ч*/
/*Fsunout0[0]=0;
Fsunout0[1]=0; ����ֵ*/
/*Fsunout[0]=0;
Fsunout[1]=0;���������*/
/*}
else
{
VectorQuatFrameTrans(SunInertia,Quatib,Sunb);  ���㱾������ϵ�µ�̫��ʸ������*/
/*InvMatrix(PbFsun,PFsunb);                    ��������̫��������1�������Ҿ���*/
/*VectorCosMatrixFrameTrans(Sunb,PFsunb,Suns);  �����������ϵ�µ�̫��ʸ������*/
/*Fsunout0[0]=atan(Suns[0]/Suns[2]);
Fsunout0[1]=atan(Suns[1]/Suns[2]); ����ֵ*/
/*srand((unsigned)time(0)+tF++);
Fsunout[0]=Fsunout0[0]+((double)(rand()%3)-1)*FsunR[0];
Fsunout[1]=Fsunout0[1]+((double)(rand()%3)-1)*FsunR[1];���������*/
/*if((My_abs(Fsunout[0])>FsunLimit[0])||(My_abs(Fsunout[1])>FsunLimit[1]))
{
*FsunFlag=0;�����ӳ���Χ,�����Ч*/
/*  }
else
{
*FsunFlag=1;
}
}
}*/

void Sun01(double SunInertia[3],double Quatib[4],double Psun01b[3][3],unsigned char EshadowFlag,unsigned char Sun01out[5])
{
	/*0-1ʽ̫������������ģ��*/
	/*����̫��ʸ��
	���̫������������0/1�ź�*/
	double Sunb[3],Suns[3];
	/*	double Psun01b[3][3];*/
	double Zs[3]={0,0,1};
	double Xs[3]={1,0,0};
	double SunXYs[3];
	double Zangle,Xangle;
	double SL;
	int i;

	VectorQuatFrameTrans(SunInertia,Quatib,Sunb);  /*���㱾������ϵ�µ�̫��ʸ������*/
	/*	InvMatrix(Pbsun01,Psun01b);     */               /*��������̫��������1�������Ҿ���*/
	VectorCosMatrixFrameTrans(Sunb,Psun01b,Suns);  /*�����������ϵ�µ�̫��ʸ������*/
	Zangle=acos(DotProduct(Suns,Zs)/LengthVector(Suns));  /*����̫��ʸ�����������ϵZ��ļн�*/
	SunXYs[0]=Suns[0];
	SunXYs[1]=Suns[1];
	SunXYs[2]=0;
	SL=LengthVector(SunXYs);
	if (SL==0)
	{
		Xangle=0;
		printf("sun01 : Length_SunXYs = 0");
	}
	else
	    Xangle=acos(DotProduct(SunXYs,Xs)/LengthVector(SunXYs));  /*����̫��ʸ���ڲ�������ϵXYƽ���ͶӰ��X��ļн�*/
	for (i=0;i<5;i++)
	{
		Sun01out[i]=0;
	}
	if ((Suns[2]>0)&&(EshadowFlag==0))
	{
		if (Zangle<ACD_PI*45.0/180.0)
		{
			Sun01out[4]=1;
		}
		if ((Zangle>ACD_PI*20.0/180)&&(Zangle<ACD_PI/2.0))
		{
			if (((Suns[1]>0)&&(Xangle<ACD_PI*95.0/180.0))||((Suns[1]<0)&&(Xangle<ACD_PI*5.0/180.0))||((Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*105.0/180.0))||((Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*15.0/180.0)))
			{
				Sun01out[0]=1;
			}
			if (((Suns[1]>0)&&(Xangle>ACD_PI*85.0/180.0))||((Suns[1]<0)&&(Xangle>ACD_PI*175.0/180.0))||((Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*75.0/180.0))||((Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*165.0/180.0)))
			{
				Sun01out[3]=1;
			}
			if (((Suns[1]>0)&&(Xangle>ACD_PI*175.0/180.0))||((Suns[1]<0)&&(Xangle>ACD_PI*85.0/180.0))||((Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*165.0/180.0))||((Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*75.0/180.0)))
			{
				Sun01out[2]=1;
			}
			if (((Suns[1]>0)&&(Xangle<ACD_PI*5.0/180.0))||((Suns[1]<0)&&(Xangle<ACD_PI*95.0/180.0))||((Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*15.0/180.0))||((Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*105.0/180.0)))
			{
				Sun01out[1]=1;
			}
			/*

			mtmp1 =  ((Suns[1]>0)&&(Xangle<ACD_PI*95.0/180.0));
			mtmp2 = ((Suns[1]<0)&&(Xangle<ACD_PI*5.0/180.0));
			mtmp3 = ((Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*105.0/180.0));
			mtmp4 = ((Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*15.0/180.0));
			if(mtmp1||mtmp2||mtmp3||mtmp4)
			Sun01out[0]=1;

			mtmp5 =  ((Suns[1]>0)&&(Xangle>ACD_PI*85.0/180.0));
			mtmp6 = (Suns[1]<0)&&(Xangle>ACD_PI*175.0/180.0);
			mtmp7 = (Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*75.0/180.0);
			mtmp8 = (Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*165.0/180.0);
			if(mtmp5||mtmp6||mtmp7||mtmp8)
			Sun01out[3]=1;
			mtmp9 =  (Suns[1]>0)&&(Xangle>ACD_PI*175.0/180.0);
			mtmp10 = (Suns[1]<0)&&(Xangle>ACD_PI*85.0/180.0);
			mtmp11 = (Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*165.0/180.0);
			mtmp12 = (Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle>ACD_PI*75.0/180.0);
			if(mtmp9||mtmp10||mtmp11||mtmp12)
			Sun01out[2]=1;
			mtmp13=  (Suns[1]>0)&&(Xangle<ACD_PI*5.0/180.0);
			mtmp14= (Suns[1]<0)&&(Xangle<ACD_PI*95.0/180.0);
			mtmp15 = (Suns[1]>0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*15.0/180.0);
			mtmp16 = (Suns[1]<0)&&(Zangle<ACD_PI*45.0/180.0)&&(Xangle<ACD_PI*105.0/180.0);
			if(mtmp13||mtmp14||mtmp15||mtmp16)
			Sun01out[1]=1;
			*/
		}
	}
	else
	{;}
	if (Sun01out[0]==1&&Sun01out[2]==1)
	{
		printf("Sun01out A = C");
	}
	if (Sun01out[3]==1&&Sun01out[1]==1)
	{
		printf("Sun01out B = D");
	}
}

double DotProduct(double R1[],double R2[])
/*ʸ�����*/
{
	return(R1[0]*R2[0]+R1[1]*R2[1]+R1[2]*R2[2]);
}

void GravityTorque(double Y[6],double Quatib[4],double Isat[3][3],double Td_gravity[3])
{
	/*�����ݶ�����ģ��*/
	/*��������λ���ٶ�
	�������ϵ����ϵ��̬��Ԫ��
	��������ת������
	��������ݶ�����*/
	double nine[9],Worbity,Quatio[4],Quatoi[4],Quatob[4];
	double icb1[3],icb[3],Rc,a[3],out[3];
	double ico[3]={0,0,-1};
	int i;

	PV2Six(Y,nine,&Worbity);  /*λ���ٶ�ת��������ģ�ͣ���λ�ú��ٶȼ�������Ҫ��*/
	Iner2orbitquat(nine[3],nine[2],nine[6],Quatio);  /*�������ϵ��Թ�������ϵ����̬��Ԫ��ģ�ͣ�����DW,SI��SU,����Qio*/
	QuatConjuction(Quatio,Quatoi);
	QuatProduct(Quatoi,Quatib,Quatob);  /*������ϵ����ϵ��̬��Ԫ��*/
	snormal(Quatob);
	VectorQuatFrameTrans(ico,Quatob,icb1);
	for (i=0;i<3;i++)
	{
		icb[i]=icb1[i]/LengthVector(icb1);     /*������ϵ�µ����������������߷���ĵ�λʸ��*/
	}
	Rc=sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]);  /*������ľ�*/
	MatrixProductVector(Isat,icb,a);
	CrossProduct(icb,a,out);
	for (i=0;i<3;i++)
	{
		Td_gravity[i]=3*ACD_MIUE*out[i]/Rc/Rc/Rc;  /*���������ݶ�����*/
	}
}

void AeroTorque(double Y[6],double Quatib[4],double SaeroC,double SaeroS,double Cd,double Cp[3],double F_aero[3],double Td_aero[3])
{
	/*��������ģ��*/
	/*��������λ���ٶ�ʸ��
	������̬��Ԫ��
	��������ڹ���ϵ�µ��ٶ�ʸ��
	������Чӭ�����
	��������ϵ��
	�����������ĵ�����ѹ�ĵ�ʸ��
	�������ϵ��������
	�������ϵ����������*/
	double Vi[3],Vb[3],Vaeroi[3],nine[9],Worbity,Quatio[4];
	double density;   /*�����ܶ�*/
	double Saero;/*ӭ�����*/
	int i;
	Vaero(Y,Vaeroi);   /*�����ٶ�ʸ��ģ��*/
	for(i=0;i<3;i++)
	{
		Vi[i]=Vaeroi[i]-Y[i+3];
	}
	VectorQuatFrameTrans(Vi,Quatib,Vb);  /*������ϵ�µ������ٶ�ʸ��*/
	Density(Y,ACS_density0,ACS_r0,ACS_H0,ACS_miu0,&density);   /*�����ܶ�ģ��*/
	PV2Six(Y,nine,&Worbity);
	Iner2orbitquat(nine[3],nine[2],nine[6],Quatio);  /*�������ϵ��Թ�������ϵ����̬��Ԫ��ģ��*/
	SAero(Quatib,Quatio,SaeroC,SaeroS,&Saero);/*ӭ���������ģ��*/
	for(i=0;i<3;i++)
	{
		F_aero[i]=density*Cd*Saero*LengthVector(Vb)*Vb[i]/2;   /*���㱾��ϵ��������*/
	}
	CrossProduct(Cp,F_aero,Td_aero);   /*���㱾��ϵ����������*/
}

void SunTorque(double SunInertia[3],double Quatib[4],double Psunb[3][3],double niu,double Ssun,double Psun,double Cps[3],double F_Sun[3],double Td_Sun[3])
{
	/*̫����ѹ����ģ��*/
	/*�������ϵ��̫��ʸ��
	������̬��Ԫ��
	����̫�����尲װ����
	������淴��ϵ��
	�����������
	����̫����ѹ
	����̫����ѹѹ�ĵ��������ĵ�ʸ��
	���̫����ѹ��
	���̫����ѹ����*/
	double Sunb[3],Suns0[3],Suns[3],F_Suns[3],Pbsun[3][3];
	double n[3]={0,0,1};
	double tao0[3],tao[3];
	double costheta,sintheta;
	double Ltao0;
	int i;
	VectorQuatFrameTrans(SunInertia,Quatib,Sunb);   /*���㱾��ϵ��̫��ʸ��*/
	VectorCosMatrixFrameTrans(Sunb,Psunb,Suns0);
	for (i=0;i<3;i++)
	{
		Suns[i]=Suns0[i]/LengthVector(Suns0);    /*����̫�����尲װϵ�µ�λ̫��ʸ��*/
	}
	tao0[0]=Suns[0];
	tao0[1]=Suns[1];
	tao0[2]=0;
	Ltao0=LengthVector(tao0);
	if (Ltao0!=0)
	{
	  for (i=0;i<3;i++)
	  {
		tao[i]=tao0[i]/LengthVector(tao0);    /*����̫�����尲װϵ��̫������������λʸ��*/
	  }
	}
	else
		printf("Suntorque : Length_tao0 = 0");
	costheta=-DotProduct(n,Suns);
	sintheta=sqrt(1-costheta*costheta);
	for (i=0;i<3;i++)
	{
		F_Suns[i]=-Psun*(1+niu)*Ssun*costheta*costheta*n[i]+Psun*(1-niu)*Ssun*sintheta*costheta*tao[i];
		/*����̫�����尲װϵ��̫����ѹ��*/
	}
	InvMatrix(Psunb,Pbsun);
	VectorCosMatrixFrameTrans(F_Suns,Pbsun,F_Sun);   /*���㱾��ϵ��̫����ѹ��*/
	CrossProduct(Cps,F_Sun,Td_Sun);   /*���㱾��ϵ��̫����ѹ����*/
}

void Acceler(double F_Th[3],double F_aero[3],double F_Sun[3],double Msat,double Psb[3][3],double AccelerR[3],double AccelerRC[3],double as0[3],double as[3])
{
	/*���ٶȼƲ���ģ��*/
	/*����������������
	������̬����������
	����������
	����̫����ѹ��
	������ٶȼư�װ����
	��������߼��ٶ�*/
	double F_Ths[3],F_aeros[3],F_Suns[3];
	static int ta=0;
	int i;
	VectorCosMatrixFrameTrans(F_Th,Psb,F_Ths);
	VectorCosMatrixFrameTrans(F_aero,Psb,F_aeros);
	VectorCosMatrixFrameTrans(F_Sun,Psb,F_Suns);
	srand((unsigned)time(0)+ta++);
	for (i=0;i<3;i++)
	{
		as0[i]=(F_Ths[i]+F_aeros[i]+F_Suns[i])/Msat;
		as[i]=(F_Ths[i]+F_aeros[i]+F_Suns[i])/Msat+AccelerRC[i]+((double)(rand()*2.0/RAND_MAX-1))*AccelerR[i];
	}
}

void Density(double Y[6],double density0,double r0,double H0,double miu0,double *density)
{
	/*�����ܶ�ģ��*/
	/*��������λ���ٶ�ʸ��
	��������ܶ�*/
	double r;
	r=sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]);
	*density=density0*exp(-(r-r0)/(H0+miu0*(r-r0)));
}

void Vaero(double Y[6],double Vaeroi[3])
{
	/*�����ٶ�ģ��*/
	/*��������λ���ٶ�ʸ��
	�������ϵ�´����ٶ�ʸ��*/
	double S1_Position[3];
	double Ridus;
	double costheta, sintheta, cosChijing, sinChijing;
	double VaeroBDD[3];
	int i;
	for(i=0;i<3;i++)
	{
		S1_Position[i]=Y[i];
	}
	/* �������ǵ��ľ� */
	Ridus=LengthVector(S1_Position);
	if (Ridus<ACD_Rearth)
	{
		Ridus=ACD_Rearth;
	}
	else{;}
	/* �������ǳྭ����γ */
	costheta=S1_Position[2]/Ridus;
	sintheta=sqrt(1.0-costheta*costheta);

	if (sintheta<0.001)
	{
		sintheta=0.001;
	}
	else{;}
	cosChijing=S1_Position[0]/Ridus/sintheta;
	sinChijing=S1_Position[1]/Ridus/sintheta;
	/* ����ϵ�任 Bci=Az(-Chijing)Ay(180-costheta)Bneg */
	VaeroBDD[0]=VaeroBDD[2]=0.0;
	VaeroBDD[1]=1.5*ACD_Wearth*Ridus;
	Vaeroi[0]=-cosChijing*costheta*VaeroBDD[0]-sinChijing*VaeroBDD[1]-cosChijing*sintheta*VaeroBDD[2];
	Vaeroi[1]=-sinChijing*costheta*VaeroBDD[0]+cosChijing*VaeroBDD[1]-sinChijing*sintheta*VaeroBDD[2];
	Vaeroi[2]=sintheta*VaeroBDD[0]-costheta*VaeroBDD[2];
}

void Magnet(double Bci[3],double Quatib[4],double PbBB[3][3],double MagoutR[3],double Magout0[3],double Magout[3])
{
	/*��ǿ��ģ��*/
	/*������ϵ�´ų�ǿ��ʸ��
	�����ǿ�ư�װ������
	�����ǿ�Ʋ������*/
	double Bcb[3],Bcs[3];
	double PBBb[3][3];
	static double tM=0;
	int i;
	VectorQuatFrameTrans(Bci,Quatib,Bcb);
	InvMatrix(PbBB,PBBb);
	VectorCosMatrixFrameTrans(Bcb,PBBb,Bcs);
	srand((unsigned int)time(0)+tM++);
	for(i=0;i<3;i++)
	{
		Magout0[i]=Bcs[i];
		Magout[i]=Magout0[i]+((double)(rand()*2.0/RAND_MAX-1))*MagoutR[i];
	}
}

void EarthVectori(double Y[6],double Earthi[3])
{
	/*����ϵ�µ���ʸ��ģ��*/
	/*��������λ���ٶ�ʸ��
	�������ʸ��*/
	double YY[3];
	int i;

	for(i=0;i<3;i++)
	{
		YY[i]=Y[i];
	}
	for(i=0;i<3;i++)
	{
		Earthi[i]=-YY[i]/LengthVector(YY);  /*��һ��*/
	}
}

double Angle2Double(int degrees,int minutes,double seconds)
{
	/* ���ȷ����ʽ�ĽǶ�ת���ɸ����ʽ�ĽǶ� */
	double deg;
	double min;
	double sec;
	double decimal;

	deg = degrees;
	min = minutes/60.0;
	sec = seconds/3600.0;
	decimal=sec+min+deg;
	return(decimal);
}

double Deg2Rad(double Deg)
{
	/*��ת������ģ��*/
	double Rad;
	Rad=Deg*ACD_PI/180;
	return(Rad);
}

//��������ʱ�䣬���������Ӻ���ʱ�ǣ�����λ���ٶ�[6]���������[9]��������ٶ�y���������GPS���ջ�����ϵ��[14]��GPS����ѡ���־
void GPS(double Time,double Jingdu_GreenWitch_inertia,double Y[6],double nine[9],double Worbity,double GPSR[14],unsigned char GPSSixFlag, double GPSout0[14], double GPSout[14])
{
	double EverSixE[6]; // ƽ�����������
	double Cwi[3][3];
	double Position[3],DPosition[3],PWGS[3],DPWGS[3],PMDp[3];
	double YR[6],nineR[9],WorbityR;
	static int tGPS=0;
	double Wearth[3]={0,0,7.292115855e-5};/*˲ʱ���ٶ� 7.272205217e-5,������ת���ٶȣ�ƽ̫���գ�ƽ������ת�� 7.292115855e-5*/
	int i;

	srand((unsigned)time(0)+tGPS++);
	for(i=0;i<6;i++)
	{
		YR[i]=Y[i]+((double)(rand()*2.0/RAND_MAX-1))*GPSR[i];
	}
	PV2Six(YR,nineR,&WorbityR);
	if(GPSSixFlag==1)
	{
		/*���ƽ��������*/
		SimSix2EverSix(nine,EverSixE);/*����ƽ�����������*/
		GPSout0[0]=EverSixE[0];
		GPSout0[1]=EverSixE[1];
		GPSout0[2]=EverSixE[2];
		GPSout0[3]=EverSixE[3];
		GPSout0[4]=EverSixE[4];
		GPSout0[5]=EverSixE[5];
		SimSix2EverSix(nineR,EverSixE);/*����ƽ�����������*/
		GPSout[0]=EverSixE[0];
		GPSout[1]=EverSixE[1];
		GPSout[2]=EverSixE[2];
		GPSout[3]=EverSixE[3];
		GPSout[4]=EverSixE[4];
		GPSout[5]=EverSixE[5];
	}
	else
	{
		/*���˲ʱ������*/
		GPSout0[0]=nine[0];
		GPSout0[1]=nine[1];
		GPSout0[2]=nine[2];
		GPSout0[3]=nine[3];
		GPSout0[4]=nine[4];
		GPSout0[5]=nine[5];
		/*���˲ʱ������*/
		GPSout[0]=nineR[0];
		GPSout[1]=nineR[1];
		GPSout[2]=nineR[2];
		GPSout[3]=nineR[3];
		GPSout[4]=nineR[4];
		GPSout[5]=nineR[5];
	}

	for (i=3;i<6;i++)
	{

		if(GPSout0[i]<0)
			GPSout0[i]=GPSout0[i]+2*ACD_PI;
		if(GPSout[i]<0)
			GPSout[i]=GPSout[i]+2*ACD_PI;
	}
	GPSout0[6]=-Worbity;
	GPSout[6]=-WorbityR;
	GPSout0[7]=(sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2])-6378140.0)/1000.0;/**************20100204***********************/
	Iner2WGS(Time,Jingdu_GreenWitch_inertia,Cwi);
	for(i=0;i<3;i++)
	{
		Position[i]=Y[i];
		DPosition[i]=Y[i+3];
	}
	VectorCosMatrixFrameTrans(Position,Cwi,PWGS);/*����WGS84ϵ�µ�λ��ʸ��*/
	VectorCosMatrixFrameTrans(DPosition,Cwi,DPWGS);/*����WGS84ϵ�µ��ٶ�ʸ��*/

	CrossProduct(Wearth,PWGS,PMDp);/*����ǣ���ٶ�*/
	for (i=0;i<3;i++)
	{
	   	DPWGS[i] = DPWGS[i]-PMDp[i];

	}

	for(i=0;i<3;i++)
	{
		GPSout0[i+8]=PWGS[i];
		GPSout0[i+11]=DPWGS[i];
	}


	for(i=0;i<3;i++)
	{
		Position[i]=Y[i]+((double)(rand()*2.0/RAND_MAX-1))*GPSR[i+8];
		DPosition[i]=Y[i+3]+((double)(rand()*2.0/RAND_MAX-1))*GPSR[i+11];
	}
	VectorCosMatrixFrameTrans(Position,Cwi,PWGS);/*����WGS84ϵ�µ�λ��ʸ��*/
	VectorCosMatrixFrameTrans(DPosition,Cwi,DPWGS);/*����WGS84ϵ�µ��ٶ�ʸ��*/

	CrossProduct(Wearth,PWGS,PMDp);/*����ǣ���ٶ�*/
	for (i=0;i<3;i++)
	{
	   	DPWGS[i] = DPWGS[i]-PMDp[i];

	}

	for(i=0;i<3;i++)
	{
		GPSout[i+8]=PWGS[i];
		GPSout[i+11]=DPWGS[i];
	}





}

void all2Xupdate(double AOCS_Xupdata[19],double *Time,double Wwheel[4],double ThTime[14])
{
	int i;

	*Time=AOCS_Xupdata[0];
	for(i=0;i<4;i++)
		Wwheel[i]=AOCS_Xupdata[i+1]; /*����ת��*/
	for(i=0;i<14;i++)
		ThTime[i]=AOCS_Xupdata[i+5];  /*����������ʱ��*/
}

void Xupdate2all(double Time,double Wwheel[4],double ThTime[14],double AOCS_Xupdata[19])
{
	int i;
	AOCS_Xupdata[0]=Time;
	for(i=0;i<4;i++)
		AOCS_Xupdata[i+1]=Wwheel[i]; /*����ת��*/
	for(i=0;i<14;i++)
		AOCS_Xupdata[i+5]=ThTime[i];  /*����������ʱ��*/
}

void all2Xinitial(double AOCS_Xinitial[14],double Y[6],double Wib[3],double Quatib[4],double *Msat)
{
	int i;
	for(i=0;i<6;i++)
		Y[i]=AOCS_Xinitial[i]; /*����ϵ������λ���ٶ�ʸ��*/
	for(i=0;i<3;i++)
		Wib[i]=AOCS_Xinitial[i+6];  /*����ϵ�±���ϵ��Թ���ϵ��̬���ٶ�*/
	for(i=0;i<4;i++)
		Quatib[i]=AOCS_Xinitial[i+9]; /*����ϵ������ϵ��̬��Ԫ��*/
	*Msat=AOCS_Xinitial[13]; /*��������*/
}

void Xinitial2all(double Y[6],double Wib[3],double Quatib[4],double Msat,double AOCS_Xinitial[14])
{
	int i;
	for(i=0;i<6;i++)
		AOCS_Xinitial[i]=Y[i]; /*����ϵ������λ���ٶ�ʸ��*/
	for(i=0;i<3;i++)
		AOCS_Xinitial[i+6]=Wib[i];  /*����ϵ�±���ϵ��Թ���ϵ��̬���ٶ�*/
	for(i=0;i<4;i++)
		AOCS_Xinitial[i+9]=Quatib[i]; /*����ϵ������ϵ��̬��Ԫ��*/
	AOCS_Xinitial[13]=Msat; /*��������*/
}

void MsatUPDATE(double *Msat,double DMsat[14],double ThTime[14])
{
	double M;
	int i;

	M=*Msat;
	for(i=0;i<14;i++)
		if(ThTime[i]>ACD_EPS)
			M=M-ThTime[i]*DMsat[i]/1000.0;

	if (M<ACS_MsatMin)
		M = ACS_MsatMin;
	*Msat=M;
}

void SAero(double Quatib[4],double Quatio[4],double SaeroC,double SaeroS,double *Saero)
{
	/*����ӭ�����*/
	double Quatoi[4],Quatob[4];
	double Cbo[3][3];
	double theta;
	QuatConjuction(Quatio,Quatoi);
	QuatProduct(Quatoi,Quatib,Quatob);/*������ϵ������ϵ��̬��Ԫ��*/
	snormal(Quatob);
	Quat2Tatti(Quatob,Cbo);/*���㱾��ϵ��Թ��ϵ��̬����*/
	/*������ZYX˳��ת��*/
	if (My_abs(Cbo[0][2])>1)
	{
		Cbo[0][2]=My_Sign(Cbo[0][2]);
	}
	theta=asin(-Cbo[0][2]);/*���㸩����-pi/2~pi/2*/
	*Saero=SaeroC+SaeroS*cos(theta);/*ӭ�����=������+������*/
}

void EarthShadow(double Y[6],double SunInertia[3],unsigned char *EshadowFlag)
{
	double a,b,c;
	double Position[3];
	int i;
	for(i=0;i<3;i++)
	{
		Position[i]=Y[i];
	}
	b=DotProduct(Position,SunInertia)/LengthVector(Position)/LengthVector(SunInertia);
	if(My_abs(b)>1)
	{
		b=My_Sign(b);
	}
	a=acos(b);/*����λ��ʸ����̫��ʸ���н�*/
	c=LengthVector(Position)*sin(a);
	if(c<ACD_Rearth)
	{
		if(a>ACD_PI/2)
		{
			*EshadowFlag=1;/*�����Ӱ*/
		}
		else
		{
			*EshadowFlag=0;
		}
	}
	else
	{
		*EshadowFlag=0;
	}
}

void Fthrust(double ThTime[14],double KFthurst[14],double cosattis[14][3],double Cattis[14][3],double F_Th[3],double T_Th[3])
{
	double Ft[14],F_Thi[3],T_Thi[3],Catti[3];
	int i,j;
	for(i=0;i<3;i++)
	{
		F_Th[i]=0;
		T_Th[i]=0;
	}
	for(i=0;i<14;i++)
	{
		if(ThTime[i]>ACD_EPS)
		{
			Ft[i]=KFthurst[i]*ThTime[i]; /*��������������*/
			for(j=0;j<3;j++)
			{
				F_Thi[j]=-Ft[i]*cosattis[i][j];/*cos(cosattis[i][j]*ACD_PI/180);*//*���㱾��ϵ�µ�������������*/
				Catti[j]=Cattis[i][j];
			}
			CrossProduct(Catti,F_Thi,T_Thi);    /*���㱾��ϵ�µ�������������������*/
		}
		else
		{
			for(j=0;j<3;j++)
			{
				F_Thi[j]=0;
				T_Thi[j]=0;/*������δ����,����������Ϊ0*/
			}
		}
		for(j=0;j<3;j++)
		{
			F_Th[j]=F_Th[j]+F_Thi[j];		/*���㱾��ϵ��������������*/
			T_Th[j]=T_Th[j]+T_Thi[j];		/*���㱾��ϵ��������������*/
		}
	}
}

void MoonVectori(double Time,double Y[6],double Mooni[3])
{
	/*����ʸ��ģ��*/
	/*����ʱ��
	�������ϵ������λ���ٶ�ʸ��
	�������ʸ��*/
	static double Time0=0;/*126158400.0;*//* Beijing Time 2004010100.00->UTC(s): (365.25*4-0.5)*86400.0-8*3600.0 */
	double JD2000;
	double iL,LL,TL,DWL;
	double ML,Ais,u;
	double C1,C2,C3;
	double REL,iLP,LLP,DWLP;
	double K[3]={2.223,1.161,0.324};/*,0,0.103,0.1,0.093,0.08,0.072,0.061,0.053,0.027};*/
	double Q[3]={1.002,0.825,0.012};/*,0,0.009,0.042,0.090,0.056,0.034,0.029,0.028,0};*/
	double DW[3]={0.082,-0.23,0.283};/*,-0.213,0.045,0,-0.008,-0.012,0,0,0,-2.724};*/
	double J[3]={0,-0.019,0};/*,-0.019,0,0,0,0,0,0,0,-0.241};*/
	double afa[3];
	double D,Ms,F;
	double a1=0,a2=0,a3=0,a4=0;
	double PL[3],QL[3],RL[3],Moontmp[3];
	double Cx[3][3],Z[3];
	int j;

	JD2000=(Time+Time0)/86400.0/36525.0;	/* ��ʼʱ��s need update !*/
	iL=2.0*asin(0.044751305);  /*����ư׽���*/
	LL=218.0+18.0/60.0+59.96/3600.0+(481267.0+52.0/60.0+52.833/3600.0)*JD2000-(4.787/3600.0)*JD2000*JD2000;
	LL=Deg2Rad(LL);          /*���㼸��ƽ�ƾ�*/
	TL=83.0+21.0/60.0+11.67/3600.0+(4069.0+49.36/3600.0)*JD2000-(37.165/3600.0)*JD2000*JD2000;
	TL=Deg2Rad(TL);          /*������ص�ƽ�ƾ�*/
	DWL=125.0+2.0/60.0+40.40/3600.0-(1934.0+8.0/60.0+10.266/3600.0)*JD2000+(7.476/3600.0)*JD2000*JD2000;
	DWL=Deg2Rad(DWL);        /*����������ƽ�ƾ�*/
	ML=LL-TL;                /*��������ƽ�����*/
	Ais=23.0+26.0/60.0+21.448/3600.0-(46.815/3600.0)*JD2000-(0.00059/3600.0)*JD2000*JD2000;
	Ais=Deg2Rad(Ais);        /*����Ƴཻ��*/
	D=297.0+51.0/60+0.74/3600.0+(445267.0+6.0/60.0+41.469/3600.0)*JD2000-(5.882/3600.0)*JD2000*JD2000;
	D=Deg2Rad(D);            /*����������̫����ƽ��ǣ��������*/
	Ms=357.0+31.0/60.0+44.76/3600.0+(129596581.04/3600.0)*JD2000-(0.562/3600.0)*JD2000*JD2000;
	Ms=Deg2Rad(Ms);          /*����̫����ƽ�����*/
	F=93.0+16.0/60.0+19.56/3600.0+(483202.0+1.0/60.0+3.099/3600.0)*JD2000-(12.254/3600.0)*JD2000*JD2000;
	F=Deg2Rad(F);
	/*����ϵ��*/
	C1=0.10976*sin(ML)+0.00373*sin(2.0*ML)+0.00017*sin(3.0*ML);
	C2=0.0545*cos(ML)+0.00297*cos(2.0*ML)+0.00018*cos(3.0*ML);
	C3=0.00047*sin(ML);
	afa[0]=2.0*D-ML;
	afa[1]=2.0*D;
	afa[2]=Ms+ACD_PI;
	/*afa[3]=2.0*F+ACD_PI;
	afa[4]=2*ML-2.0*D+ACD_PI;
	afa[5]=2.0*D-ML-Ms;
	afa[6]=2.0*D+ML;
	afa[7]=2.0*D-Ms;
	afa[8]=ML-Ms;
	afa[9]=D+ACD_PI;
	afa[10]=ML+Ms+ACD_PI;
	afa[11]=2.0*F-2.0*D+ACD_PI;*/

	for(j=0;j<3;j++)
	{
		a1=a1+Q[j]*cos(afa[j])*0.01;
		a2=a2+K[j]*sin(afa[j])*0.01;
		a3=a3+DW[j]*sin(afa[j])*0.01;
		a4=a4+J[j]*cos(afa[j])*0.01;
	}
	REL=60.26820751*ACD_Rearth/(1.0+C2+a1);
	LLP=LL+C1+a2;
	DWLP=DWL+C3+a3;
	iLP=iL+a4;
	u=LLP-DWLP;

	PL[0]=cos(DWLP);
	PL[1]=sin(DWLP);
	PL[2]=0;
	QL[0]=-sin(DWLP)*cos(iLP);
	QL[1]=cos(DWLP)*cos(iLP);
	QL[2]=sin(iLP);
	XCosMatrix(-Ais,Cx);
	for(j=0;j<3;j++)
	{
		Z[j]=REL*cos(u)*PL[j]+REL*sin(u)*QL[j];
	}
	MatrixProductVector(Cx,Z,RL);
	for(j=0;j<3;j++)
	{
		Moontmp[j]=RL[j]-Y[j];  /*��������ʸ��*/
	}
	for(j=0;j<3;j++)
	{
		Mooni[j]=Moontmp[j]/LengthVector(Moontmp);  /*��һ������ʸ��*/
	}
}

void XCosMatrix(double a, double q[3][3])
{
	/*x�᷽�����̬���Ҿ���ģ��*/
	q[1][1]=cos(a);
	q[2][2]=cos(a);
	q[1][2]=sin(a);
	q[2][1]=-sin(a);
	q[0][0]=1;
	q[0][1]=q[0][2]=q[1][0]=q[2][0]=0.0;
}

void SimSix2EverSix(double nine[9],double EverSixE[6])
{
	/*
	˲ʱ���������ƽ�����������ת��
	����˲ʱ������ double SimSixE[]
	���ƽ�������� double EverSixE[]
	*/
	double a,e,i,Omiga,omiga,M,f,u,p,ecosomiga,esinomiga,lamda;
	double am,em,im,Omigam,omigam,lamdam,Mm,pm;
	double tmp,tmp2,tmp3,tmp4,tmp5,esinomigam,ecosomigam,sini,cosi,sinu,cosu;

	/* ������ȡ */
	a=nine[0];
	e=nine[1];
	i=nine[2];
	Omiga=nine[3];
	omiga=nine[4];
	M=nine[5];

	if(e>(1.0-ACD_EPS))
	{
		e=1.0-ACD_EPS;
	}
	else{;}
	if(e<0.0)
	{
		e=0.0;
	}
	else{;}
	if(a<ACD_Rearth)
	{
		a=ACD_Rearth;
	}
	else{;}
	f=M2f(M,e);
	/* ���������Ǻ�γ�ȷ��� */
	u=omiga+f;
	ecosomiga=e*cos(omiga);
	esinomiga=e*sin(omiga);

	p=a*(1.0-e*e);

	tmp=ACD_Rearth/p*ACD_Rearth/p*ACD_J2;
	sini=sin(i);
	cosi=cos(i);
	sinu=sin(u);
	cosu=cos(u);

	lamda=M+omiga;

	/* see page 195��6-249�� Book ���������ȷ�� */
	am=a-1.5*tmp*p/a*p*sini*sini*(1.0-2.0*sinu*sinu);
	if(am<ACD_Rearth)
	{
		am=ACD_Rearth;
	}
	else{;}

	Omigam=Omiga-0.75*cosi*tmp*2.0*sinu*cosu;
	im=i-0.375*2.0*sini*cosi*tmp*(1.0-2.0*sinu*sinu);

	tmp2=tmp*0.25*(6.0-10.5*sini*sini);
	tmp3=tmp*0.25*(6.0-7.5*sini*sini);
	tmp4=tmp*0.875*sini*sini;
	tmp5=tmp*0.375*(3.0-5.0*cosi*cosi);

	esinomigam=esinomiga-tmp2*sinu-tmp4*sinu*(3.0-4.0*sinu*sinu);
	ecosomigam=ecosomiga-tmp3*cosu-tmp4*cosu*(4.0*cosu*cosu-3.0);
	em=sqrt(esinomigam*esinomigam+ecosomigam*ecosomigam);

	if(em>(1.0-ACD_EPS))
	{
		em=1.0-ACD_EPS;
	}
	else{;}
	if(em<0.0)
	{
		em=0.0;
	}
	else{;}
	if((My_abs(esinomigam)<ACD_EPS)&&(My_abs(ecosomigam)<ACD_EPS))
	{
		ecosomigam=1.0;
	}
	else{;}
	omigam=atan2(esinomigam,ecosomigam);
	lamdam=lamda-tmp5*2.0*sinu*cosu;
	Mm=lamdam-omigam;

	pm=am*(1.0-em*em);
	if(am<ACD_Rearth)
	{
		am=ACD_Rearth;
	}
	else{;}
	if(pm<ACD_Rearth)
	{
		pm=ACD_Rearth;
	}
	else{;}

	/* see page 195��6-248�� Book ���������ȷ�� */
	/*	dOmigam=-tmpm*nm*cosim;
	domigam=tmpm*nm*(2.0-2.5*sinim*sinim);
	dlamdam=tmpm*nm*(2.0-2.5*sinim*sinim+(1.0-1.5*sinim*sinim)*sqrt(1.0-em*em));

	dam=0.0;
	dem=0.0;*/

	EverSixE[0]=am;
	EverSixE[1]=em;
	EverSixE[2]=AngleWrap(im);
	EverSixE[3]=AngleWrap(Omigam);
	EverSixE[4]=omigam;
	EverSixE[5]=AngleWrap(Mm);

	/*EverSixE[6]=dam;
	EverSixE[7]=dem;
	EverSixE[8]=dOmigam;
	EverSixE[9]=domigam;
	EverSixE[10]=dlamdam;
	EverSixE[11]=nm;*/
}

void Iner2WGS(double Time,double Jingdu_GreenWitch_inertia,double Cwi[3][3])
{
	double JD2000;
	double P[3][3],N[3][3],R[3][3],Z[3][3];
	double PR1[3][3],PR2[3][3],PR3[3][3],PRtmp[3][3];
	double kexi,theta,z;
	double NR1[3][3],NR2[3][3],NR3[3][3],NRtmp[3][3];
	double Ais,dAis,dFai,DWL,F,D,Msun,Lsun,Lmoon,SG;
	double Ctmp1[3][3],Ctmp2[3][3];
	static double Time0=0.0;/* Beijing Time 2010010100.00->UTC(s): (365.25*4-0.5)*86400.0-8*3600.0 */
	JD2000=(Time+Time0)/86400.0/36525.0;	/* ��ʼʱ��s need update !*/

	/*����������P*/
	kexi=0.64061613889*JD2000+0.0000838556*JD2000*JD2000+(4.9994e-6)*JD2000*JD2000*JD2000;
	/*kexi=0.6406161*JD2000+0.0000839*JD2000*JD2000+0.000015*JD2000*JD2000*JD2000;*/
	kexi=kexi*ACD_d2r;
	/*theta=0.556753*JD2000-0.0001185*JD2000*JD2000-0.0000116*JD2000*JD2000*JD2000;*/
	theta=0.556753023*JD2000-0.00011851389*JD2000*JD2000-0.000011620278*JD2000*JD2000*JD2000;
	theta=theta*ACD_d2r;

	/*z=0.6406161*JD2000+0.0003040778*JD2000*JD2000+0.0000051*JD2000*JD2000*JD2000;*/
	z=0.64061613889*JD2000+0.0003040778*JD2000*JD2000+0.0000050564*JD2000*JD2000*JD2000;
	z=z*ACD_d2r;
	ZCosMatrix(-z,PR1);
	YCosMatrix(theta,PR2);
	ZCosMatrix(-kexi,PR3);
	MatrixProduct(PR1,PR2,PRtmp);
	MatrixProduct(PRtmp,PR3,P);

	/*�����¶�����N*/
	Ais=Angle2Double(23,26,21.448)-Angle2Double(0,0,46.815)*JD2000-Angle2Double(0,0,0.00059)*JD2000*JD2000;
	/*Ais=Deg2Rad(Ais);        ����Ƴཻ��*/
	Ais=Ais*ACD_d2r;

	DWL=Angle2Double(125,2,40.280)-Angle2Double(1934,8,10.539)*JD2000+Angle2Double(0,0,7.455)*JD2000*JD2000;
  /*	DWL=Deg2Rad(DWL);        ��������������ƽ�ƾ�*/
	DWL=DWL*ACD_d2r;

	F=Angle2Double(93,16,18.877)+Angle2Double(483202,1,3.137)*JD2000-Angle2Double(0,0,13.257)*JD2000*JD2000;
	/*F=Deg2Rad(F);            ��������ƽ���Ǿ�*/
	F=F*ACD_d2r;

	D=Angle2Double(297,51,1.307)+Angle2Double(445267,6,41.328)*JD2000-Angle2Double(0,0,6.891)*JD2000*JD2000;
	/*D=Deg2Rad(D);            ����������̫����ƽ��ǣ��������*/
	D=D*ACD_d2r;

	Msun=Angle2Double(357,31,39.804)+Angle2Double(35999,3,1.224)*JD2000-Angle2Double(0,0,0.577)*JD2000*JD2000;
	/*Msun=357+31.0/60+44.76/3600+129596581.04/3600*JD2000-0.562/3600*JD2000*JD2000;*/
	Msun=Msun*ACD_d2r;/* ����̫��ƽ����� */

	/*
	Mmoon =Angle2Double(134,57,46.733)+Angle2Double(477198,52,2.633)*JD2000+Angle2Double(0,0,31.310)*JD2000*JD2000;
	Mmoon = Mmoon*ACD_d2r;��������ƽ����� */

	/*IAU1980*/
	Lmoon = F+DWL;/*����ƽ�ƾ�*/
	Lsun = Lmoon-D;/*̫��ƽ�ƾ�*/
	dAis = (9.2025+0.00089*JD2000)*cos(DWL)+(0.5736-0.00031*JD2000)*cos(2.0*Lsun)+(0.0977-0.00005*JD2000)*cos(2.0*Lmoon)+(-0.0895+0.00005*JD2000)*cos(2*DWL)+(0.0054-0.00001*JD2000)*cos(Msun);
	dAis = dAis/3600*ACD_d2r;

	dFai = (-17.1996-0.01742*JD2000)*sin(DWL)+(-1.3187-0.00016*JD2000)*sin(2.0*Lsun)+(-0.2274-0.00002*JD2000)*sin(2.0*Lmoon)+(0.2026-0.00002*JD2000)*sin(2*DWL)+(0.1426-0.00034*JD2000)*sin(Msun);
	dFai = dFai/3600*ACD_d2r;
	XCosMatrix(-dAis,NR1);
	YCosMatrix(dFai*sin(Ais),NR2);
	ZCosMatrix(-dFai*cos(Ais),NR3);
	MatrixProduct(NR1,NR2,NRtmp);
	MatrixProduct(NRtmp,NR3,N);
	/*
	dAis=9.2025/3600*cos(DWL)+0.5736/3600*cos(2*F-2*D+2*DWL)+0.0977/3600*cos(2*F+2*DWL)-0.0895/3600*cos(2*DWL);
	dAis=Deg2Rad(dAis);
	dFai=-17.1996/3600*sin(DWL)-1.3187/3600*sin(2*F-2*D+2*DWL)-0.2274/3600*sin(2*F+2*DWL)+0.2062/3600*sin(2*DWL)+0.1426/3600*sin(Msun);
	dFai=Deg2Rad(dFai);

	XCosMatrix(-Ais-dAis,NR1);
	ZCosMatrix(-dFai,NR2);
	XCosMatrix(Ais,NR3);
	MatrixProduct(NR1,NR2,NRtmp);
	MatrixProduct(NRtmp,NR3,N);
	*/
	/*���������ת����R*/
	SG = Jingdu_GreenWitch_inertia+dFai*cos(Ais);/*���������ʱ*/
	ZCosMatrix(SG,R);
	/*���㼫�ƾ���Z*/
	Z[0][0]=1;
	Z[0][1]=0;
	Z[0][2]=0;/*Xp*/
	Z[1][0]=0;
	Z[1][1]=1;
	Z[1][2]=0;/*-Yp*/
	Z[2][0]=0;/*-Xp*/
	Z[2][1]=0;/*Yp*/
	Z[2][2]=1;

	/*MatrixProduct(P,N,Ctmp1);
	MatrixProduct(Ctmp1,R,Ctmp2);
	MatrixProduct(Ctmp2,Z,Cwi);*/
	/*�������ϵ��WGS84����ϵ������ת������Cwi*/
	MatrixProduct(Z,R,Ctmp1);
	MatrixProduct(Ctmp1,N,Ctmp2);
	MatrixProduct(Ctmp2,P,Cwi);
}
void YCosMatrix(double a, double q[3][3])
{
	/*y�᷽�����̬���Ҿ���ģ��*/
	q[0][0]=cos(a);
	q[2][2]=cos(a);
	q[0][2]=-sin(a);
	q[2][0]=sin(a);
	q[1][1]=1;
	q[0][1]=q[1][2]=q[1][0]=q[2][1]=0.0;
}

void ZCosMatrix(double a, double q[3][3])
{
	/*z�᷽�����̬���Ҿ���ģ��*/
	q[0][0]=cos(a);
	q[1][1]=cos(a);
	q[0][1]=sin(a);
	q[1][0]=-sin(a);
	q[2][2]=1;
	q[0][2]=q[1][2]=q[2][0]=q[2][1]=0.0;
}

void Wheel(double Wwheel[4],double Jwheel[4],double Pwheelb[4][3],double Twheelb[3])
{
	/*���㱾��ϵ�·��������������*/
	static double preWwheel[4]={0,0,0,0};
	double DWheel[4],Twheel[4];
	int i,j;
	for(i=0;i<4;i++)
	{
		DWheel[i]=(Wwheel[i]-preWwheel[i])/0.05;/*����ת������*/
		Twheel[i]=-Jwheel[i]*DWheel[i];/*���㵥�����ֲ�������*/
		if(My_abs(Twheel[i])>ACS_Twheel_MAX)
			Twheel[i]=ACS_Twheel_MAX*My_Sign(Twheel[i]);
	}
	for(i=0;i<3;i++)
	{
		Twheelb[i]=0;
		for(j=0;j<4;j++)
		{
			Twheelb[i]=Twheelb[i]+Twheel[j]*Pwheelb[j][i];
		}
	}
	for(i=0;i<4;i++)
	{
		preWwheel[i]=Wwheel[i];
	}

      printf(" Tx :%3f , Ty: %3f , Tz : %3f, Ts :%3f \n",My_abs(Twheel[0]),Twheel[1],Twheel[2],Twheel[3] );
}

void Quat2Euler(double Quatib[4],double Quatio[4],double Euler[3])
{
	double Quatoi[4],Quatob[4];
	double Cbo[3][3],cosfai,sintheta,costheta,sinpesi,cospesi;


	QuatConjuction(Quatio,Quatoi); //ȡ����
	QuatProduct(Quatoi,Quatib,Quatob); //��Ԫ���˷�����õ�Quatob
	snormal(Quatob);
	Quat2Tatti(Quatob,Cbo); //��Ԫ�����������Ҿ���ı任
	if (My_abs(Cbo[1][2])>1)
	{
		Cbo[1][2]=My_Sign(Cbo[1][2]);
	}
	/*��z,x,y˳��ת��*/
	Euler[0]=asin(Cbo[1][2]); /*-pi/2~pi/2*/
	cosfai=sqrt(1-Cbo[1][2]*Cbo[1][2]);
	if (My_abs(cosfai)>ACD_EPS)
	{
		sintheta=-Cbo[0][2]/cosfai;
		costheta=Cbo[2][2]/cosfai;
		sinpesi=-Cbo[1][0]/cosfai;
		cospesi=Cbo[1][1]/cosfai;
		Euler[1]=dftan(costheta,sintheta); /*0~2pi*/
		Euler[2]=dftan(cospesi,sinpesi); /*0~2pi*/
	}
	else
	{
		cosfai=ACD_EPS*(My_Sign(cosfai));
		sintheta=-Cbo[0][2]/cosfai;
		costheta=Cbo[2][2]/cosfai;
		sinpesi=-Cbo[1][0]/cosfai;
		cospesi=Cbo[1][1]/cosfai;
		Euler[1]=dftan(costheta,sintheta); /*0~2pi*/
		Euler[2]=dftan(cospesi,sinpesi); /*0~2pi*/
	}
	Euler[0] = AngleWrap(Euler[0]);
	Euler[1] = AngleWrap(Euler[1]);
	Euler[2] = AngleWrap(Euler[2]); /*******20100204********/
	/*
	Euler[0]=asin(Cbo[1][2]);
	Euler[1]=atan(-Cbo[0][2]/Cbo[2][2]);
	Euler[2]=atan(-Cbo[1][0]/Cbo[1][1]);
	*/
	/*	Euler[0]=Euler[0]/10.0;
	Euler[1]=Euler[1]/10.0;
	Euler[2]=Euler[2]/10.0;*/
}


void Magnetic(double Bci[3],double Quatib[4],double Mags[3],double PMagb[3][3],double Tmag[3])
{

	double Magb[3],Bcb[3],PbMag[3][3];
	InvMatrix(PMagb,PbMag);
	VectorCosMatrixFrameTrans(Mags,PbMag,Magb);
	VectorQuatFrameTrans(Bci,Quatib,Bcb);
	CrossProduct(Magb,Bcb,Tmag);
}


//����ϵͳ�����������[100]��ϵͳ���������־[6]
void Output(double AOCS_Y0[100],unsigned char AOCS_Flag0[6],double AOCS_Y[32],unsigned char AOCS_Flag[9])
{
	/*ϵͳ���*/
	int i,Nput,Mput;

	/*�������*/

	for(i=0;i<6;i++)
	{
		outputLAN.strY[i]=(float)Y[i];  /*�������ϵ������λ���ٶ�ʸ��*/
	}
	for(i=0;i<9;i++)
	{
		outputLAN.strNine[i]=(float)nine[i];  /*����������*/
	}
	for(i=0;i<4;i++)
	{
		outputLAN.strQuatib[i]=(float)Quatib[i];  /*�������ϵ��Թ���ϵ��̬��Ԫ��*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strWib[i]=(float)Wib[i];  /*�������ϵ�±���ϵ��Թ���ϵ��̬���ٶ�*/
	}
	for(i=0;i<1;i++)
	{
		outputLAN.strMsat=(float)Msat;  /*�����������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strF_Th[i]=(float)F_Th[i];  /*�������ϵ�������������������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strT_Th[i]=(float)T_Th[i];  /*�������ϵ���������������������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strTwheelb[i]=(float)Twheelb[i];  /*�������ϵ�·����������������*/
	}
	for(i=0;i<1;i++)
	{
		outputLAN.strFsun1out0=(float)Fsun1out0;  /*�������̫��������1�������*/
	}
	for(i=0;i<1;i++)
	{
		outputLAN.strFsun2out0=(float)Fsun2out0;  /*�������̫��������2�������*/
	}
	for(i=0;i<4;i++)
	{
		outputLAN.strQuatsiA0[i]=(float)QuatsiA0[i];  /*�����������A�������*/
	}
	for(i=0;i<4;i++)
	{
		outputLAN.strQuatsiB0[i]=(float)QuatsiB0[i];  /*�����������B�������*/
	}
	for(i=0;i<4;i++)
	{
		outputLAN.strWgyro10[i]=(float)Wgyro10[i];  /*���ƽ̨���ݲ������*/
	}
	for(i=0;i<4;i++)
	{
		outputLAN.strWgyro20[i]=(float)Wgyro20[i];  /*����������ݲ������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.stras0[i]=(float)MagU[i];  /*��������������Ƶ�ѹ*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strMagout0[i]=(float)Tmag[i];  /*�������������������*/
	}
	for(i=0;i<14;i++)
	{
		outputLAN.strGPSout0[i]=(float)GPSout0[i];  /*���GPS�������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strSunInertia[i]=(float)SunInertia[i];  /*�������ϵ��̫��ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strMooni[i]=(float)Mooni[i];  /*�������ϵ������ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strEarthi[i]=(float)Earthi[i];  /*�������ϵ�µ���ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strBci[i]=(float)Bci[i];  /*�������ϵ����������λ�õĴų�ǿ��*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strTd_geomagnetic[i]=(float)Td_geomagnetic[i];  /*�������ϵ�µشŸ�������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strTd_gravity[i]=(float)Td_gravity[i];  /*�������ϵ�������ݶȸ�������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strTd_aero[i]=(float)Td_aero[i];  /*�������ϵ��������������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strTd_Sun[i]=(float)Td_Sun[i];  /*�������ϵ��̫����ѹ��������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strF_aero[i]=(float)F_aero[i];  /*�������ϵ������������*/
	}
	for(i=0;i<3;i++)
	{
		outputLAN.strF_Sun[i]=(float)F_Sun[i];  /*�������ϵ��̫����ѹ������*/
	}
	outputLAN.fEuler[0] =(float)(Euler[0]*57.3);
	outputLAN.fEuler[1] =(float)(Euler[1]*57.3);
	outputLAN.fEuler[2] =(float)(Euler[2]*57.3);
	outputLAN.fHbody[0] =(float)(Hbody[0]);
	outputLAN.fHbody[1] =(float)(Hbody[1]);
	outputLAN.fHbody[2] =(float)(Hbody[2]);

	outputLAN.strEshadowFlag=EshadowFlag;/*�����Ӱ��־*/
	outputLAN.strFsun1Flag=Fsun1Flag;/*����̫��������1��־λ*/
	outputLAN.strFsun2Flag=Fsun2Flag;/*����̫��������2��־λ*/

	outputLAN.strStarAFlag=StarAFlag;/*��������A��־λ*/
	outputLAN.strStarBFlag=StarBFlag;/*��������B��־λ*/



	/*�������*****************************************************/
	Nput=0;
	for(i=0;i<6;i++)
	{
		AOCS_Y0[Nput++]=Y[i];  /*�������ϵ������λ���ٶ�ʸ��*/
	}
	for(i=0;i<9;i++)
	{
		AOCS_Y0[Nput++]=nine[i];  /*����������*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y0[Nput++]=Quatib[i];  /*�������ϵ��Թ���ϵ��̬��Ԫ��*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Wib[i];  /*�������ϵ�±���ϵ��Թ���ϵ��̬���ٶ�*/
	}
	for(i=0;i<1;i++)
	{
		AOCS_Y0[Nput++]=Msat;  /*�����������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=F_Th[i];  /*�������ϵ�������������������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=T_Th[i];  /*�������ϵ���������������������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Twheelb[i];  /*�������ϵ�·����������������*/
	}
	for(i=0;i<1;i++)
	{
		AOCS_Y0[Nput++]=Fsun1out0;  /*�������̫��������1�������*/
	}
	for(i=0;i<1;i++)
	{
		AOCS_Y0[Nput++]=Fsun2out0;  /*�������̫��������2�������*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y0[Nput++]=QuatsiA0[i];  /*�����������A�������*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y0[Nput++]=QuatsiB0[i];  /*�����������B�������*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y0[Nput++]=Wgyro10[i];  /*���ƽ̨���ݲ������*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y0[Nput++]=Wgyro20[i];  /*����������ݲ������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=MagU[i];  /*��������������Ƶ�ѹ*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Tmag[i];  /*���������������*/
	}
	for(i=0;i<14;i++)
	{
		AOCS_Y0[Nput++]=GPSout0[i];  /*���GPS�������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=SunInertia[i];  /*�������ϵ��̫��ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Mooni[i];  /*�������ϵ������ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Earthi[i];  /*�������ϵ�µ���ʸ��*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Bci[i];  /*�������ϵ����������λ�õĴų�ǿ��*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Td_geomagnetic[i];  /*�������ϵ�µشŸ�������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Td_gravity[i];  /*�������ϵ�������ݶȸ�������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Td_aero[i];  /*�������ϵ��������������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=Td_Sun[i];  /*�������ϵ��̫����ѹ��������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=F_aero[i];  /*�������ϵ������������*/
	}
	for(i=0;i<3;i++)
	{
		AOCS_Y0[Nput++]=F_Sun[i];  /*�������ϵ��̫����ѹ������*/
	}
	Mput=0;
	AOCS_Flag0[Mput++]=EshadowFlag;/*�����Ӱ��־*/
	AOCS_Flag0[Mput++]=Fsun1Flag;/*����̫��������1��־λ*/
	AOCS_Flag0[Mput++]=Fsun2Flag;/*����̫��������2��־λ*/

	AOCS_Flag0[Mput++]=StarAFlag;/*��������A��־λ*/
	AOCS_Flag0[Mput++]=StarBFlag;/*��������B��־λ*/

	/*�����ģ����*/
	Nput=0;
	for(i=0;i<1;i++)
	{
		AOCS_Y[Nput++]=Fsun1out;  /*�������̫��������1�������[0]*/
	}
	for(i=0;i<1;i++)
	{
		AOCS_Y[Nput++]=Fsun2out;  /*�������̫��������2�������[1]*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y[Nput++]=Wgyro1[i];  /*���ƽ̨���ݲ������[2--5]*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y[Nput++]=Wgyro2[i];  /*����������ݲ������[6--9]*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y[Nput++]=QuatsiA[i];  /*�����������A�������[10--13]*/
	}
	for(i=0;i<4;i++)
	{
		AOCS_Y[Nput++]=QuatsiB[i];  /*�����������B�������[14--17]*/
	}

    for(i=0;i<1;i++)
	{
		AOCS_Y[Nput++] = Distance;  /*�������״�Ĳ�������*/
	}
	for(i=0;i<1;i++)
	{
		AOCS_Y[Nput++] = Theta_yz;  /*����������Ĳ����Ƕ�*/
	}
    for(i=0;i<1;i++)
	{
		AOCS_Y[Nput++] = Theta_yx;  /*����������Ĳ����Ƕ�*/
	}

	Mput=0;
	AOCS_Flag[Mput++]=Fsun1Flag;/*����̫��������1��־λ*/
	AOCS_Flag[Mput++]=Fsun2Flag;/*����̫��������2��־λ*/
	AOCS_Flag[Mput++]=StarAFlag;/*��������A��־λ*/
	AOCS_Flag[Mput++]=StarBFlag;/*��������B��־λ*/

	int k;
	for(k=0; k<21; k++)
        printf("AOCS_Y[%d] is %lf\n", k, AOCS_Y[k]);
    for(k=0; k<4; k++)
        printf("AOCS_Flag[%d] is %d\n", k, AOCS_Flag[k]);
}

void Distance_in_Radar(double Sat_r1[6], double Sat_r2[6], double *Distance)
{
	double delta_r[3];
	int i;
	double temp = 0;
	for(i=0; i<3; i++)
		delta_r[i] = Sat_r1[i] - Sat_r2[i];
	for(i=0; i<3; i++)
		temp += delta_r[i]*delta_r[i];
	*Distance = sqrt(temp);
}

void Angel_in_Camera(double Sat_r1[3], double Sat_r2[3],double Starsb[4],double Quatib[4],double *Theta_yz,double *Theta_xy)
{
	/* �������������������С */
	double Starbs[4];	//���嵽����ϵ����Ԫ����
	//double StarZ[3]={0,0,1};/*������������*/
	double Delta_ri[3], Delta_rb[3], Delta_rs[3];
	int i;
	for(i=0; i<3; i++)
        Delta_ri[i] = Sat_r2[i] - Sat_r1[i];

	/*������������Ч��*/
	QuatConjuction(Starsb,Starbs);
	VectorQuatFrameTrans(Delta_ri,Quatib,Delta_rb);
	VectorQuatFrameTrans(Delta_rb,Starbs,Delta_rs);

	*Theta_yz = ACD_PI - atan(Delta_rs[2]/Delta_rs[1]);
	*Theta_xy = ACD_PI - atan(Delta_rs[0]/Delta_rs[1]);
	printf("%lf\t, %lf", *Theta_xy, *Theta_yz);
	//*StarFlag=0;

}
