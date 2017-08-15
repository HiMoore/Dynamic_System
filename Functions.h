#ifndef __Functions_h__
#define __Functions_h__

extern void Xupdate2all(double Time,double Wwheel[4],double ThTime[14],double AOCS_Xupdata[19]);
extern void all2Xupdate(double AOCS_Xupdata[19],double *Time,double Wwheel[4],double ThTime[14]);
extern void Xinitial2all(double Y[6],double Wib[3],double Quatib[4],double Msat,double AOCS_Xinitial[14]);
extern void all2Xinitial(double AOCS_Xinitial[14],double Y[6],double Wib[3],double Quatib[4],double *Msat);

extern void Output(double AOCS_Y0[100],unsigned char AOCS_Flag0[6],double AOCS_Y[32],unsigned char AOCS_Flag[9]);

extern void GeomagneticTorque(double Bci[3],double Quatib[4],double Mdis[3],double Td_geomagnetic[3]);
extern void GravityTorque(double Y[6],double Quatib[4],double Isat[3][3],double Td_gravity[3]);
extern void AeroTorque(double Y[6],double Quatib[4],double SaeroC,double SaeroS,double Cd,double Cp[3],double F_aero[3],double Td_aero[3]);
extern void SunTorque(double SunInertia[3],double Quatib[4],double Psunb[3][3],double niu,double Ssun,double Psun,double Cps[3],double F_Sun[3],double Td_Sun[3]);
extern void Magnetic(double Bci[3],double Quatib[4],double Mags[3],double PMagb[3][3],double Tmag[3]);

extern void Geomagnetic(double Y[6],double Jingdu_GreenWitch_inertia,double gnm[7][7],double hnm[7][7],unsigned char Number_mag,double Bci[3]);
extern void SunVectori(double Time,double *Yousun,double *Aisun,double *Jingdu_GreenWitch_inertia,double SunInertia[],double Quatis[]);
extern void EarthVectori(double Y[6],double Earthi[3]);
extern void MoonVectori(double Time,double Y[6],double Mooni[3]);
extern void Density(double Y[6],double density0,double r0,double H0,double miu0,double *density);
extern void Vaero(double Y[6],double Vaeroi[3]);
extern void EarthShadow(double Y[6],double SunInertia[3],unsigned char *EshadowFlag);

extern void orbJ2_diff(double Y[6],double ai[3],double DY[6]);
extern void orbJ2_rk4(double Y[6],double H,double Fb[3],double Quatib[4],double Msat);
extern void PV2Six(double Y[6],double nine[9],double *Worbity);/*nine[9]:A,Ex,SI,DW,SW,DM,SU,theta,n*/
extern void PV2Six_OLD(double Y[6],double nine[9],double *Worbity);/*nine[9]:A,Ex,SI,DW,SW,DM,SU,theta,n*/
extern void SimSix2EverSix(double nine[9],double EverSixE[6]);
extern void Quat2Euler(double Quatib[4],double Quatio[4],double Euler[3]);

extern void BiasattiDynamic_rk4(double H,double Hbody[3],double Tb[3],double I[3][3],double Wib[3]);
extern void BiasattiDynamic_diff(double T[3],double I[3][3],double Hbody[3],double W[3],double DW[3]);
extern void QuatTatti_diff(double Q[4],double W_gyro[3],double DQ[4]);
extern void QuatTatti_rk4(double H,double W_gyro[3],double Q[4]);

extern void quatstar(double SunInertia[3],double Earthi[3],double Mooni[3],double StarLimit[3],double Starsb[4],double Quatib[4],double quatstarR[3],unsigned char *StarFlag,double Quatsi0[4],double Quatsi[4]);
extern void Gyro(double w[3],double Pgyrob[3][3],unsigned char NoGyro[3],double GyroR[3],double GyroRC[3],double Wgyro0[3],double Wgyro[3]);
extern void Fsun(double SunInertia[3],double Quatib[4],double PbFsun[3][3],double FsunLimit,double FsunR,unsigned char EshadowFlag,unsigned char *FsunFlag,double *Fsunout0,double *Fsunout);
extern void Sun01(double SunInertia[3],double Quatib[4],double Pbsun01[3][3],unsigned char EshadowFlag,unsigned char Sun01out[5]);
extern void GPS(double Time,double Jingdu_GreenWitch_inertia,double Y[6],double nine[9],double Worbity,double GPSR[14],unsigned char GPSSixFlag,double GPSout0[14],double GPSout[14]);
extern void Acceler(double F_Th[3],double F_aero[3],double F_Sun[3],double Msat,double Psb[3][3],double AccelerR[3],double AccelerRC[3],double as0[3],double as[3]);
extern void Magnet(double Bci[3],double Quatib[4],double PbBB[3][3],double MagoutR[3],double Magout0[3],double Magout[3]);

extern void MsatUPDATE(double *Msat,double DMsat[14],double ThTime[14]);

extern void Wheel(double Wwheel[4],double Jwheel[4],double Pwheelb[4][3],double Twheelb[3]);
extern void Fthrust(double ThTime[14],double KFthurst[14],double cosattis[14][3],double Cattis[14][3],double F_Th[3],double T_Th[3]);

extern void SAero(double Quatib[4],double Quatio[4],double SaeroC,double SaeroS,double *Saero);
extern double dftan(double X,double Y);
extern void Iner2orbitquat(double Omega,double i,double u,double Quatio[4]);
extern double M2f(double M, double e);
extern void Iner2WGS(double Time,double Jingdu_GreenWitch_inertia,double Cwi[3][3]);
extern void WheelMomentum(double Jwheel[4],double Wwheel[4],double MomentumWheel[4]);
extern void WheelSysToThree(double WheelSys[4],double Pwheelb[4][3],double WheelThree[3]);
extern void formHB(double Hwheel[3],double Wbody[3],double IBf[3][3],double Hbody[3]);
extern void XQuat(double a,double q[4]);
extern void YQuat(double a,double q[4]);
extern void ZQuat(double a,double q[4]);
extern void QuatProduct(double l[4],double p[4],double q[4]);
extern void snormal(double Quat[4]);
extern void MatrixProductVector(double b[3][3],double c[3],double a[3]);
extern void InvMatrix(double Mat[3][3],double InvMat[3][3]);
extern void MatrixProduct(double b[3][3],double c[3][3],double a[3][3]);
extern void antisymMatrix(double a[3],double a_mat[3][3]);
extern double LengthVector(double Vector[]);
extern void CrossProduct(double Vin1[],double Vin2[],double Vout[]);
extern double DotProduct(double R1[3],double R2[3]);
extern void VectorQuatFrameTrans(double Rin[],double QuatRin2Rout[],double Rout[]);
extern void VectorCosMatrixFrameTrans(double Rin[],double CosMatrixOutIn[3][3],double Rout[]);
extern void QuatConjuction(double Quat[4],double QuatConj[4]);
extern void Quat2Tatti(double Quat[], double Tatti[3][3]);
extern void Tatti2Quat(double Tatti[3][3], double Quat[4]);
extern double My_abs(double x);
extern int My_Sign(double x);
extern double AngleWrap(double u);
extern double Angle2Double(int degrees,int minutes,double seconds);
extern double Deg2Rad(double Deg);
extern void XCosMatrix(double a, double q[3][3]);
extern void ZCosMatrix(double a, double q[3][3]);
extern void YCosMatrix(double a, double q[3][3]);

extern void Distance_in_Radar(double Sat_r1[6], double Sat_r2[6], double *Distance);
extern void Angel_in_Camera(double Sat_r1[3], double Sat_r2[3],double Starsb[4],double Quatib[4],double *Theta_yz,double *Theta_xy);

#endif


