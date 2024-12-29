#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <omp.h>
const int E=20, UN=200, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel;

//void ecfcal(int,double[]);

int  N, trim, check, checkl, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, it, nf,  block[3], blockd, jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct, rul, sla, veei, veec, SHA, ptr, ptrr, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt, gv, bvs, bvsi, flag, varho, vart, iflag, flagx;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, cp, cph, cpsh, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9], f[3], bl1, bl2, bl3, bl[3], f1, f2, f3, ar, per; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, fe, tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON, eb, teb, tbl, pbc,  dif, difu, difus, difl, difusl, adefus, dive, divea, divei, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm, hintgrl, ene, vol, sene, nene, tnene, intgrl1, intgrl2, ene2, sene2, inv[11][2], liml, limu, inv1[10][2], inv2[10][2], inv3[10][2], inv4[11][2], inv5[10][2], fintgrl, svali[100][2], hyna1[10][2], hyna2[10][2], hyna3[10][2], viper1[10][2], spider1[10][2], spider2[10][2], spider3[10][2], spider4[10][2],  nematicd[3], rdistance, ofact, range_1, range_2;  

//double subroutine_integral(double, double, double, double, double, int);

double subroutine_integral(double liml, double limu, double inv[][2], double* fintgrl, double svali[][2], int iflag);





int main(void)
{

	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	fscanf(bu,"%lf",&bdm);
	N=(int)round(bdm);
	//if (N==0) printf("raam tere desh men \n");
	// Definition of universal constants starts!
	
	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;
	
	
	fscanf(bu,"%lf",&bdm);
	step=bdm;
	fscanf(bu,"%lf",&bdm);
	tsim=(int)round(bdm);
	
	//tsim=100000;
	ivf=0.004;
	
	
	
	
	//fvf= 0.02 ;
	
	fscanf(bu,"%lf",&bdm);
	fvf=bdm;
	
	//fscanf(bu,"%lf",&fvf);
	
	
	//fvf= 0.02 ;
	
	
	fscanf(bu,"%lf",&bdm);
	nematicd[2]=bdm;
	
	
	//fscanf(bu,"%lf",&fvf);
	
	
	
	
	
	
	int war;
	
	//fscanf(bu,"%d",&war);
	fscanf(bu,"%lf",&bdm);
	war=(int)round(bdm);
	
	//war=3;
	
	 
	double* ar=malloc(war*sizeof(double));
	
	double* sar=malloc(war*sizeof(double));
	
	double* per=malloc(war*sizeof(double));
	
	double* sp=malloc(war*sizeof(double));
	
	int* nt=malloc(war*sizeof(int));
	
	int* sco=malloc(war*sizeof(int));
	
	
	int** pop=malloc(war*sizeof(int*));
	for (i=0; i<war; i++)
	{
		pop[i]=malloc(war*sizeof(int));
	}
	
	// printf("raam tere desh men \n");
	
	
	
	
	for (i=0; i<war; i++)
	{
		sco[i]=0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			sco[i]=(int)round(bdm);
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		ar[i]=1.0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			ar[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		sar[i]=1.0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			sar[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	//per[0]=1.0;
	for (i=0; i<war-1; i++)
	{
		per[i]=0.0;
	}
	
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			per[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	red=0.0;
	for (i=0; i<war-1; i++)
	{
		red=red+per[i];
	}
	
	per[war-1]=1.0-red;
	
	
	for (i=0; i<war  ; i++)
	{
		nt[i]=(int)round((float)N*per[i]);
		//printf("%d  \n", nt[i]);
	}
	
	
	//Last fraction        !!!!!!!

	
	
	
	// Shape parameters 0 and < 1.0 ; Oblate and 1 and > 1.0 ; Prolate          !!!!!!!!!!!!!! 
	
	SHA = 1 ;
	
	
	
	
	// Definition of simulation constants ends        !!!!!!!!!!!
	
	
	
	// Definition of compression parameter starts          !!!!!!!!!!
	
	veec=200;
	vees=0.005;    // Compression steps!
	vel=(int)(100.0*(fvf/0.5));
	//vel=10;
	printf("%ld  \n", vel);
	
	// Definition of compression parameter ends            !!!!!!!!!!
	
	

	// Declaration of core bond variable starts          !!!!!!!!!!
	
	double epsi;
	epsi=0.0;
	
	// Declaration of core bond variable ends         !!!!!!!!!!
	
	
	
	// Evaluation of bond variable starts           !!!!!!!!!
	
	
	
	bn=40;
	bnll=40;
	sst=1;
	
	
	//double epsi1, epsi2, epsi3, omega2, beta1, beta2, beta3, betai;
	
	// Declaration of core bond variable ends!
	
	
	// Evaluation of bond variable starts!
	


	fclose(bu);
	remove("Base_data.txt");



	
	
	
	
	
	
	
	
	
	// Evaluation of bond variable ends!
	
	printf("Number of particles-    %d \n", N);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Step size-   %lf \n", step);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of simulation step-   %ld \n", tsim);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Starting volume fraction-   %lf\n", ivf);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Volume fraction-    %lf \n", fvf);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Type of particles-    %d \n", war);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Shape factor (aspect ratios) are given by \n");
	for (i=0; i<war; i++)
	{
	
		printf("Aspect ratio for particle type %d is given by-   %lf\n", i+1, ar[i]);
		
		printf("Size factor for particle type %d is given by-   %lf\n", i+1, sar[i]);
	
	}
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Particles fractions are given by \n");
	for (i=0; i<war; i++)
	{
		printf("Fraction for particle type %d is given by-   %lf\n", i+1, per[i]);
		
	}
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of particles for each type (aspect ratios) are given by \n");
	for (i=0; i<war; i++)
	{
		printf("Number of particles for particle type %d is given by-   %d\n", i+1, nt[i]);
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("Potential parameters are given by:- \n");
	printf("\n");
			 printf("\n");
			 printf("\n");
	
	
	printf("Selected Epsi to harnesh the particle in the range %lf \n", epsi);
	//printf("Values corresponding to the particular parameter is given by\n");
	
	
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	
	double* d=malloc(war*sizeof(double));
	
	double* r=malloc(war*sizeof(double));
	
	double* rsm=malloc(war*sizeof(double));
	
	double* rm=malloc(war*sizeof(double));
	
	double* rs=malloc(war*sizeof(double));
	
	double** e=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		e[i]=malloc(3*sizeof(double)); 
	
	}
	
	double** ds=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		ds[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** dss=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		dss[i]=malloc(war*sizeof(double)); 
	
	}
	
	int** pps=malloc(war*sizeof(int*));
	for(i=0; i<war; i++)
	{
		pps[i]=malloc(war*sizeof(int)); 
	
	}
	
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			if (ar[i]==1.0 && ar[j]==1.0) pps[i][j]=1;
			else pps[i][j]=0;
		
		}
	
	}
	
	
	
	for (i=0; i<war; i++)
	{
		if (sco[i]==0)
		{
			
			d[i]=sar[i]*pow(ar[i],(2.0/3.0));
			r[i]=d[i]/2.0;
			rs[i]=r[i]*r[i];
			rsm[i]=rs[i]/(ar[i]*ar[i]);
			rm[i]=r[i]/ar[i];
			
			e[i][0]=rm[i]+(epsi/2.0);
			e[i][1]=rm[i]+(epsi/2.0);
			e[i][2]=r[i]+(epsi/2.0);
			
			if(i==0)
			{
				d[i]=2.0*sar[i]*ar[i];
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=rm[i]+(epsi/2.0);
				e[i][1]=rm[i]+(epsi/2.0);
				e[i][2]=r[i]+(epsi/2.0);
			
			}
			
		}
		else
		{
			if (ar[i]>1.0)
			{
				d[i]=sar[i]*1.0;
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=rm[i]+(epsi/2.0);
				e[i][1]=rm[i]+(epsi/2.0);
				e[i][2]=r[i]+(epsi/2.0);
			}
			else
			{
				d[i]=sar[i]*ar[i];
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=rm[i]+(epsi/2.0);
				e[i][1]=rm[i]+(epsi/2.0);
				e[i][2]=r[i]+(epsi/2.0);
			
			}
		}
		
	}	
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			if (rm[i]<r[i]) 
			{ 
				ds[i][j]=rm[i];
				dss[i][j]=r[i];
			}
			else
			{
				ds[i][j]=r[i];
				dss[i][j]=rm[i];
			}
			
			if (rm[j]<r[j]) 
			{ 
				ds[i][j]=ds[i][j]+rm[j];
				dss[i][j]=dss[i][j]+r[j];
			}
			else
			{
				ds[i][j]=ds[i][j]+r[j];
				dss[i][j]=dss[i][j]+rm[j];
			}
			
			dss[i][j]=dss[i][j]+epsi;
			
			dss[i][j]=dss[i][j]*dss[i][j];
			ds[i][j]=ds[i][j]*ds[i][j];
			
			
		}
		
		
	}
	printf("dsdsdsdsds%lf\n", ds[0][0]);
	fe=0.0;
	for (i=0; i<war; i++)
	{
		fe=fe+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*fvf));
	}
	fe=pow(fe,(1.0/3.0));
	
	f[0]=0.0;
	
	for (i=0; i<war; i++)
	{
		f[0]=f[0]+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*ivf));
	}
	
	f[0]=pow(f[0],(1.0/3.0));
	f[1]=f[0];
	f[2]=f[0];
	
	ble=1.0;
	for (i=0; i<war; i++)
	{
		if (rm[i]<r[i])
		{
			if (d[i]>ble) ble=d[i];
		}
		else
		{
			if ((d[i]/ar[i])>ble) ble=d[i]/ar[i];
		}
	}
	
	ble=ble+epsi;
	
	
	
	
	
	
	bl[0]=ble;
	bl[1]=ble;
	bl[2]=ble;
	
	block[0]=(int)floor(f[0]/bl[0]);
	bl[0]=f[0]/(float)block[0];
	f[0]=bl[0]*(float)block[0];
	
	block[1]=(int)floor(f[1]/bl[1]);
	bl[1]=f[1]/(float)block[1];
	f[1]=bl[1]*(float)block[1];
	
	block[2]=(int)floor(f[2]/bl[2]);
	bl[2]=f[2]/(float)block[2];
	f[2]=bl[2]*(float)block[2];
	
	bl1=bl[0];
	bl2=bl[1];
	bl3=bl[2];
	
	blockd=block[0]+20;
	f1=f[0];
	f2=f[1];
	f3=f[2];
	
	printf("                                ******:::::::   Calculated parameters are given by   :::::::******\n");
	printf("\n");
	printf("Starting edge is given by-  %lf\n", f[0]);
	printf("\n");
	printf("\n");
	printf("Aimed edge is given by-  %lf\n", fe);
	printf("\n");
	printf("\n");
	printf("Starting number of shell at each edge is given by-  %d\n", block[0]);
	printf("\n");
	printf("\n");
	printf("Shell length is given by-   %lf \n\n\n", ble);
	
	pbc=ble+0.5;
	
	for (i=0; i<war; i++)
	{
		printf("For family of particle %d-\n", i);
		printf("Diameter-  %lf\n", d[i]);
		printf("Length of semi-symmetry-axis-  %lf\n", r[i]);
		printf("Length of semi-axis perpendicular to the symmetry-axis-  %lf \n \n \n \n", rm[i]);
	
	}

	// Evaluation of simulation variables ends!
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Declaration of core variables starts!
	
	int pit[UN], wpi;
	
	
	double** li=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		li[i]=malloc(9*sizeof(double));
	}
	
	
	double** el=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		el[i]=malloc(9*sizeof(double));
	}
	
	double** ami=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		ami[i]=malloc(3*sizeof(double));
	}
	int* pid=malloc(N*sizeof(int));
	
	
	double ldcxyz[N][3], ldcxy[N][2], ldcz[N], lom[N][3];
	
	
	
	int** nrdi=malloc(N*sizeof(int*));
	int* wnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		nrdi[i]=malloc(UN*sizeof(int));
	}
	
	
	
	
	
	
	// Declaration of compresser variables starts!
	
	int vee;
	
	// Declaration of compresser variables ends!
	
	// Declaration of core variables ends!
	
	
	
	
	
	
	// Declaration of bonding variable starts!
	
	
	
	
	int** bobb=malloc(N*sizeof(int*));
	int* wbobb=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		bobb[k]=malloc(E*sizeof(int));
	
	}	
	
	
	int** cbobb=malloc(N*sizeof(int*));
	int* wcbobb=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		cbobb[k]=malloc(E*sizeof(int));
	
	}	
	
	
	int wpbobb;
	int psbobb[E];

	int hold[N][E], whold[N], bumbque[E];
	

	int** bobbi=malloc(N*sizeof(int*));
	int* wbobbi=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		bobbi[k]=malloc(E*sizeof(int));
	
	}	
		
	
	
	int** cbobbi=malloc(N*sizeof(int*));
	int* wcbobbi=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		cbobbi[k]=malloc(E*sizeof(int));
	
	}	
		
	
	int wpbobbi; 
	int psbobbi[E];

	int holdi[N][E], wholdi[N], bumbquei[E], wbuqi, tnkx[27][3];
	
		
		
	

	
	// Declaration of bonding variable ends! 
	
	
	
	
	
	
	
	// Initialization of core variables start!
	
	
	adefus=0.0;
	adefusl=0.0;
	for (i=0; i<N; i++)
	{
		
		
		ldcxyz[i][0]=0.0;
		ldcxyz[i][1]=0.0;
		ldcxyz[i][2]=0.0;
		ldcxy[i][0]=0.0;
		ldcxy[i][1]=0.0;
		ldcz[i]=0.0;
		lom[i][0]=0.0;
		lom[i][1]=0.0;
		lom[i][2]=0.0;
		
		
		wnrd[i]=0;
		for (j=0; j<UN; j++)
		{
			nrdi[i][j]=6000;
		
		}			
		
	}
	
	
	
	
	l=-1;
	jt=0;
	for (i=0; i<3; i++)
	{
		m=-1;
		for (j=0; j<3; j++)
		{
			n=-1;
			for (k=0; k<3; k++)
			{
				tnkx[jt][0]=l;
				tnkx[jt][1]=m;
				tnkx[jt][2]=n;
				n=n+1;
				jt=jt+1;
			}
			m=m+1;
		
		}
		l=l+1;
	}
	
	
	
		
	
	
	
	for (i=0; i<N; i++)
	{
		jt=0;
		jtt=0;
		for (j=0; j<war; j++)
		{
			jtt=jtt+(int)round(per[j]*(float)N);
			if (i>=jtt)    jt=jt+1;		
		}
		pid[i]=jt;
		
		
		//printf("%d  %d\n", i, pid[i]);
	}
	
	
	
	
	double* rdtheta=malloc(war*sizeof(double));
	double* grxy=malloc(war*sizeof(double));
	double* grz=malloc(war*sizeof(double));
	double* grth=malloc(war*sizeof(double));
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Dynamical parameters are given as- \n");
	
	for (jt=0; jt<war; jt++)
	{
		if(rm[jt]<r[jt])
		{
			asp=ar[jt];
			if (asp==1.0)
			{
				
				grz[jt]=1.0;
				grxy[jt]=1.0;
				grth[jt]=1.0;
			
			}
			else
			{
				q1=asp+pow((pow(asp,2.0)-1.0),0.5);
				q2=asp-pow((pow(asp,2.0)-1.0),0.5);
				q3=pow(asp,2.0)-1.0;
				q4=pow(q3,(3.0/2.0));
				
				q5=1.0-(asp*asp);
				q6=2.0*(asp*asp);
				q7=q6-1.0;
				grz[jt]=(8.0/3.0)*(1.0/(((2.0*asp)/q5)+((q7/q4)*log(q1/q2))));
				grxy[jt]=(8.0/3.0)*(1.0/((asp/q3)+(((q6-3.0)/q4)*log(q1))));
				grth[jt]=(2.0/3.0)*((pow(asp,4.0)-1.0)/(asp*(((q7/pow(q3,0.5))*log(q1))-asp)));
				
				/*q11=rm[jt]/0.5;
				grz[jt]=pow((1.0/(grz[jt]*q11)),0.5);
				grxy[jt]=pow((1.0/(grxy[jt]*q11)),0.5);
				grth[jt]=pow((1.0/grth[jt]),0.5);*/
			}
			
			q11=rm[jt]/0.5;
			q22=(rm[jt]*rm[jt]*r[jt])/(0.5*0.5*0.5);	
			grz[jt]=pow((1.0/(grz[jt]*q11)),0.5);
			grxy[jt]=pow((1.0/(grxy[jt]*q11)),0.5);
			grth[jt]=pow((1.0/(grth[jt]*q22)),0.5);
			
		}
		else
		{
			asp=ar[jt];
			if (asp==1.0)
			{
				grz[jt]=1.0;
				grxy[jt]=1.0;
				grth[jt]=1.0;
			
			}
			else
			{
				q1=1.0-(2.0*asp*asp);
				q2=1.0-(asp*asp);
				q3=pow(q2,0.5);
				q4=3.0-(2.0*asp*asp);
				q5=-q1;
				q6=(1.0/q3)*acos(asp);
				grz[jt]=(6.0/(8.0*q2))*(((q1*acos(asp))/q3)+asp);
				grxy[jt]=(3.0/(8.0*q2))*(((q4*acos(asp))/q3)-asp);
				grth[jt]=(3.0/2.0)*((((q5*q6)-asp)*asp)/(pow(asp,4.0)-1.0));
				
				/*q11=rm[jt]/0.5;
				grz[jt]=pow((grz[jt]/q11),0.5);
				grxy[jt]=pow((grxy[jt]/q11),0.5);
				grth[jt]=pow(grth[jt],0.5);*/
			}
			q22=(rm[jt]*rm[jt]*r[jt])/(0.5*0.5*0.5);			
			q11=rm[jt]/0.5;
			grz[jt]=pow((grz[jt]/q11),0.5);
			grxy[jt]=pow((grxy[jt]/q11),0.5);
			grth[jt]=pow((grth[jt]/q22),0.5);
			
		}
		printf("For family of particles %d:     grz= %lf,  grxy=%lf, and  grth=%lf\n", jt, grz[jt], grxy[jt], grth[jt]); 
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	for (jt=0; jt<war; jt++)
	{
		
		rdtheta[jt]=(pow(2.0,0.5))*step*grth[jt];
		printf("For family of particles %d:     S_xy = %lf,  S_z = %lf and  S_r = %lf\n", jt, grxy[jt]*step, grz[jt]*step, rdtheta[jt]);
		
	}
	
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	// Initialization of core variables ends!
	
	
	
	
	
	
	
	
	
	// Initialization of bond variables starts!
	
	
	for (i=0; i<N; i++)
	{
		wbobb[i]=0;
		wcbobb[i]=0;
		
		wbobbi[i]=0;
		wcbobbi[i]=0;
		for (j=0; j<E; j++)
		{
			bobb[i][j]=6000;
			cbobb[i][j]=6000;
			
			bobbi[i][j]=6000;
			cbobbi[i][j]=6000;
			
		}			
	}
		
	wpbobbi=0;
	wpbobb=0;
	for(i=0; i<E; i++)
	{
		psbobb[i]=6000;
		
		psbobbi[i]=6000;
	}
		
	

	
	
	// Initialization of bond variables ends!
	
	
	
	
	
	
	
	// Setting initial coordinates and direction starts!
	
	ix=block[0];
	iy=block[1];
	iz=block[2];
	
	FILE* cfu;
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	hintgrl=1000;
	
	
	double** vali=malloc(hintgrl*sizeof(double*));
	for (i=0; i<hintgrl; i++)
	{
		vali[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<hintgrl; i++)
	{
		vali[i][0]=0.0;
		vali[i][1]=0.0;
	}
	
	FILE*  raj;
	
	raj=fopen("Fitted_curve.txt", "r");
	for (i=0; i<10000; i++)
	{
		fscanf(raj, "%lf", &vali[i][0]);
		fscanf(raj, "%lf", &vali[i][1]);
		
		if ( vali[i][1]==0.0) {printf("File have been read\n");  break;}
		if (vali[i][1]<0.0) vali[i][1]=0.0;
	}

	fclose(raj);
	
	
	
	
	printf("Initial configuration within the simulation box have been started generating: \n \n");
	for(i=0; i<N; i++)	
    	{	
		
		
		ra1=drand48();
        	theta=acos(1.0-(2.0*ra1));
		ra2=drand48();
      	 	phi=2.0*piee*ra2;
		
		
		pdx=sin(theta)*cos(phi);
		pdy=sin(theta)*sin(phi);
		pdz=cos(theta);
		
		li[i][6]=0.0;
		li[i][7]=0.0;
		li[i][8]=1.0;
		
		//printf("number filled %d\n", i);
	}
	
	if (ar[1]>1.0)
	{
		
		ofact=1.0;
	}
	else
	{
		
		ofact=-1.0;
	}
	
		//printf("Go yourself\n");
		
	cow=0;
	sitts=0;
	quao=0.0;
	quaoi=0.0;
	
	char filename[25] = {0};

	sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

	
	cfu=fopen(filename,"w");
	fclose(cfu);
	printf("Initial configuration generated successfully!!!!!!! \n \n");
	
	
	
	
	for (k=0; k<N; k++)
	{
		rpid=pid[k];
		
		
		ami[k][2]=20.0;
		ami[k][0]=20.0;
		ami[k][1]=20.0;
		
		theta=acos(li[k][8]);
		//printf("%lf\n", theta);
		if (fabs(li[k][8]) > 0.99999)
		{
		dx=1.0;
		dy=0.0;
		dz=0.0;
		}
		else
		{
		dx=-li[k][7]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
		dy=li[k][6]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
		dz=0.0;
		}
		q1=cos(theta/2.0);
		q2=dx*sin(theta/2.0);
		q3=dy*sin(theta/2.0);
		q4=dz*sin(theta/2.0);
		q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
		q21=2.0*((q2*q3)+(q4*q1));
		q31=2.0*((q2*q4)-(q3*q1));
		q12=2.0*((q2*q3)-(q4*q1));
		q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
		q32=2.0*((q3*q4)+(q2*q1));
		q13=2.0*((q2*q4)+(q3*q1));
		q23=2.0*((q3*q4)-(q2*q1)); 
		q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));

		li[k][0]=1.0*q11;
		li[k][1]=1.0*q21;
		li[k][2]=1.0*q31;
		li[k][3]=1.0*q12;
		li[k][4]=1.0*q22;
		li[k][5]=1.0*q32;
		li[k][6]=1.0*q13;
		li[k][7]=1.0*q23;
		li[k][8]=1.0*q33;
		
		//printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \n", q11, q21, q31, q12, q22, q32, q13, q23, q33);
	}
	
	
	for (i=0; i<N; i++)
	{
		rpid=pid[i];
		if (rpid>=war) printf("raaam naam\n");
		el[i][0]=(rsm[rpid]*li[i][0]*li[i][0])+(rsm[rpid]*li[i][3]*li[i][3])+(rs[rpid]*li[i][6]*li[i][6]);
		el[i][1]=(rsm[rpid]*li[i][0]*li[i][1])+(rsm[rpid]*li[i][3]*li[i][4])+(rs[rpid]*li[i][6]*li[i][7]);
		el[i][2]=(rsm[rpid]*li[i][0]*li[i][2])+(rsm[rpid]*li[i][3]*li[i][5])+(rs[rpid]*li[i][6]*li[i][8]);
		el[i][3]=(rsm[rpid]*li[i][0]*li[i][1])+(rsm[rpid]*li[i][3]*li[i][4])+(rs[rpid]*li[i][6]*li[i][7]);
		el[i][4]=(rsm[rpid]*li[i][1]*li[i][1])+(rsm[rpid]*li[i][4]*li[i][4])+(rs[rpid]*li[i][7]*li[i][7]);
		el[i][5]=(rsm[rpid]*li[i][1]*li[i][2])+(rsm[rpid]*li[i][4]*li[i][5])+(rs[rpid]*li[i][7]*li[i][8]);
		el[i][6]=(rsm[rpid]*li[i][0]*li[i][2])+(rsm[rpid]*li[i][3]*li[i][5])+(rs[rpid]*li[i][6]*li[i][8]);
		el[i][7]=(rsm[rpid]*li[i][1]*li[i][2])+(rsm[rpid]*li[i][4]*li[i][5])+(rs[rpid]*li[i][7]*li[i][8]);
		el[i][8]=(rsm[rpid]*li[i][2]*li[i][2])+(rsm[rpid]*li[i][5]*li[i][5])+(rs[rpid]*li[i][8]*li[i][8]);	
	}
	
	
	// Setting initial coordinates and direction ends!
	


	
	
     	// Randomization starts here!
     	
     	srand48(time(NULL));
	srand(time(NULL));
	//srand48(6);
	//srand(6);

	static unsigned long long  x=123456789,y=987654321,z=43219876,cc=6543217; 
	unsigned int JKISS()
	{ 
		unsigned long long t;
		x=314527869*x+1234567; 
		y^=y <<5;y^=y>>7;y^=y<<22;
		t = 4294584393ULL*z+cc; cc = t>>32; z= t;
		return x+y+z; 
	}

	// Randomization ends here!






	//Preparing for simulation!
     
     
     	
     
     	//FILE* dwu;
     	//dwu=fopen("Data_cord_dire.txt","w");
     	
     
     
     
     
     
     
     
     	int* trac=malloc(war*sizeof(int));
     	
     	for (i=0; i<war; i++)
     	{
     		trac[i]=0;
     	}
     	
     	int* tl=malloc(war*sizeof(int));
     	
     	for (i=0; i<war; i++)
     	{
     		tl[i]=(int)((float)nt[i]*2.0);
     	}
     
     	
     	
     	
     	
     	
     	
     	
	lsd=(int)floor((1.0/(step*step)));
	iltt=(int)floor((float)tsim/(float)lsd);
	rdd=(float)lsd/10.0;
	kd=0.0;
	istep=1.0/(step*step);
	
	//cow=0;
	titi=0;
	ri=0;
	nlies=2*N;
	teb=f[2]-0.01;
	eb=teb;
	rul=1;
	sla=0;
	ve=0;
	ptr=1;
	veei=1;
	vee=0;
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("t_sim required to achive unit physical time for the given step length-    ");
	printf("%d \n", lsd);
	printf("\n");
	printf("\n");

	printf("Maximum possible physical time with given step length and tsim is given by-   ");
	printf("%d\n",iltt);
	printf("\n");
	printf("\n");

	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("                                       *** Shree Ganeshaay Namah ***  :::\n");
	//printf("aimed edge, f[0] %lf %lf %d %d %d\n", f[0], fe, block[0], blockd, iz);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("      *  ::  Be aware of yourself. Rest I will take care. By the way, I am running for you  ::  *\n");
	printf("\n");
	printf("\n");
	printf("\n");
	
	//Preparations done!
	
	//printf("%d\n", pid[299]);
	
	
	
	
	rpid=pid[N-1];
	
	
	difus=2.0;
	difusl=2.0*difus/(float)(N-3);

	hintgrl=100;
	boss=piee/hintgrl;
	
	j=(int)(1000/hintgrl);
	
	
	
	k=j/2;
	qua=0.0;
	bdm=boss*0.5;
	for(i=0; i<hintgrl; i++)
	{
		
		qua=qua+boss*vali[k][1]*sin(bdm);
		
		bdm=bdm+boss;
		//printf("%lf \n", vali[k][1]);
	
		k=k+j;
	}
	
	
	k=j/2;
	for(i=0; i<hintgrl; i++)
	{
	
		svali[i][0]=(float)i*boss+0.5*boss;
		svali[i][1]=vali[k][1]/(qua*2.0*piee);
		k=k+j;
	}
	
	qua=0.0;
	bdm=0.0;
	
	for (i=0; i<hintgrl; i++)
	{
		bdm=bdm+svali[i][1]*sin(svali[i][0]);
		
		//printf("%lf  %lf \n", svali[i][0],  svali[i][1]);
	}
	
	printf("here is the problem \n\n %lf \n\n", tpiee* piee*  bdm/(float)hintgrl);
	
	
	
	FILE* hulu; 
	hulu=fopen("New_fitt.txt", "w");
	for (i=0; i<hintgrl; i++)
	{
		fprintf(hulu, "%lf   %lf\n", 180*svali[i][0]/piee, svali[i][1]);
	
	}
	fclose(hulu);
	
	qua=0.0;
	for(i=0; i<hintgrl; i++)
	{
		qua=qua+boss*svali[i][1];
	
	}
	
	
	
	//printf("Here is the problem %lf \n", qua);

	nematicd[0]=sqrt(1.0-nematicd[2]*nematicd[2]);
	nematicd[1]=0.0;
	
	
	//Simulation starts here!
	for (i=0; i<N; i++) printf("%d\n", pid[i]);
	
	
	FILE* hu;
	hu=fopen("ZM.txt","w");
	double rdthp[2];
	sitts=hintgrl*0.5;
	
	
	intgrl2=0.0;
	printf("Budhdhha was born to decore the world !!!!!\n");
	
	double sphereadd, nsphereadd;
	sphereadd=0.0;
	nsphereadd=0.0;
	for (i=0 ;i <10; i++) {  inv5[i][0]= tpiee*(float)i/10.0;   inv5[i][1]=0.0; }
	
	
	if (ar[1]>=1.0)
	{
		range_1=rm[0]+rm[1]; 
		range_2=3.0*r[1]-rm[1];
	}
	else
	{
		range_1=rm[0]+r[1]; 
		range_2=3.0*rm[1]-r[1];
	}
	
	
	
		
	for(varho=0; varho<10; varho++)
	{
		
		
		
		
		for (i=0 ;i <11; i++) {  inv4[i][0]= rm[0] + rm[rpid1] + (3.*r[rpid1]-rm[rpid1])*(float)i/10.0;   inv4[i][1]=0.0; }
		
		printf("Budhdhha was born to decore the world !!!!!\n");
		for(vart=0; vart<11; vart++)
		{
			
			rpid1=pid[N-2];
			rdistance=(range_1 + 0.000001+ range_2*(float)vart/10.0);
			
			ami[N-2][1]=20.0+sin(tpiee*((float)varho/10.0))*rdistance;
			
			ami[N-2][0]=20.0+cos(tpiee*((float)varho/10.0))*rdistance;
			
			
			for (i=0 ;i <10; i++) {  inv3[i][0]= piee*(float)i/10.0;   inv3[i][1]=0.0; }
			for(sitt=0; sitt<10; sitt++)
			{	
				
				//printf("Budhdhha was born to decore the world !!!!!\n");
				
				//printf("raam\n");
				
				
				k=N-1;
				
				rpid1=pid[k];
				
				
				
				rdthp[0]=piee*(float)sitt/10.0;
				qua=0.0;
				
				
				if (ar[rpid1]>1.0)  kd=d[rpid1];
				else kd=2.0*rm[rpid1];
				
				
				
				for (i=0 ;i <10; i++) {  inv2[i][0]= tpiee*(float)i/10.0;   inv2[i][1]=0.0; }
				
				for(ittt=0; ittt<10; ittt++)
				{	
				
				
				//	printf("Budhdhha was born to decore the world !!!!!\n");
					k=N-2;
					
					rpid1=pid[k];
				
					
					qua=0.0;
					sene=0.0;
					
					sene2=0.0;
					
					rpid1=pid[k];
					
					
			    		
					
					rdphi=tpiee*(float)ittt/10.0;
					
					
					rdth=rdthp[0];
					
					
					
					xf=sin(rdth)*cos(rdphi);
					yf=sin(rdth)*sin(rdphi);
					zf=cos(rdth);
					
					if (fabs(zf) > 0.99999)
					{
						dx=1.0;
						dy=0.0;
						dz=0.0;
					}
					else
					{
						dx=(-yf)*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
						dy=xf*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
						dz=0.0;
							
						
					}
					
					
					
					
					
					rdom[0]=dx;
					rdom[1]=dy;
					rdom[2]=dz;
					q1=cos(rdth/2.0);
					q2=dx*sin(rdth/2.0);
					q3=dy*sin(rdth/2.0);
					q4=dz*sin(rdth/2.0);
					q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
					q12=2.0*((q2*q3)+(q4*q1));
					q13=2.0*((q2*q4)-(q3*q1));
					q21=2.0*((q2*q3)-(q4*q1));
					q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
					q23=2.0*((q3*q4)+(q2*q1));
					q31=2.0*((q2*q4)+(q3*q1));
					q32=2.0*((q3*q4)-(q2*q1)); 
					q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
					rdli[0]=q11;
					rdli[1]=q12;
					rdli[2]=q13;
					rdli[3]=q21;
					rdli[4]=q22;
					rdli[5]=q23;
					rdli[6]=q31;
					rdli[7]=q32;
					rdli[8]=q33;
					
					
					theta=acos(nematicd[2]);
					if (fabs(nematicd[2]) > 0.99999)
					{
						dx=1.0;
						dy=0.0;
						dz=0.0;
					}
					else
					{
						
						dx=-nematicd[1]*(1.0/pow(((nematicd[0]*nematicd[0])+(nematicd[1]*nematicd[1])),0.5));
						dy=nematicd[0]*(1.0/pow(((nematicd[0]*nematicd[0])+(nematicd[1]*nematicd[1])),0.5));
						dz=0.0;
						
					}
					q1=cos(theta/2.0);
					q2=dx*sin(theta/2.0);
					q3=dy*sin(theta/2.0);
					q4=dz*sin(theta/2.0);
					q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
					q21=2.0*((q2*q3)+(q4*q1));
					q31=2.0*((q2*q4)-(q3*q1));
					q12=2.0*((q2*q3)-(q4*q1));
					q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
					q32=2.0*((q3*q4)+(q2*q1));
					q13=2.0*((q2*q4)+(q3*q1));
					q23=2.0*((q3*q4)-(q2*q1)); 
					q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
					
					li[k][0]=rdli[0]*q11+rdli[1]*q12+rdli[2]*q13;
					li[k][1]=rdli[0]*q21+rdli[1]*q22+rdli[2]*q23;
					li[k][2]=rdli[0]*q31+rdli[1]*q32+rdli[2]*q33;
					li[k][3]=rdli[3]*q11+rdli[4]*q12+rdli[5]*q13;
					li[k][4]=rdli[3]*q21+rdli[4]*q22+rdli[5]*q23;
					li[k][5]=rdli[3]*q31+rdli[4]*q32+rdli[5]*q33;
					li[k][6]=rdli[6]*q11+rdli[7]*q12+rdli[8]*q13;
					li[k][7]=rdli[6]*q21+rdli[7]*q22+rdli[8]*q23;
					li[k][8]=rdli[6]*q31+rdli[7]*q32+rdli[8]*q33;
					
					
					rst=rs[rpid1];
					rsmt=rsm[rpid1];
					
					//printf("%lf   %lf %d\n", rst, rsmt, rpid1);
					el[k][0]=(rsmt*li[k][0]*li[k][0])+(rsmt*li[k][3]*li[k][3])+(rst*li[k][6]*li[k][6]);
					el[k][1]=(rsmt*li[k][0]*li[k][1])+(rsmt*li[k][3]*li[k][4])+(rst*li[k][6]*li[k][7]);
					el[k][2]=(rsmt*li[k][0]*li[k][2])+(rsmt*li[k][3]*li[k][5])+(rst*li[k][6]*li[k][8]);
					el[k][3]=(rsmt*li[k][0]*li[k][1])+(rsmt*li[k][3]*li[k][4])+(rst*li[k][6]*li[k][7]);
					el[k][4]=(rsmt*li[k][1]*li[k][1])+(rsmt*li[k][4]*li[k][4])+(rst*li[k][7]*li[k][7]);
					el[k][5]=(rsmt*li[k][1]*li[k][2])+(rsmt*li[k][4]*li[k][5])+(rst*li[k][7]*li[k][8]);
					el[k][6]=(rsmt*li[k][0]*li[k][2])+(rsmt*li[k][3]*li[k][5])+(rst*li[k][6]*li[k][8]);
					el[k][7]=(rsmt*li[k][1]*li[k][2])+(rsmt*li[k][4]*li[k][5])+(rst*li[k][7]*li[k][8]);
					el[k][8]=(rsmt*li[k][2]*li[k][2])+(rsmt*li[k][5]*li[k][5])+(rst*li[k][8]*li[k][8]);
					
	
					
					k=N-2;
				
					w1=li[k][0];
					w2=li[k][1];
					w3=li[k][2];
					w4=li[k][3];
					w5=li[k][4];
					w6=li[k][5];
					w7=li[k][6];
					w8=li[k][7];
					w9=li[k][8];
					
					del[0]=el[k][0];
					del[1]=el[k][1];
					del[2]=el[k][2];
					del[3]=el[k][3];
					del[4]=el[k][4];
					del[5]=el[k][5];
					del[6]=el[k][6];
					del[7]=el[k][7];
					del[8]=el[k][8];
					
					
					
					
					
					check=0;
	
					
					
					
					mxi=ami[k][0];
		    			myi=ami[k][1];
		    			mzi=ami[k][2];
			
					difus=1.0;
					
					
					
					flag=0;
					
					for (ip=0; ip<N-2; ip++)
					{
						
						puc=ip;
						rpid2=pid[ip];
						
						if (puc<N-2)
						{
							ami[puc][2]=ami[k][2]+difus;
				
						}
						
						
						rpid2=pid[ip];
						
						u1=li[puc][0];
						u2=li[puc][1];
						u3=li[puc][2];
						u4=li[puc][3];
						u5=li[puc][4];
						u6=li[puc][5];
						u7=li[puc][6];
						u8=li[puc][7];
						u9=li[puc][8];
						re1=ami[puc][0]-mxi;
						re2=ami[puc][1]-myi;
						re3=ami[puc][2]-mzi;
						
						
						
						mr=(re1*re1+re2*re2+re3*re3);
						mrp=pow(mr,0.5);
						
					
						
						
						ree[0]=(re1*w1+re2*w2+re3*w3);
						ree[1]=(re1*w4+re2*w5+re3*w6);
						ree[2]=(re1*w7+re2*w8+re3*w9);
						
						rtm[0][0]=(u1*w1+u2*w2+u3*w3);
						rtm[1][0]=(u1*w4+u2*w5+u3*w6);
						rtm[2][0]=(u1*w7+u2*w8+u3*w9);
						rtm[0][1]=(u4*w1+u5*w2+u6*w3);
						rtm[1][1]=(u4*w4+u5*w5+u6*w6);
						rtm[2][1]=(u4*w7+u5*w8+u6*w9);
						rtm[0][2]=(u7*w1+u8*w2+u9*w3);
						rtm[1][2]=(u7*w4+u8*w5+u9*w6);
						rtm[2][2]=(u7*w7+u8*w8+u9*w9);
						
											
						dive=0.0;
						for (i = 0; i < 3; i++)
						{	
							for (j = 0; j < 3; j++) abrtm[i][j] = fabs(rtm[i][j]) + EPSILON;
						}
						tri=0;
						
						for (i = 0; i < 3; i++) 
						{
							ra = e[rpid1][i];
							rb = e[rpid2][0] * abrtm[i][0] + e[rpid2][1] * abrtm[i][1] + e[rpid2][2] * abrtm[i][2];
							if (fabs(ree[i]) > ra + rb) tri=1;
						}
						
						
						for (i = 0; i < 3; i++) 
						{
							ra = e[rpid1][0] * abrtm[0][i] + e[rpid1][1] * abrtm[1][i] + e[rpid1][2] * abrtm[2][i];
							rb = e[rpid2][i];
							if (fabs(ree[0] * rtm[0][i] + ree[1] * rtm[1][i] + ree[2] * rtm[2][i]) > ra + rb) tri=1;
						}
						
						if (tri==1) goto label17;
						
						ra= e[rpid1][1] * abrtm[2][0] + e[rpid1][2] * abrtm[1][0];
						rb= e[rpid2][1] * abrtm[0][2] + e[rpid2][2] * abrtm[0][1];
						if (fabs(ree[2] * rtm[1][0] - ree[1] * rtm[2][0]) > ra + rb) goto label17;
						
						ra= e[rpid1][1] * abrtm[2][1] + e[rpid1][2] * abrtm[1][1];
						rb= e[rpid2][0] * abrtm[0][2] + e[rpid2][2] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][1] - ree[1] * rtm[2][1]) > ra + rb) goto label17;
						
						ra= e[rpid1][1] * abrtm[2][2] + e[rpid1][2] * abrtm[1][2];
						rb= e[rpid2][0] * abrtm[0][1] + e[rpid2][1] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][2] - ree[1] * rtm[2][2]) > ra + rb) goto label17;
						
						ra = e[rpid1][0] * abrtm[2][0] + e[rpid1][2] * abrtm[0][0];
						rb = e[rpid2][1] * abrtm[1][2] + e[rpid2][2] * abrtm[1][1];
						if (fabs(ree[0] * rtm[2][0] - ree[2] * rtm[0][0]) > ra + rb) goto label17;
					
						ra= e[rpid1][0] * abrtm[2][1] + e[rpid1][2] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[1][2] + e[rpid2][2] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][1] - ree[2] * rtm[0][1]) > ra + rb) goto label17;
						
						ra= e[rpid1][0] * abrtm[2][2] + e[rpid1][2] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[1][1] + e[rpid2][1] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][2] - ree[2] * rtm[0][2]) > ra + rb) goto label17;
						
						ra= e[rpid1][0] * abrtm[1][0] + e[rpid1][1] * abrtm[0][0];
						rb= e[rpid2][1] * abrtm[2][2] + e[rpid2][2] * abrtm[2][1];
						if (fabs(ree[1] * rtm[0][0] - ree[0] * rtm[1][0]) > ra + rb) goto label17;
							
						ra= e[rpid1][0] * abrtm[1][1] + e[rpid1][1] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[2][2] + e[rpid2][2] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][1] - ree[0] * rtm[1][1]) > ra + rb) goto label17;
						
						ra= e[rpid1][0] * abrtm[1][2] + e[rpid1][1] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[2][1] + e[rpid2][1] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][2] - ree[0] * rtm[1][2]) > ra + rb) goto label17;
						
						if (tri==0)
						{
								
							
							del2[0]=el[puc][0];
							del2[1]=el[puc][1];
							del2[2]=el[puc][2];
							del2[3]=el[puc][3];
							del2[4]=el[puc][4];
							del2[5]=el[puc][5];
							del2[6]=el[puc][6];
							del2[7]=el[puc][7];
							del2[8]=el[puc][8];
							
							a=0.1;
							b=1.0;
							c=(b+a)/2.0;
							
							trig=1.0;
							
							cfun[0]=0.0;
							cfun[1]=0.0;
							tri=0;
							
							dive=0.0;
							bvs=0;
							bvsi=0;
							divea=0.0;
							divei=0.0;
							
							
							
							
							while (trig>0.0001)
							{
								tri=tri+1;
								
								
								lemm=1.0-c;
								lem=c;
								lemms=lemm*lemm;
								lems=lem*lem;
								wt=0;
								
								for (i=0; i<3; i++)
								{
									for (j=0; j<3; j++)
									{
										lt[i][j]=del[wt]*lemm+del2[wt]*lem;
										wt=wt+1;
									}
								}
								
								wt=0;
								
								for (i=0; i<3; i++)
								{
									for (j=0; j<3; j++)
									{
										lts[i][j]=del[wt]*lemms-del2[wt]*lems;
										wt=wt+1;
									}
								}
								lt[0][3]=re1;
								lt[1][3]=re2;
								lt[2][3]=re3;
								cont=1;
								for (it=0; it<2; it++) 
								{
								    	dummy=0;
								    	if  (lt[it][it]==0.0)
									{
										for (jt=0; jt<3; jt++) 
										{
									    		if (dummy == 0)	
											{
												if (lt[jt][it] != 0.0) 
												{
										    			mt=jt;
										    			dummy=dummy+1;
										    			dump[0]=lt[it][0];
													dump[1]=lt[it][1];
													dump[2]=lt[it][2];
													dump[3]=lt[it][3];
										    			lt[it][0]=lt[mt][0];
													lt[it][1]=lt[mt][1];
													lt[it][2]=lt[mt][2];
													lt[it][3]=lt[mt][3];
										    			lt[mt][0]=dump[0];
													lt[mt][1]=dump[1];
													lt[mt][2]=dump[2];
													lt[mt][3]=dump[3];
												}
											}
										}
										if (dummy==0) continue; 
									}
								    	for (jt=cont; jt<3; jt++)
									{
										dumm=lt[jt][it];
										for (mt=0; mt<4; mt++)	lt[jt][mt]=lt[jt][mt]-((lt[it][mt]/lt[it][it])*dumm);
										
									}
								    	cont=cont+1;
								}
								//NOW CHECK EXISTENCE OF SOLUTION OF EQUATION
								sol[0]=0.0; sol[1]=0.0; sol[2]=0.0; kt=3;
								for (it=2; it>=0; it--)
								{    	sol[it]=lt[it][kt];
								    	for (jt=0; jt<3; jt++)
									{
										if (it==jt) continue;
										sol[it]=sol[it]-(lt[it][jt]*sol[jt]);
									}
								   	sol[it]=sol[it]/lt[it][it];
								}
								ec=0.0;
								for (i=0; i<3; i++)	ec=ec+(lts[i][0]*sol[0]+lts[i][1]*sol[1]+lts[i][2]*sol[2])*sol[i];
								
													
								if (ec>0.0)
								{
									a=c;
									cfun[1]=ec;
								}
			    					else
			    					{
									b=c;
									cfun[2]=ec;
								}
								c=(a+b)/2.0;	
			    					trig=fabs(b-a);
			    					ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
			    					
			    					ep=ec/mr;
			    					cm=mrp-dive;
								ec=ep*cm*cm;
				    					
					    			if (ec>1.0) trig=0.000001;
			    					
				    				//if (ec>1.0) trig=0.000001;	
								
							}
							
							
							
							
							ep=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
							
							
							if (ec<=1.0)
							{
								
								if (ep<1.0)
								{
							
									check=1;
									jet=jet+1;
								}
							
								
				
							}
							
							
						}
						
						label17:
						
						if (check==1)
						{
							
							break;
							
						
						}
						
						difus=difus-1;
					
					}
					
					flagx=0;
					if (check==1)   		flagx=1;	
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					for (i=0 ;i <10; i++) {  inv1[i][0]= piee*(float)i/10.0;   inv1[i][1]=0.0; }
					
					
					
					
					
					
					for(lo=0; lo<10; lo++)
					{	
				
			    			
			    			rdthp[1]=piee*(float)lo/10.0;
			    			
			    			
			    			pf=0;
			    			
			    			//printf("Budhdhha was born to live a great life!!!!!!\n");
			    			
			    			
			    			for(no=0; no<10; no++)
			    			{
			    			
				    			
							
							k=N-1;
							rpid1=pid[k];
							rdphi=tpiee*(float)no/10.0;
							rdth=rdthp[1];
							xf=sin(rdth)*cos(rdphi);
							yf=sin(rdth)*sin(rdphi);
							zf=cos(rdth);
							
							if (fabs(zf) > 0.99999)
							{
								dx=1.0;
								dy=0.0;
								dz=0.0;
							}
							else
							{
								dx=(-yf)*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
								dy=xf*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
								dz=0.0;
									
								
							}
							
							
							rdom[0]=dx;
							rdom[1]=dy;
							rdom[2]=dz;
							q1=cos(rdth/2.0);
							q2=dx*sin(rdth/2.0);
							q3=dy*sin(rdth/2.0);
							q4=dz*sin(rdth/2.0);
							q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
							q12=2.0*((q2*q3)+(q4*q1));
							q13=2.0*((q2*q4)-(q3*q1));
							q21=2.0*((q2*q3)-(q4*q1));
							q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
							q23=2.0*((q3*q4)+(q2*q1));
							q31=2.0*((q2*q4)+(q3*q1));
							q32=2.0*((q3*q4)-(q2*q1)); 
							q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
							rdli[0]=q11;
							rdli[1]=q12;
							rdli[2]=q13;
							rdli[3]=q21;
							rdli[4]=q22;
							rdli[5]=q23;
							rdli[6]=q31;
							rdli[7]=q32;
							rdli[8]=q33;
							
							
							theta=acos(nematicd[2]);
							if (fabs(nematicd[2]) > 0.99999)
							{
								dx=1.0;
								dy=0.0;
								dz=0.0;
							}
							else
							{
								
								dx=-nematicd[1]*(1.0/pow(((nematicd[0]*nematicd[0])+(nematicd[1]*nematicd[1])),0.5));
								dy=nematicd[0]*(1.0/pow(((nematicd[0]*nematicd[0])+(nematicd[1]*nematicd[1])),0.5));
								dz=0.0;
								
							}
							q1=cos(theta/2.0);
							q2=dx*sin(theta/2.0);
							q3=dy*sin(theta/2.0);
							q4=dz*sin(theta/2.0);
							q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
							q21=2.0*((q2*q3)+(q4*q1));
							q31=2.0*((q2*q4)-(q3*q1));
							q12=2.0*((q2*q3)-(q4*q1));
							q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
							q32=2.0*((q3*q4)+(q2*q1));
							q13=2.0*((q2*q4)+(q3*q1));
							q23=2.0*((q3*q4)-(q2*q1)); 
							q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
							
							li[k][0]=rdli[0]*q11+rdli[1]*q12+rdli[2]*q13;
							li[k][1]=rdli[0]*q21+rdli[1]*q22+rdli[2]*q23;
							li[k][2]=rdli[0]*q31+rdli[1]*q32+rdli[2]*q33;
							li[k][3]=rdli[3]*q11+rdli[4]*q12+rdli[5]*q13;
							li[k][4]=rdli[3]*q21+rdli[4]*q22+rdli[5]*q23;
							li[k][5]=rdli[3]*q31+rdli[4]*q32+rdli[5]*q33;
							li[k][6]=rdli[6]*q11+rdli[7]*q12+rdli[8]*q13;
							li[k][7]=rdli[6]*q21+rdli[7]*q22+rdli[8]*q23;
							li[k][8]=rdli[6]*q31+rdli[7]*q32+rdli[8]*q33;
							
							
							rst=rs[rpid1];
							rsmt=rsm[rpid1];
							
						//	printf("%lf   %lf %d\n", rst, rsmt, rpid1);
							el[k][0]=(rsmt*li[k][0]*li[k][0])+(rsmt*li[k][3]*li[k][3])+(rst*li[k][6]*li[k][6]);
							el[k][1]=(rsmt*li[k][0]*li[k][1])+(rsmt*li[k][3]*li[k][4])+(rst*li[k][6]*li[k][7]);
							el[k][2]=(rsmt*li[k][0]*li[k][2])+(rsmt*li[k][3]*li[k][5])+(rst*li[k][6]*li[k][8]);
							el[k][3]=(rsmt*li[k][0]*li[k][1])+(rsmt*li[k][3]*li[k][4])+(rst*li[k][6]*li[k][7]);
							el[k][4]=(rsmt*li[k][1]*li[k][1])+(rsmt*li[k][4]*li[k][4])+(rst*li[k][7]*li[k][7]);
							el[k][5]=(rsmt*li[k][1]*li[k][2])+(rsmt*li[k][4]*li[k][5])+(rst*li[k][7]*li[k][8]);
							el[k][6]=(rsmt*li[k][0]*li[k][2])+(rsmt*li[k][3]*li[k][5])+(rst*li[k][6]*li[k][8]);
							el[k][7]=(rsmt*li[k][1]*li[k][2])+(rsmt*li[k][4]*li[k][5])+(rst*li[k][7]*li[k][8]);
							el[k][8]=(rsmt*li[k][2]*li[k][2])+(rsmt*li[k][5]*li[k][5])+(rst*li[k][8]*li[k][8]);
							
							
							k=N-1;
							ene=0.0;
							nene=0.0;
							tnene=0.0;
							uti=50;
							
							rpid1=pid[k];
							
							
							
						
							for (mo=0; mo<1000; mo++)
							{
							
								//printf("raam %d\n", mo);
								w1=li[k][0];
								w2=li[k][1];
								w3=li[k][2];
								w4=li[k][3];
								w5=li[k][4];
								w6=li[k][5];
								w7=li[k][6];
								w8=li[k][7];
								w9=li[k][8];
								
								del[0]=el[k][0];
								del[1]=el[k][1];
								del[2]=el[k][2];
								del[3]=el[k][3];
								del[4]=el[k][4];
								del[5]=el[k][5];
								del[6]=el[k][6];
								del[7]=el[k][7];
								del[8]=el[k][8];
								
								
								
								
								
								check=0;
				
								ra1 = JKISS() / 4294967296.0;
								rdth = acos(1.0-(2.0*ra1));
								ra2 = JKISS() / 4294967296.0;
							       	phi = 2.0*piee*ra2;
								di = kd*pow((drand48()), 1.0/3.0);
								dxi = (di * sin(rdth) * cos(phi))+ami[N-2][0];
								dyi = (di * sin(rdth) * sin(phi))+ami[N-2][1];
								dzi = (di * cos(rdth))+ami[N-2][2]; 
								
								ami[k][0]=dxi;
								ami[k][1]=dyi;
								ami[k][2]=dzi; 
								
								
								mxi=ami[k][0];
					    			myi=ami[k][1];
					    			mzi=ami[k][2];
						
								difus=1.0;
								
								
								
								flag=0;
								
								for (ip=0; ip<N-1; ip++)
								{
									
									check=0;
									puc=ip;
									rpid2=pid[ip];
									
									if (puc<N-2)
									{
										ami[puc][2]=ami[k][2]+difus;
							
									}
									
									
									rpid2=pid[ip];
									
									u1=li[puc][0];
									u2=li[puc][1];
									u3=li[puc][2];
									u4=li[puc][3];
									u5=li[puc][4];
									u6=li[puc][5];
									u7=li[puc][6];
									u8=li[puc][7];
									u9=li[puc][8];
									re1=ami[puc][0]-mxi;
									re2=ami[puc][1]-myi;
									re3=ami[puc][2]-mzi;
									
									
									
									mr=(re1*re1+re2*re2+re3*re3);
									mrp=pow(mr,0.5);
									
								
									
									
									ree[0]=(re1*w1+re2*w2+re3*w3);
									ree[1]=(re1*w4+re2*w5+re3*w6);
									ree[2]=(re1*w7+re2*w8+re3*w9);
									
									rtm[0][0]=(u1*w1+u2*w2+u3*w3);
									rtm[1][0]=(u1*w4+u2*w5+u3*w6);
									rtm[2][0]=(u1*w7+u2*w8+u3*w9);
									rtm[0][1]=(u4*w1+u5*w2+u6*w3);
									rtm[1][1]=(u4*w4+u5*w5+u6*w6);
									rtm[2][1]=(u4*w7+u5*w8+u6*w9);
									rtm[0][2]=(u7*w1+u8*w2+u9*w3);
									rtm[1][2]=(u7*w4+u8*w5+u9*w6);
									rtm[2][2]=(u7*w7+u8*w8+u9*w9);
									
														
									dive=0.0;
									for (i = 0; i < 3; i++)
									{	
										for (j = 0; j < 3; j++) abrtm[i][j] = fabs(rtm[i][j]) + EPSILON;
									}
									tri=0;
									
									for (i = 0; i < 3; i++) 
									{
										ra = e[rpid1][i];
										rb = e[rpid2][0] * abrtm[i][0] + e[rpid2][1] * abrtm[i][1] + e[rpid2][2] * abrtm[i][2];
										if (fabs(ree[i]) > ra + rb) tri=1;
									}
									
									
									for (i = 0; i < 3; i++) 
									{
										ra = e[rpid1][0] * abrtm[0][i] + e[rpid1][1] * abrtm[1][i] + e[rpid1][2] * abrtm[2][i];
										rb = e[rpid2][i];
										if (fabs(ree[0] * rtm[0][i] + ree[1] * rtm[1][i] + ree[2] * rtm[2][i]) > ra + rb) tri=1;
									}
									
									if (tri==1) goto label16;
									
									ra= e[rpid1][1] * abrtm[2][0] + e[rpid1][2] * abrtm[1][0];
									rb= e[rpid2][1] * abrtm[0][2] + e[rpid2][2] * abrtm[0][1];
									if (fabs(ree[2] * rtm[1][0] - ree[1] * rtm[2][0]) > ra + rb) goto label16;
									
									ra= e[rpid1][1] * abrtm[2][1] + e[rpid1][2] * abrtm[1][1];
									rb= e[rpid2][0] * abrtm[0][2] + e[rpid2][2] * abrtm[0][0];
									if (fabs(ree[2] * rtm[1][1] - ree[1] * rtm[2][1]) > ra + rb) goto label16;
									
									ra= e[rpid1][1] * abrtm[2][2] + e[rpid1][2] * abrtm[1][2];
									rb= e[rpid2][0] * abrtm[0][1] + e[rpid2][1] * abrtm[0][0];
									if (fabs(ree[2] * rtm[1][2] - ree[1] * rtm[2][2]) > ra + rb) goto label16;
									
									ra = e[rpid1][0] * abrtm[2][0] + e[rpid1][2] * abrtm[0][0];
									rb = e[rpid2][1] * abrtm[1][2] + e[rpid2][2] * abrtm[1][1];
									if (fabs(ree[0] * rtm[2][0] - ree[2] * rtm[0][0]) > ra + rb) goto label16;
								
									ra= e[rpid1][0] * abrtm[2][1] + e[rpid1][2] * abrtm[0][1];
									rb= e[rpid2][0] * abrtm[1][2] + e[rpid2][2] * abrtm[1][0];
									if (fabs(ree[0] * rtm[2][1] - ree[2] * rtm[0][1]) > ra + rb) goto label16;
									
									ra= e[rpid1][0] * abrtm[2][2] + e[rpid1][2] * abrtm[0][2];
									rb= e[rpid2][0] * abrtm[1][1] + e[rpid2][1] * abrtm[1][0];
									if (fabs(ree[0] * rtm[2][2] - ree[2] * rtm[0][2]) > ra + rb) goto label16;
									
									ra= e[rpid1][0] * abrtm[1][0] + e[rpid1][1] * abrtm[0][0];
									rb= e[rpid2][1] * abrtm[2][2] + e[rpid2][2] * abrtm[2][1];
									if (fabs(ree[1] * rtm[0][0] - ree[0] * rtm[1][0]) > ra + rb) goto label16;
										
									ra= e[rpid1][0] * abrtm[1][1] + e[rpid1][1] * abrtm[0][1];
									rb= e[rpid2][0] * abrtm[2][2] + e[rpid2][2] * abrtm[2][0];
									if (fabs(ree[1] * rtm[0][1] - ree[0] * rtm[1][1]) > ra + rb) goto label16;
									
									ra= e[rpid1][0] * abrtm[1][2] + e[rpid1][1] * abrtm[0][2];
									rb= e[rpid2][0] * abrtm[2][1] + e[rpid2][1] * abrtm[2][0];
									if (fabs(ree[1] * rtm[0][2] - ree[0] * rtm[1][2]) > ra + rb) goto label16;
									
									if (tri==0)
									{
											
										
										del2[0]=el[puc][0];
										del2[1]=el[puc][1];
										del2[2]=el[puc][2];
										del2[3]=el[puc][3];
										del2[4]=el[puc][4];
										del2[5]=el[puc][5];
										del2[6]=el[puc][6];
										del2[7]=el[puc][7];
										del2[8]=el[puc][8];
										
										a=0.1;
										b=1.0;
										c=(b+a)/2.0;
										
										trig=1.0;
										
										cfun[0]=0.0;
										cfun[1]=0.0;
										tri=0;
										
										dive=0.0;
										bvs=0;
										bvsi=0;
										divea=0.0;
										divei=0.0;
										
										
										
										
										while (trig>0.0001)
										{
											tri=tri+1;
											
											
											lemm=1.0-c;
											lem=c;
											lemms=lemm*lemm;
											lems=lem*lem;
											wt=0;
											
											for (i=0; i<3; i++)
											{
												for (j=0; j<3; j++)
												{
													lt[i][j]=del[wt]*lemm+del2[wt]*lem;
													wt=wt+1;
												}
											}
											
											wt=0;
											
											for (i=0; i<3; i++)
											{
												for (j=0; j<3; j++)
												{
													lts[i][j]=del[wt]*lemms-del2[wt]*lems;
													wt=wt+1;
												}
											}
											lt[0][3]=re1;
											lt[1][3]=re2;
											lt[2][3]=re3;
											cont=1;
											for (it=0; it<2; it++) 
											{
											    	dummy=0;
											    	if  (lt[it][it]==0.0)
												{
													for (jt=0; jt<3; jt++) 
													{
												    		if (dummy == 0)	
														{
															if (lt[jt][it] != 0.0) 
															{
													    			mt=jt;
													    			dummy=dummy+1;
													    			dump[0]=lt[it][0];
																dump[1]=lt[it][1];
																dump[2]=lt[it][2];
																dump[3]=lt[it][3];
													    			lt[it][0]=lt[mt][0];
																lt[it][1]=lt[mt][1];
																lt[it][2]=lt[mt][2];
																lt[it][3]=lt[mt][3];
													    			lt[mt][0]=dump[0];
																lt[mt][1]=dump[1];
																lt[mt][2]=dump[2];
																lt[mt][3]=dump[3];
															}
														}
													}
													if (dummy==0) continue; 
												}
											    	for (jt=cont; jt<3; jt++)
												{
													dumm=lt[jt][it];
													for (mt=0; mt<4; mt++)	lt[jt][mt]=lt[jt][mt]-((lt[it][mt]/lt[it][it])*dumm);
													
												}
											    	cont=cont+1;
											}
											//NOW CHECK EXISTENCE OF SOLUTION OF EQUATION
											sol[0]=0.0; sol[1]=0.0; sol[2]=0.0; kt=3;
											for (it=2; it>=0; it--)
											{    	sol[it]=lt[it][kt];
											    	for (jt=0; jt<3; jt++)
												{
													if (it==jt) continue;
													sol[it]=sol[it]-(lt[it][jt]*sol[jt]);
												}
											   	sol[it]=sol[it]/lt[it][it];
											}
											ec=0.0;
											for (i=0; i<3; i++)	ec=ec+(lts[i][0]*sol[0]+lts[i][1]*sol[1]+lts[i][2]*sol[2])*sol[i];
											
																
											if (ec>0.0)
											{
												a=c;
												cfun[1]=ec;
											}
						    					else
						    					{
												b=c;
												cfun[2]=ec;
											}
											c=(a+b)/2.0;	
						    					trig=fabs(b-a);
						    					ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
						    					
						    					ep=ec/mr;
						    					cm=mrp-dive;
											ec=ep*cm*cm;
							    					
								    			if (ec>1.0) trig=0.000001;
						    					
							    				//if (ec>1.0) trig=0.000001;	
											
										}
										
										
										
										
										ep=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
										
										
										if (ec<=1.0)
										{
											
											if (ep<1.0)
											{
										
												check=1;
												jet=jet+1;
											}
										
											
							
										}
										
										
									}
									
									label16:
									
									if (check==1)
									{
										
										qua=qua+ami[k][2];
										
										
										if(ip<N-2)
										{
											flag=1;
										}
										else
										{
											ene2=ene2+1;
										
										}

										if (ip==N-2 ) 
										{   
											
											if (flag==0 && flagx==0)
											{	
												ene=ene+1; 
												
												dx=mxi-ami[0][0];
												dy=myi-ami[0][1];
												dx=dx/sqrt(dx*dx+dy*dy);
												bdm=(float)varho*tpiee/10.0;
												
												dx=acos(dx);
												if (dy<0.0) dx=tpiee-dx; 
												
												if (bdm-0.1<dx && dx<bdm+0.1)
												{
													bdm=sqrt((ami[0][0]-mxi)*(ami[0][0]-mxi)+(ami[0][1]-myi)*(ami[0][1]-myi));
													if (uti>bdm) uti=bdm;
												}	
											}
											
											nene=nene+1;
	
										}
										
										
									
									}
									else 
									{
										if (ip==N-2 ) 
										{   
											
											if (flag==0)
											{	
												dx=mxi-ami[0][0];
												dy=myi-ami[0][1];
												dx=dx/sqrt(dx*dx+dy*dy);
												bdm=(float)varho*tpiee/10.0;
												
												dx=acos(dx);
												if (dy<0.0) dx= tpiee-dx; 
												
												if (bdm-0.1 < dx && dx < bdm+0.1)
												{
													bdm=sqrt((ami[0][0]-mxi)*(ami[0][0]-mxi)+(ami[0][1]-myi)*(ami[0][1]-myi));
													if (uti>bdm) uti=bdm; 
												}	
											}
											
										}	
									
									}
									difus=difus-1;
								
								}
								
								
								
								
								//printf("raam %d\n", mo);
							}
							
							ene=(ene-nene)*rdistance/1000.0;
							vol=(4.0/3.0)*piee*kd*kd*kd;
						
						
							inv[no][1]=ene*vol;
							
							hyna1[no][1]=0.5*(uti*uti-rm[0]*rm[0]);
							inv[no][0]=tpiee*(float)no/10.0;
							hyna1[no][0]=inv[no][0];
							
							spider1[no][0]=tpiee*(float)no/10.0;
							
							spider1[no][1]=nene*vol/1000.0;
							
							//if (fabs(li[k][8]) <0.1)
							//printf("Here is the energy value:   %lf %lf %lf %lf \n", uti, li[N-1][6], li[N-1][7],li[N-1][8]);
				
						}
						
						
						liml=0.0;
						limu=tpiee;
						iflag=0;
						subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
						
					///	printf("%lf \n", fintgrl);
						
						inv1[lo][0]=rdthp[1];
						inv1[lo][1]=fintgrl;	
						
						if (vart==0 && ittt==0 && sitt==0)
						{
							for (i=0 ;i <10; i++) {  inv[i][0]=hyna1[i][0];   inv[i][1]=hyna1[i][1]; }
							//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
							
							subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
							
							hyna2[lo][0]=rdthp[1];
							hyna2[lo][1]=fintgrl;
							
							
							/*FILE * chuu;
							char filename44[150]={0};
							sprintf(filename44, "Check_phi_%d.txt", lo);
							chuu=fopen(filename44, "w");
							for (i=0; i<10; i++)
							{
								fprintf(chuu, "%lf %lf \n", hyna1[i][0], hyna1[i][1]);
							
							}
							exit(0);
							fclose(chuu);*/
						
						}
						iflag=0;
						
						for (i=0 ;i <10; i++) {  inv[i][0]=spider1[i][0];   inv[i][1]=spider1[i][1]; }
						liml=0.0;
						limu=tpiee;
						iflag=0;
						
						//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
						subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
					///	printf("%lf \n", fintgrl);
						
						spider2[lo][0]=rdthp[1];
						spider2[lo][1]=fintgrl;	
								
						
					}	
					
					for (i=0 ;i <10; i++) {  inv[i][0]=inv1[i][0];   inv[i][1]=inv1[i][1]; }
					
					liml=0.0;
					limu=piee;
					iflag=1;
					//subroutine_integral(liml, limu, inv[11][2], fintgrl,  svali[100][2],  iflag);
					subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
					inv2[ittt][0]=tpiee*(float)ittt/10.0;;
					inv2[ittt][1]=fintgrl;	
					
					if (vart==0 && ittt==0 && sitt==0)
					{
						for (i=0 ;i <10; i++) {  inv[i][0]=hyna2[i][0];   inv[i][1]=hyna2[i][1]; }
						//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
						subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
						hyna3[varho][0]=(float)varho*tpiee/10.0;
						hyna3[varho][1]=fintgrl;
						
						printf("Here is your value whic: %lf  %lf \n", hyna3[varho][0], hyna3[varho][1]);
						
						//exit(0);
						iflag=3;
						//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
						subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
						viper1[varho][0]=(float)varho*tpiee/10.0;
						viper1[varho][1]=fintgrl;
						
					//	exit(0);
						/*FILE * chu;
						char filename4[150]={0};
						sprintf(filename4, "Check_%d.txt", varho);
						chu=fopen(filename4, "w");
						for (i=0; i<10; i++)
						{
							fprintf(chu, "%lf %lf \n", hyna2[i][0], hyna2[i][1]);
						
						}
						
						fclose(chu);*/
						
					}
					
					for (i=0 ;i <10; i++) {  inv[i][0]=spider2[i][0];   inv[i][1]=spider2[i][1]; }
					liml=0.0;
					limu=piee;
					iflag=1;
					
					//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
					subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
				///	printf("%lf \n", fintgrl);
					
					spider3[ittt][0]=tpiee*(float)ittt/10.0;;
					spider3[ittt][1]=fintgrl;	
					
						
					
				}
				
				for (i=0 ;i <10; i++) {  inv[i][0]=inv2[i][0];   inv[i][1]=inv2[i][1]; }
						
				liml=0.0;
				limu=tpiee;
				iflag=0;
				
				//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
				subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
				//printf("%lf \n", fintgrl);
				inv3[sitt][0]=rdthp[0];
				inv3[sitt][1]=fintgrl;		
				
				printf("Simulation is going on si %d %lf\n", sitt, qua/100.0);
				fprintf(hu, "%lf %lf \n", rdth, qua/100.0);
				
				for (i=0 ;i <10; i++) {  inv[i][0]=spider3[i][0];   inv[i][1]=spider3[i][1]; }
				liml=0.0;
				limu=tpiee;
				iflag=0;
				
				//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
				
				subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
			///	printf("%lf \n", fintgrl);
				
				spider4[sitt][0]=rdthp[0];
				spider4[sitt][1]=fintgrl;	
			
			}
			for (i=0 ;i <10; i++) {  inv[i][0]=inv3[i][0];   inv[i][1]=inv3[i][1]; }
			liml=0.0;
			limu=piee;
			iflag=1;
			//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
			
			subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
			printf("%lf \n", fintgrl);
			
			inv4[vart][0]=(range_1 + range_2*(float)vart/10.0);
			inv4[vart][1]=fintgrl;	
			
			for (i=0 ;i <10; i++) {  inv[i][0]=spider4[i][0];   inv[i][1]=spider4[i][1]; }
			liml=0.0;
			limu=piee;
			iflag=1;
			//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
			
			subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);	
			///	
			
			printf("Sphere_distribution_calvirial: %lf \n", fintgrl);
			sphereadd=sphereadd+fintgrl;
			nsphereadd=nsphereadd+1;
			
			
		
		}
		
		for (i=0 ;i <11; i++) {  inv[i][0]=inv4[i][0];   inv[i][1]=inv4[i][1]; }
		
		liml=range_1;
		limu=range_1+range_2;
		iflag=2;
		//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
		
		subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
		printf(" Outer most core_  %lf \n", fintgrl);
		
		inv5[varho][0]=(float)varho*tpiee/10.0;
		inv5[varho][1]=fintgrl;	
		

	}
	
	
	FILE* dudu;
	char filenamedu[150]={0};
	sprintf(filenamedu, "Energy.txt");
	dudu=fopen(filenamedu, "w");
	
	fprintf(dudu, "%lf    ", acos(nematicd[2])*180/piee);
	
	bdm=0.0;
	
	
	
	kd=piee*(range_1*range_1-rm[0]*rm[0]);
	
	kd=kd*sphereadd/nsphereadd;
	
	
	fprintf(dudu, "%lf   ", kd);	
		
	double dmue, gamov, entropy, tenergy;	
	liml=0.0;
	limu=tpiee;
	iflag=0;
	for (i=0 ;i <10; i++) {  inv[i][0]=inv5[i][0];   inv[i][1]=inv5[i][1]; }
	//subroutine_integral(liml, limu, inv[10][2], fintgrl, svali[100][2],  iflag);
	
	subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
	entropy=-0.5*(fintgrl-kd)*fvf*6.0*fvf*6.0/(piee*piee);
	
	printf("S_trans_ext:   %lf\n", entropy);

	fprintf(dudu, "%lf    ", entropy);
	
	dmue=7.0; //    9.0 is the chemical potential
	
	dmue=dmue-log(fvf*6.0/(piee*4.0*piee));
	
	liml=0.0;
	limu=tpiee;
	iflag=0;
	for (i=0 ;i <10; i++) {  inv[i][0]=hyna3[i][0];   inv[i][1]=hyna3[i][1]; }
	//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
	subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
	
	gamov=-fintgrl*fvf*6.0/piee;
	
	
	
	printf("S_trans_id:   %lf\n", gamov);
	fprintf(dudu, "%lf    ", gamov);
	
	entropy=entropy+gamov;
	tenergy=-gamov*dmue;
	
	printf("Delta-mue * gamma:   %lf \n", tenergy);
	
	
	fprintf(dudu, "%lf    ", tenergy);
	
	liml=0.0;
	limu=tpiee;
	iflag=0;
	for (i=0 ;i <10; i++) {  inv[i][0]=viper1[i][0];   inv[i][1]=viper1[i][1]; }
	//subroutine_integral(liml, limu, inv[11][2], fintgrl, svali[100][2], iflag);
	subroutine_integral(liml, limu, inv, &fintgrl, svali, iflag);
	
	printf("S_rot:   %lf\n", fintgrl*fvf*6.0/piee);
	
	
	fprintf(dudu, "%lf    ", fintgrl*fvf*6.0/piee);
	
	
	entropy = entropy + fintgrl*fvf*6.0/piee ;
	
	tenergy=tenergy-entropy;
	
	
	
	fprintf(dudu, "%lf    ", tenergy);
	

	
	printf("Here is your energy calculated:              %lf \n", tenergy);
	
	printf("Here is the value::::::::::::::::");
	
	
	
	
	
	
	fclose(dudu);
	
	
	fclose(hu);
	
	jt=0;
		
	
	
	
	
	
	



	//FINAL CORD Zone starts!

	
	
	//FINAL CORD Zone ends!







	




	
	
	
	

	//Now freedom is achieved!	


	
}



//double subroutine_integral(double, double, double, double, double, int)
//double subroutine_integral(double liml, double limu, double inv[][2], double fintgrl, double svali[][2], int iflag)

double subroutine_integral(double liml, double limu, double inv[][2], double* fintgrl, double svali[][2], int iflag)
{


	int  sri, srj;
	
	
	
	
	double srqua;
	int ptr =15, srk;
	double tinv[500][2];
	
	double ttinv[15][2];
	
	ttinv[0][0]=inv[8][0]-limu;
	ttinv[1][0]=inv[9][0]-limu;
	
	ttinv[0][1]=inv[8][1];
	ttinv[1][1]=inv[9][1];
	
	ttinv[12][0]=inv[0][0]+limu;
	ttinv[13][0]=inv[1][0]+limu;
	ttinv[14][0]=inv[2][0]+limu;
	
	ttinv[12][1]=inv[0][1];
	ttinv[13][1]=inv[1][1];
	ttinv[14][1]=inv[2][1];
	for(i=0; i<10 ; i++)
	{
		ttinv[i+2][0]=inv[i][0];
		ttinv[i+2][1]=inv[i][1];
	}
	
	
	
	double slope;
	
	
	bdm=((4.0*inv[1][0]+limu-liml)/500.0 );
	int srnt=20;
	
	if (iflag==2)
	{
		
		slope=(inv[10][1]-inv[9][1])/(inv[10][0]-inv[9][0]);
		
		
		
		ttinv[13][0]=inv[10][0]+(inv[10][0]-inv[9][0]);
		ttinv[14][0]=inv[10][0]+(inv[10][0]-inv[9][0])+(inv[10][0]-inv[9][0]);
		

		ttinv[13][1]=inv[10][1]+slope*(inv[10][0]-inv[9][0]);
		ttinv[14][1]=inv[10][1]+2.0*slope*(inv[10][0]-inv[9][0]);
		
		slope=(inv[1][1]-inv[0][1])/(inv[1][0]-inv[0][0]);
		
		
		ttinv[0][0]=inv[0][0]-2.0*(inv[1][0]-inv[0][0]);
		ttinv[1][0]=inv[0][0]-(inv[1][0]-inv[0][0]);
		
		ttinv[0][1]=inv[0][1]-2.0*(inv[1][0]-inv[0][0])*slope;
		ttinv[1][1]=inv[0][1]-(inv[1][0]-inv[0][0])*slope;
		
		for(i=0; i<11 ; i++)
		{
			ttinv[i+2][0]=inv[i][0];
			ttinv[i+2][1]=inv[i][1];
		}
		
		bdm=((3.0*(inv[1][0]-inv[0][0])+limu-liml)/500.0 );
		/*FILE * juju;
		juju=fopen("Gardlood_2.txt", "w");
		for (i=0 ; i < 15; i++) fprintf(juju, "%lf    %lf   %lf \n", ttinv[i][0],ttinv[i][1], slope ) ;
		
		fclose(juju);*/
		 srnt=2;
	
	}
	
	
	
	
	
	
	
	
	double* srsol=malloc((srnt+1)*sizeof(double));
	//printf("raaam1\n");
	
	double* srdump = malloc((srnt+2)*sizeof(double));
	//printf("raaam1\n");
	
	
	double* srbump=malloc((srnt+(2*srnt)+2)*sizeof(double));
	int srwbump;
	srwbump=0;
	
	
	double* srtrump=malloc((srnt+(2*srnt)+2)*sizeof(double));
	int srwtrump;
	srwbump=0;
	
	//printf("raaam3\n");
	for( sri=0; sri<srnt+(2*srnt)+2; sri++)
	{
		srtrump[sri]=0;
	}
	
	
	
	
	
	
	for (i=1; i<ptr; i++)
	{
		
		/*
		if (i==ptr)
		{
			ttinv[i][0]=limu;
			ttinv[i][1]=inv[0][1];
		}
		else
		{
			ttinv[i][0]=inv[i][0];
			ttinv[i][1]=inv[i][1];
		}
		*/
		
		slope=(ttinv[i][1]- ttinv[i-1][1])/ (ttinv[i][0]- ttinv[i-1][0]);
		
		ec=ttinv[0][0];
		
		for(j=0; j< 500; j++)
		{
			ec=ec+bdm;
			if (ttinv[i-1][0]< ec && ec < ttinv[i][0])
			{
				tinv[j][0]=ec;
				tinv[j][1]=ttinv[i-1][1]+(ec-ttinv[i-1][0])*slope;
				
				/*if (iflag==2)
				{
					FILE * yuyu;
		
					yuyu=fopen("Gardalood.txt", "a");				
					fprintf(yuyu, "%lf   %lf   %lf\n", tinv[j][0], tinv[j][1], bdm);
					fclose(yuyu);
				}*/
				
			}
			
			
		
		}
		
		
		
	}
	ptr=490;
	
	
	
	
	////printf("Budhdhha was born to decore the world !!!!!\n");
	////printf("raaam3\n");
	//printf("%lf %lf \n", inv[0][0], inv[0][1]);
	for (srj=0; srj< ptr; srj++)
	{
		
		//if (iflag==2)printf("bgbg hhhhhhh %d   %lf %lf \n", srj, tinv[srj][0], tinv[srj][1]);
		
		///inv[srj][0]=(float)srj;
		
		//inv[srj][1]=(float)srj;
		srwbump=0;
		srqua=tinv[srj][1];
		srbump[srwbump]=srqua;
		srwbump=srwbump+1;
		for (sri=0; sri< srnt; sri++)
		{
			srqua=srqua*tinv[srj][0]; 
			srbump[srwbump]=srqua;
			srwbump=srwbump+1;
		}	
	    	srbump[srwbump]=1.0;
	    	srwbump=srwbump+1;
	    	
	    	srqua=tinv[srj][0];
	    	
	    	srbump[srwbump]=srqua;
		srwbump=srwbump+1;
	 
		for (sri=0;   sri<((srnt*2)-1); sri++)
		{
			srqua=srqua*tinv[srj][0];
			srbump[srwbump]=srqua;
			srwbump=srwbump+1;
		}
		
		for(sri=0; sri< srnt+(2*srnt)+2; sri++)
		{
			srtrump[sri]=srtrump[sri]+srbump[sri];
		}
		
		
		
	}
	
	//printf("Budhdhha was born to decore the world 1!!!!!\n");
	
	double** srlt=malloc((srnt+1)*sizeof(double*));
	for (sri=0; sri<srnt+1; sri++ )
	{
		srlt[sri]=malloc((srnt+2)*sizeof(double));
	}
	
	
	
	//printf("Budhdhha was born to decore the world 2!!!!!\n");
	
	for ( srj=0; srj< (srnt+1); srj++)
	{
	
		for ( sri=0; sri< (srnt+1); sri++)
		{
			srk=srj+sri+srnt+1;
			srlt[srj][sri]=srtrump[srk];   
		}
		//printf("%d\n", i);
		srlt[srj][sri]=srtrump[srj];
	}
	
	//printf("Budhdhha was born to decore the world 3 !!!!!\n");
	
	double srdumm;
	     
	int srcont=1, srjtt, srdummy, srmt, srjt, srkt, srit;     
	     
	for ( srit=0; srit<srnt; srit++) 
	{
	    //	printf("Here is the gone fly!!!!!!\n");	
	    	
	    	srdummy=0;
	    	if  (srlt[srit][srit]==0.0)
		{
			for (srjtt=0; srjtt<srnt+1; srjtt++) 
			{
		    		if (srdummy == 0)	
				{
					if (srlt[srjtt][srit] != 0.0) 
					{
			    			srmt=srjt;
			    			srdummy=srdummy+1;
			    			
			    			for(sri=0; sri<srnt+2; sri++)
			    			{
			    				srdump[sri]=srlt[srit][sri];
						
			    				srlt[srit][sri]=srlt[srmt][sri];
						
			    				srlt[srmt][sri]=srdump[sri];
						
			    			}
			    			
					}
				}
			}
			if (srdummy==0) continue; 
		}
	    	for (srjtt=srcont; srjtt<srnt+1; srjtt++)
		{
			srdumm=srlt[srjtt][srit];
			for (srmt=0; srmt<srnt+2; srmt++)
			{
		    		srlt[srjtt][srmt]=srlt[srjtt][srmt]-((srlt[srit][srmt]/srlt[srit][srit])*srdumm);
			}
		}
	    	srcont=srcont+1;
	}
		
	//printf("Budhdhha was born to decore the world 4 !!!!!\n");	
		
	for (sri=0; sri<srnt+1; sri++)
	{
		srsol[sri]=0.0;
	}
	
	srkt=srnt+1, srjt;
	for (srit=srnt; srit>=0; srit--)
	{    	srsol[srit]=srlt[srit][srkt];
	    	for (srjt=0; srjt<srnt+1; srjt++)
		{
			if (srit==srjt) continue;
			srsol[srit]=srsol[srit]-(srlt[srit][srjt]*srsol[srjt]);
		}
	   	srsol[srit]=srsol[srit]/srlt[srit][srit];
	   	
	   	//if (iflag==2) printf("%lf\n", srsol[srit]);
	}
	
	//printf("Budhdhha was born to decore the world 5 !!!!!\n");
	for (sri=0; sri<srnt+1; sri++ )
	{
		free(srlt[sri]);
	}


	double boss ;
	boss=piee/100.0;
	
	
	double gum, bum, yog, nyog;
	bum=(ttinv[14][0]-ttinv[0][0])/500.0;
	
	if (iflag==0 || iflag==2)
	{
		gum=ttinv[0][0];
		yog=0.0;
		nyog=0.0;
		
		for (i=0; i<500; i++)
		{
			
			qua=0.0;
			ene=1.0;
			
			for (j=0; j< srnt+1; j++)
			{
				qua=qua+srsol[j]*ene;
				
				ene=ene*gum;
			
			}
			
			if (gum>liml && gum<limu)
			{
				
				yog=yog+qua;
				nyog=nyog+1.0;
				
				
				
				
			
			}
			gum=gum+bum;
		
		}
		
		yog=yog/nyog*(limu-liml);
	
	}
	else if (iflag==3)
	{
		gum=ttinv[0][0];
		yog=0.0;
		nyog=0.0;
		
		for (i=0; i<500; i++)
		{
			
			qua=0.0;
			ene=1.0;
			
			for (j=0; j< srnt+1; j++)
			{
				qua=qua+srsol[j]*ene;
				
				ene=ene*gum;
			
			}
			
			if (gum>liml && gum<limu)
			{
				sri=(int)(gum/boss);
				if (svali[sri][1]>0.0000001)
				{
					yog=yog+qua*svali[sri][1]*log(4*piee*svali[sri][1])*sin(gum);
				}
				
				//printf("%lf  %lf  %lf  %lf  %lf  %d\n", gum, qua, log(4*piee*svali[sri][1]), limu-liml, svali[sri][1], sri );
				nyog=nyog+1.0;
			
			}
			gum=gum+bum;
		
		}
		
		yog=yog/nyog*(limu-liml);
	}
	else
	{
		gum=ttinv[0][0];
		yog=0.0;
		nyog=0.0;
		
		for (i=0; i<500; i++)
		{
			
			qua=0.0;
			ene=1.0;
			
			for (j=0; j< srnt+1; j++)
			{
				qua=qua+srsol[j]*ene;
				
				ene=ene*gum;
			
			}
			
			if (gum>liml && gum<limu)
			{
				sri=(int)(gum/boss);
			
				yog=yog+qua*svali[sri][1]*sin(gum);
				//yog=yog+qua*(1/(4.0*piee))*sin(gum);
				nyog=nyog+1.0;
				
				//yog=yog+tpiee*svali[sri][1]*sin(gum);
				//printf("%lf  %lf  %lf  %lf  %lf  %d\n", gum, qua, tpiee, limu-liml, svali[sri][1], sri );
			
			}
			gum=gum+bum;
		
		}
		
		yog=yog*(limu-liml)/nyog;
	
	
	
	}
	
	*fintgrl=yog;
	
	
	
	
	
	

	free(srlt);
	free(srsol);
	//printf("raaam1\n");
	
	free(srdump);
	//printf("raaam1\n");
	
	
	free(srbump);
	
	
	
	free(srtrump);
	
	
	

	//printf("Result is given as:    \n");
	
}



