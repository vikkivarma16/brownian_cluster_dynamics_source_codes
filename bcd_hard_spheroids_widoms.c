#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <omp.h>
const int E=20, UN=200, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel;

//void ecfcal(int,double[]);

int  N, trim, check, checkl, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, it, nf,  block[3], blockd, jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct, rul, sla, veei, veec, SHA, ptr, ptrr, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt, gv, bvs, bvsi, count, vfrap, probation, wpi, initia_p, widid, wi, wnlies;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, cp, cph, cpsh, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9], f[3], bl1, bl2, bl3, bl[3], f1, f2, f3, ar, per; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, fe, tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON, eb, teb, tbl, pbc,  dif, difu, difus, difl, difusl, adefus, dive, divea, divei, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm, widoms, twidoms;  



int main(void)
{

	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;

	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	
	fscanf(bu,"%lf",&bdm);
	initia_p=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	vfrap=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	N=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	step=bdm;
	
	
	fscanf(bu,"%lf",&bdm);
	tsim=(int)round(bdm);
	
	
	fscanf(bu,"%lf",&bdm);
	probation=(int)round(bdm);
	//tsim=100000;
	
	
	// Definition of global constants ends       !!!!!!!!!!!!!!!!
	
	ivf=0.004;			
	// Correlation function        !!!!!!!!!!!!!!!!
	
	
	
	
	// Correlation function declaration end         !!!!!!!!!!!!!!!!			
	
	
	
	
  	// Definition of simulation constants starts        !!!!!!!!!!!!!!!!
	
	
	//fvf= 0.02 ;
	fscanf(bu,"%lf",&bdm);
	fvf=bdm;
	//fscanf(bu,"%lf",&fvf);
	
	int war;
	
	//fscanf(bu,"%d",&war);
	fscanf(bu,"%lf",&bdm);
	war=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	widid=(int)round(bdm);
	
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
	
	
	
	
	double epsi;
	epsi=EPSILON;
	

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
	if (vfrap==1)  printf("Volume fraction-    %lf \n", fvf);
	else printf("Single particle diffusivity \n");
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Printing things ends now !!!!!!!
	
	// Declaration of simulation variables starts!
	
	
	// Declaration of simulation variables ends!
	
	
	
	
	// Evaluation of simulation variables starts! 
	
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
	
	fe=0.0;
	for (i=0; i<war; i++)
	{
		fe=fe+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*fvf));
	}
	fe=pow(fe,(1.0/3.0));
	
	f[0]=0.0;
	
	for (i=0; i<war; i++)
	{
		f[0]=f[0]+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*fvf));
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
	
	///printf("Hey there is %lf \n", ble);
	
	
	
	
	
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
	
	int** pit=malloc(UN*sizeof(int*));
	for (i=0; i<UN; i++)
	{
		pit[i]=malloc(2*sizeof(int));
	}
	
	
	
	double** li=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<(N+1);++i)
	{
		li[i]=malloc(9*sizeof(double));
	}
	
	
	double** el=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<(N+1);++i)
	{
		el[i]=malloc(9*sizeof(double));
	}
	
	double** ami=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<(N+1);++i)
	{
		ami[i]=malloc(3*sizeof(double));
	}
	int* pid=malloc((N+1)*sizeof(int));
	
	
	double ldcxyz[N][3], ldcxy[N][2], ldcz[N], lom[N][3];
	
	
	
	int*** nrdi=malloc(N*sizeof(int**));
	int* wnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		nrdi[i]=malloc(UN*sizeof(int*));
		for (j=0; j<UN; j++)
		{
			nrdi[i][j]=malloc(2*sizeof(int));
		}
		
	}
	
	
	
	int**** dri=malloc(blockd*sizeof(int***));
	int*** wdr=malloc(blockd*sizeof(int**));
		
	for (j=0; j<blockd; j++)
	{
		dri[j]=malloc(blockd*sizeof(int**));
		wdr[j]=malloc(blockd*sizeof(int*));
		for (i=0; i<blockd; i++)
		{
			dri[j][i]=malloc(blockd*sizeof(int*));
			wdr[j][i]=malloc(blockd*sizeof(int));
	
			for (k=0; k<blockd; k++)
			{
				dri[j][i][k]=malloc(80*sizeof(int));
			}
		}	
	}
	
	
	
	
	

	int tnkx[27][3];
	
		
		
	

	
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
			nrdi[i][j][0]=6000;
			nrdi[i][j][1]=13;
		
		}			
		
	}
	
	for (i=0; i<blockd; i++)
	{
		
		for (j=0; j<blockd; j++)
		{
			for (k=0; k<blockd; k++)
			{
				wdr[i][j][k]=0;
				for (l=0; l<80; l++)
				{
					dri[i][j][k][l]=6000;
				}
			}	
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
	
	
	pid[N]=widid;
	
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
	
	
	

	
	
	// Initialization of bond variables ends!
	
	
	
	
	
	
	
	// Setting initial coordinates and direction starts!
	
	ix=block[0];
	iy=block[1];
	iz=block[2];
	
	FILE* cfu;
	
	FILE* vu;
	vu=fopen("Relaxed_system.txt","r");
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	if (vu==NULL || initia_p==0)
	{
		///printf("Initial configuration within the simulation box have been started generating: \n \n");
		for(i=0; i<N; i++)	
	    	{	
			//printf("%d\n", i);
			ra1=drand48();
	        	theta=acos(1.0-(2.0*ra1));
			ra2=drand48();
	  	     	phi=2.0*piee*ra2;
			
			pdx=sin(theta)*cos(phi);
			pdy=sin(theta)*sin(phi);
			pdz=cos(theta);
			
			li[i][6]=pdx;
			li[i][7]=pdy;
			li[i][8]=pdz;
				
			xi=(f[0]-0.5)*drand48();
	    		yi=(f[1]-0.5)*drand48();
	    		zi=(f[2]-0.5)*drand48();
	    		xk=(int)floor(xi/bl[0]);
	    		yk=(int)floor(yi/bl[1]);
	    		zk=(int)floor(zi/bl[2]);
	    		///printf("%lf   %lf   %lf   %d   %d   %d\n", bl[0], bl[1], bl[2], xk, yk, zk);
			dum=wdr[xk][yk][zk];
			dri[xk][yk][zk][dum]=i;
			wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
			ami[i][0]=xi;
			ami[i][1]=yi;
			ami[i][2]=zi;
			
			
		}
		///printf("Go yourself\n");
		
		cow=0;
		sitts=0;
		
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		fclose(cfu);
		printf("Initial configuration generated successfully!!!!!!! \n \n");
		
		
		for (k=0; k<N; k++)
		{
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
		
	}
	else
	{
		printf("A pre-prepared system found and intial configuration will be set accordingly. \n");
		printf("Initial configuration is being read from the available file!!! \n \n  ");
		cow=1;
		
		fscanf(vu, "%d", &sitts);
	
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		
		fclose(cfu);
			
		fscanf(vu,"%lf", &f[0]);
		fscanf(vu,"%lf", &f[1]);
		fscanf(vu,"%lf", &f[2]);
		
		for (i=0; i<N; i++)
		{
			
			for(l=0;l<3;l++)
			{
				fscanf(vu,"%lf", &ami[i][l]);
			}
			
			fscanf(vu,"%d", &wnrd[i]);
			for (j=0; j<wnrd[i]; j++) 
			{
				fscanf(vu,"%d", &nrdi[i][j][0]);
				fscanf(vu,"%d", &nrdi[i][j][1]);
			}
			
			
			
			for(l=0;l<9;l++)
			{
				fscanf(vu,"%lf", &li[i][l]);
				
			}
		
		}
		
		
		for ( i=0; i<3; i++)
		{
			
			block[i]=(int)(floor(f[i]/ble));
			bl[i]=f[i]/(float)block[i];
			f[i]=bl[i]*(float)block[i];
			//printf("ptr %d \n", ptr);
		}
		
			
		
		for (i=0; i<blockd; i++)
		{
			
			for (j=0; j<blockd; j++)
			{
				for (k=0; k<blockd; k++)
				{
					wdr[i][j][k]=0;
					for (l=0; l<80; l++)
					{
						dri[i][j][k][l]=6000;
					}
				}	
			}
		}
		
		for(i=0; i<N; i++)	
	    	{	
			
			xi=ami[i][0];
			yi=ami[i][1];
			zi=ami[i][2];
	    		xk=(int)floor(xi/bl[0]);
	    		yk=(int)floor(yi/bl[1]);
	    		zk=(int)floor(zi/bl[2]);
	    		
			dum=wdr[xk][yk][zk];
			dri[xk][yk][zk][dum]=i;
			wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
			
			
		}
		cow=1;
		
		printf("Initial configuration have been copied from the file and implemented into the system successfully!!!!!!! \n \n  ");
		
			//printf("Go and yourself\n");
		
		fclose(vu);
		probation=0;
	}
	
	//printf("raam\n");
	
	
	
	
	
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
     
     
     	
     	FILE* lwu;
     	lwu=fopen("Diffusion_coefficient_el_l.txt","w");
     	fclose (lwu);
     	//FILE* dwu;
     	//dwu=fopen("Data_cord_dire.txt","w");
     	FILE* wu;
	
	wu=fopen("Widoms_ratio.txt", "w");
	
	fclose(wu);
	
	widoms=0.0;
	twidoms=0.0;
     
     
     
     
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
     
     	
     	
     	if (vfrap==0)
     	{
     		for (i=0; i<blockd; i++)
		{
			
			for (j=0; j<blockd; j++)
			{
				for (k=0; k<blockd; k++)
				{
					wdr[i][j][k]=0;
					for (l=0; l<80; l++)
					{
						dri[i][j][k][l]=6000;
					}
				}	
			}
		}
	     	
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
	wnlies=nlies;
	
	eb=teb;
	rul=1;
	sla=0;
	ve=0;
	ptr=1;
	
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

	printf("Maximum possible physical time with given step length and tsim is given by (tphy unit)-   ");
	printf("%d\n",iltt);
	printf("\n");
	printf("\n");
	printf("Probation period to generate the data is given by-  ");
	printf("%d\n", probation);
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("                                       ***  from the ruin of Humanity  ***  :::\n");
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
	
	
	
	
	int* flag=malloc(N*sizeof(int));
	
	if ( cow==0)
	{
	
		for(i=0; i<N; i++)
		{
			flag[i]=1;
		}
	}
	else
	{
		for(i=0; i<N; i++)
		{
			flag[i]=0;
		}
		
	}
	flag[N]=0;

	double* displacer=malloc(3*sizeof(double));
	int wdisplacer;
	double rotator;
	
	///printf(" %lf   %lf   %d\n", dss[0][0], bl[0], block[0]);
	//Simulation starts here!
	probation=probation+sitts;
	
	for(sitt=sitts; sitt<iltt; sitt++)
	{	
		
		
		
		if (sitt==probation)
		{
			nlies=nlies+1;
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
			}
		
		}
		
		for(ittt=0; ittt<lsd; ittt++)
		{	
		
		
			jet=0;
			
			
			if (cow==0)
			{
				ittt=0;
				cow=1;
				count=0;
				for (i=0; i<N; i++)
				{
					
					if (flag[i]==1)		{cow=0; count=count+1;}
					
						
					
				}
				////printf("You got %d number of particles left behind the line.\n", count);
				
				
				if(cow==1)
				{
					
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
						
					}
					
					
					FILE* ku;
     					ku=fopen("Compressed_system.txt","w");
     					fprintf(ku,"%d\n", sitt);
					fprintf(ku,"%lf\n", f[0]);
					fprintf(ku,"%lf\n", f[1]);
					fprintf(ku,"%lf\n", f[2]);
					for (i=0; i<N; i++)
					{
						fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
						
						
						
						
						fprintf(ku, "%d\n", wnrd[i]);
						for (j=0; j<wnrd[i]; j++) 
						{
							fprintf(ku,"   %d  %d", nrdi[i][j][0], nrdi[i][j][1]);
						}
						fprintf(ku, "\n");
						
						
						
						fprintf(ku,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", li[i][0], li[i][1], li[i][2], li[i][3], li[i][4], li[i][5], li[i][6], li[i][7], li[i][8]);
					}
					
					
					
					fclose(ku);
					printf("Compression happend!!! Now the actual simulation will start!!! %lf %lf %lf\n", f[0], f[1], f[2]);
					
					FILE* fu;
				     	fu=fopen("Initial_Cord.vtk","w");
				     	
					fprintf(fu,"# vtk DataFile Version 3.0\n");
					fprintf(fu,"Random data to test tensors\n");
					fprintf(fu,"ASCII\n");
					fprintf(fu,"DATASET POLYDATA\n");
					fprintf(fu,"POINTS %d float\n", N);
					
					for (jt=0; jt<N; jt++)
					{
						for(i=0;i<3;i++)
						{
					  		fprintf(fu,"%lf   ",ami[jt][i]);
						}
						fprintf(fu,"\n");
					}

					fprintf(fu,"\n");
					fprintf(fu,"POINT_DATA %d\n", N);
					fprintf(fu,"\n");
					fprintf(fu,"\n");
					fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
					fprintf(fu,"\n");
					
					for (i=0; i<N; i++)
					{
						rst=r[pid[i]];
						rsmt=rm[pid[i]];
						del[0]=(rsmt*li[i][0]*li[i][0])+(rsmt*li[i][3]*li[i][3])+(rst*li[i][6]*li[i][6]);
						del[1]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
						del[2]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
						del[3]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
						del[4]=(rsmt*li[i][1]*li[i][1])+(rsmt*li[i][4]*li[i][4])+(rst*li[i][7]*li[i][7]);
						del[5]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
						del[6]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
						del[7]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
						del[8]=(rsmt*li[i][2]*li[i][2])+(rsmt*li[i][5]*li[i][5])+(rst*li[i][8]*li[i][8]);
						
						fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
						fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
						fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
						fprintf(fu,"\n");
					}
					
					
					
					fclose(fu);
					
					
				}
						 
			}
			
			
			//printf("raammmmmmmmmdsfdfdsfdsfmm\n");
			
			
			
			//printf("raam\n");
			
			for (i=0; i<war; i++)
		     	{
		     		trac[i]=0;
		     	}
			wi=0;
			for(mo=0; mo<nlies; mo++)
			{
				
				
				
				
				if(mo==wnlies)
            			{
            				wi=1;
            				k=N;
            				rpid1=pid[k];
					ra2 = JKISS() / 4294967296.0;
	       				mxi=f[0]*ra2;
	       				myi=f[1]*drand48();
	       				mzi=f[2]*drand48();
            				twidoms=twidoms+1;
            				boss=0.4;
            				ami[k][0]=mxi;
            				ami[k][1]=myi;
            				ami[k][2]=mzi;
            				theta=acos(1.0-(2.0*drand48()));
					//ra2=drand48();
					ra2 = JKISS() / 4294967296.0;
	       				phi=2.0*piee*ra2;
	       				li[k][6]=sin(theta)*cos(phi);
	       				li[k][7]=sin(theta)*sin(phi);
	       				li[k][8]=cos(theta);
	       				
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
					
					rpid=pid[k];
					
					el[k][0]=(rsm[rpid]*li[k][0]*li[k][0])+(rsm[rpid]*li[k][3]*li[k][3])+(rs[rpid]*li[k][6]*li[k][6]);
					el[k][1]=(rsm[rpid]*li[k][0]*li[k][1])+(rsm[rpid]*li[k][3]*li[k][4])+(rs[rpid]*li[k][6]*li[k][7]);
					el[k][2]=(rsm[rpid]*li[k][0]*li[k][2])+(rsm[rpid]*li[k][3]*li[k][5])+(rs[rpid]*li[k][6]*li[k][8]);
					el[k][3]=(rsm[rpid]*li[k][0]*li[k][1])+(rsm[rpid]*li[k][3]*li[k][4])+(rs[rpid]*li[k][6]*li[k][7]);
					el[k][4]=(rsm[rpid]*li[k][1]*li[k][1])+(rsm[rpid]*li[k][4]*li[k][4])+(rs[rpid]*li[k][7]*li[k][7]);
					el[k][5]=(rsm[rpid]*li[k][1]*li[k][2])+(rsm[rpid]*li[k][4]*li[k][5])+(rs[rpid]*li[k][7]*li[k][8]);
					el[k][6]=(rsm[rpid]*li[k][0]*li[k][2])+(rsm[rpid]*li[k][3]*li[k][5])+(rs[rpid]*li[k][6]*li[k][8]);
					el[k][7]=(rsm[rpid]*li[k][1]*li[k][2])+(rsm[rpid]*li[k][4]*li[k][5])+(rs[rpid]*li[k][7]*li[k][8]);
					el[k][8]=(rsm[rpid]*li[k][2]*li[k][2])+(rsm[rpid]*li[k][5]*li[k][5])+(rs[rpid]*li[k][8]*li[k][8]);	
					
					
            				
            			}
            			else
            			{
            			
		    			//k=rand() % N;
		    			boss=drand48();
		    			
		    			//boss=0.4;
		    			
		    			ntt=N;
		    			/*for (i=0; i<war; i++)
		    			{
		    				if(trac[i]>tl[i])	ntt=ntt-nt[i];
		    			}*/
		    			k=rand()  % ntt;
		    			
		    			/*loo=0;
		    			for (i=0; i<war; i++)
		    			{
		    				if(trac[i]>tl[i])
		    				{
		    					if (k>loo) k=k+nt[i];
		    				}
		    				loo=loo+nt[i];
		    			}*/
		    			
		    			rpid1=pid[k];
		    			trac[rpid1]=trac[rpid1]+1;
		    			
		    			pf=0;
		    			loo=0;
		    				
		    			mxi=ami[k][0];
		    			myi=ami[k][1];
		    			mzi=ami[k][2];
		    		}
				
	    			l=(int)floor(mxi/bl[0]);
	    			m=(int)floor(myi/bl[1]);
	    			n=(int)floor(mzi/bl[2]);
				//printf("raamtlllllllllllllllll%d\n",k);
			
				
		    			
				displacer[0]=0.0;
				displacer[1]=0.0;
				displacer[2]=0.0;
				wdisplacer=0;
				rotator=0.0;
				
            			if (boss >= 0.5)
				{
					
					rdphi=tpiee*drand48();
					rdth=rdtheta[rpid1];
					
					xf=sin(rdth)*cos(rdphi);
					yf=sin(rdth)*sin(rdphi);
					zf=cos(rdth);
					dx=(-yf)*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
					dy=xf*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
					dz=0.0;
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
					
					
					theta=acos(li[k][8]);
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
					
					dli[0]=rdli[0]*q11+rdli[1]*q12+rdli[2]*q13;
					dli[1]=rdli[0]*q21+rdli[1]*q22+rdli[2]*q23;
					dli[2]=rdli[0]*q31+rdli[1]*q32+rdli[2]*q33;
					dli[3]=rdli[3]*q11+rdli[4]*q12+rdli[5]*q13;
					dli[4]=rdli[3]*q21+rdli[4]*q22+rdli[5]*q23;
					dli[5]=rdli[3]*q31+rdli[4]*q32+rdli[5]*q33;
					dli[6]=rdli[6]*q11+rdli[7]*q12+rdli[8]*q13;
					dli[7]=rdli[6]*q21+rdli[7]*q22+rdli[8]*q23;
					dli[8]=rdli[6]*q31+rdli[7]*q32+rdli[8]*q33;
					
					
					
					rst=rs[rpid1];
					rsmt=rsm[rpid1];
					del[0]=(rsmt*dli[0]*dli[0])+(rsmt*dli[3]*dli[3])+(rst*dli[6]*dli[6]);
					del[1]=(rsmt*dli[0]*dli[1])+(rsmt*dli[3]*dli[4])+(rst*dli[6]*dli[7]);
					del[2]=(rsmt*dli[0]*dli[2])+(rsmt*dli[3]*dli[5])+(rst*dli[6]*dli[8]);
					del[3]=(rsmt*dli[0]*dli[1])+(rsmt*dli[3]*dli[4])+(rst*dli[6]*dli[7]);
					del[4]=(rsmt*dli[1]*dli[1])+(rsmt*dli[4]*dli[4])+(rst*dli[7]*dli[7]);
					del[5]=(rsmt*dli[1]*dli[2])+(rsmt*dli[4]*dli[5])+(rst*dli[7]*dli[8]);
					del[6]=(rsmt*dli[0]*dli[2])+(rsmt*dli[3]*dli[5])+(rst*dli[6]*dli[8]);
					del[7]=(rsmt*dli[1]*dli[2])+(rsmt*dli[4]*dli[5])+(rst*dli[7]*dli[8]);
					del[8]=(rsmt*dli[2]*dli[2])+(rsmt*dli[5]*dli[5])+(rst*dli[8]*dli[8]);
					
							
						
					check=0;
					
					
					
					w1=dli[0];
					w2=dli[1];
					w3=dli[2];
					w4=dli[3];
					w5=dli[4];
					w6=dli[5];
					w7=dli[6];
					w8=dli[7];
					w9=dli[8];
					
 					for (ip=0; ip<wnrd[k]; ip++)
					{
						
						puc=nrdi[k][ip][0];
						jtt=nrdi[k][ip][1];
						
						rpid2=pid[puc];
						
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
						
						re1=re1+(float)tnkx[jtt][0]*f[0];
						re2=re2+(float)tnkx[jtt][1]*f[1];
						re3=re3+(float)tnkx[jtt][2]*f[2];
						
						
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
						
						if (tri==1) goto label15;
						
						ra= e[rpid1][1] * abrtm[2][0] + e[rpid1][2] * abrtm[1][0];
						rb= e[rpid2][1] * abrtm[0][2] + e[rpid2][2] * abrtm[0][1];
						if (fabs(ree[2] * rtm[1][0] - ree[1] * rtm[2][0]) > ra + rb) goto label15;
						
						ra= e[rpid1][1] * abrtm[2][1] + e[rpid1][2] * abrtm[1][1];
						rb= e[rpid2][0] * abrtm[0][2] + e[rpid2][2] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][1] - ree[1] * rtm[2][1]) > ra + rb) goto label15;
						
						ra= e[rpid1][1] * abrtm[2][2] + e[rpid1][2] * abrtm[1][2];
						rb= e[rpid2][0] * abrtm[0][1] + e[rpid2][1] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][2] - ree[1] * rtm[2][2]) > ra + rb) goto label15;
						
						ra = e[rpid1][0] * abrtm[2][0] + e[rpid1][2] * abrtm[0][0];
						rb = e[rpid2][1] * abrtm[1][2] + e[rpid2][2] * abrtm[1][1];
						if (fabs(ree[0] * rtm[2][0] - ree[2] * rtm[0][0]) > ra + rb) goto label15;
					
						ra= e[rpid1][0] * abrtm[2][1] + e[rpid1][2] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[1][2] + e[rpid2][2] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][1] - ree[2] * rtm[0][1]) > ra + rb) goto label15;
						
						ra= e[rpid1][0] * abrtm[2][2] + e[rpid1][2] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[1][1] + e[rpid2][1] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][2] - ree[2] * rtm[0][2]) > ra + rb) goto label15;
						
						ra= e[rpid1][0] * abrtm[1][0] + e[rpid1][1] * abrtm[0][0];
						rb= e[rpid2][1] * abrtm[2][2] + e[rpid2][2] * abrtm[2][1];
						if (fabs(ree[1] * rtm[0][0] - ree[0] * rtm[1][0]) > ra + rb) goto label15;
							
						ra= e[rpid1][0] * abrtm[1][1] + e[rpid1][1] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[2][2] + e[rpid2][2] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][1] - ree[0] * rtm[1][1]) > ra + rb) goto label15;
						
						ra= e[rpid1][0] * abrtm[1][2] + e[rpid1][1] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[2][1] + e[rpid2][1] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][2] - ree[0] * rtm[1][2]) > ra + rb) goto label15;
						
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
						
							
							if (pps[rpid1][rpid2]==1) { ec=1.0; ep=mr/ds[rpid1][rpid2]; goto label13; }	
							
							
							
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
							
							label13:
							if (ec<=1.0)
							{
								
								if (ep<1.0)
								{
							
									if (flag[k]==0 || flag[puc]==0)
									{
										check=1;
										jet=jet+1;
									}
									else
									{
										rotator=rotator+fabs(w7*u7+w8*u8+w9*u9)-fabs(li[k][6]*u7+li[k][7]*u8+li[k][8]*u9);
										
										
									}
									
								}
							
								
				
							}
							
							
						}
						
						label15:
						
						if (check==1)
						{
							break;
						}
					
					}
	 					
					if (check==0) 
					{	
						if (flag[k]==1)
						{
							//if (fabs(rotator)<0.0000001) flag[k]=0;
							if (rotator<0.00000001)	check=1;
						}
					}	
							
					
	
				
			
					if (check==0) 
					{
						
						// Variables to store the cross product (omega vector)
                                               

                                                // Calculate the cross product (li x dli)
                                                dx = li[k][7] * dli[8] - li[k][8] * dli[7];
                                                dy = li[k][8] * dli[6] - li[k][6] * dli[8];
                                                dz = li[k][6] * dli[7] - li[k][7] * dli[6];

                                                // Calculate the magnitude of the cross product vector (omega magnitude)
                                                bdm = sqrt(dx * dx + dy * dy + dz * dz);
                                                
                                                dx = dx/bdm;
                                                dy = dy/bdm;
                                                dz = dz/bdm;
                                                
                                        
                                                bdm = li[k][6] * dli[6] + li[k][7] * dli[7] + li[k][8] * dli[8];
                                                
                                                if (bdm <1.0) bdm = acos(bdm);
                                                else bdm  = 0.0;
                                                
                                                dx = dx*bdm;
                                                dy = dy*bdm;
                                                dz = dz*bdm;

						
						
						
						
						for (jt=0; jt<9; jt++)
						{
							li[k][jt]=dli[jt];
							el[k][jt]=del[jt];
						}
						
						
						lom[k][0]=lom[k][0]+dx;
						lom[k][1]=lom[k][1]+dy;
						lom[k][2]=lom[k][2]+dz;
					}
						
				}
                           	
                		else
				{
					          			
					trim=2;
					
	        			theta=acos(1.0-(2.0*drand48()));
					//ra2=drand48();
					ra2 = JKISS() / 4294967296.0;
	       				phi=2.0*piee*ra2;
					
					
					
					if(cow==0)
					{
						bdm=0.001;
						if (count<500 && fvf>0.64) bdm=0.0001;
						if (count<470 && fvf>0.64) bdm=0.00001;
						bdm=0.005;
						dx=bdm*sin(theta)*cos(phi);
						dy=bdm*sin(theta)*sin(phi);
						dz=bdm*cos(theta);
						
					}
				
					else
					{
						//printf("zppppppppppp\n");
						pdx=sin(theta)*cos(phi);
						pdy=sin(theta)*sin(phi);
						pdz=cos(theta);
						
						tpart=(li[k][6]*pdx+li[k][7]*pdy+li[k][8]*pdz)*step*grz[rpid1];
						x1per=(li[k][0]*pdx+li[k][1]*pdy+li[k][2]*pdz)*step*grxy[rpid1];
						y1per=(li[k][3]*pdx+li[k][4]*pdy+li[k][5]*pdz)*step*grxy[rpid1];
						
						//tpart=(li[k][6]*pdx+li[k][7]*pdy+li[k][8]*pdz)*step;
						//x1per=(li[k][0]*pdx+li[k][1]*pdy+li[k][2]*pdz)*step;
						//y1per=(li[k][3]*pdx+li[k][4]*pdy+li[k][5]*pdz)*step;
						
						dx=li[k][6]*tpart+li[k][0]*x1per+li[k][3]*y1per;
						dy=li[k][7]*tpart+li[k][1]*x1per+li[k][4]*y1per;
						dz=li[k][8]*tpart+li[k][2]*x1per+li[k][5]*y1per; 
						
						
						
					}
					
					
					xi=dx+ami[k][0];	
					yi=dy+ami[k][1];
					zi=dz+ami[k][2];
					
	        			if (xi >= f[0]) {xi=xi-f[0];}
	        			else if (xi < 0.0) {xi=xi+f[0];}
					if (yi >= f[1]) {yi=yi-f[1];}
	        			else if (yi < 0.0) {yi=yi+f[1];}
					if (zi >= f[2]) {zi=zi-f[2];}
	        			else if (zi < 0.0) {zi=zi+f[2];}
	        			
	        			
	        			
	        			
	        			xk=(int)floor(xi/bl[0]);
	        			yk=(int)floor(yi/bl[1]);
	        			zk=(int)floor(zi/bl[2]);
					//printf("raaaam2 %d %d %d %lf %lf %lf\n",xk,yk,zk,dxi,dyi,bl);
						
	    				for (j=0; j<UN; j++)	{pit[j][0]=6000; pit[j][1]=0;}
					wpi=0;
					
					if (vfrap==1)
					{
					
						for(jt=0;jt<27; jt++)
						{
							
							nxk=xk+tnkx[jt][0];
							nyk=yk+tnkx[jt][1];
							nzk=zk+tnkx[jt][2];
							
							
							if (nxk >= block[0])
							{
				    				nxk=0;
				   	 			cpsh=1;
							}
							
							else if (nxk==-1)
							{
				    				nxk=nxk+block[0];
				    				cpsh=-1;
							}
							else 
							{	
								cpsh=0;
							}
							
							
							if (nyk >= block[1])
			       				{    
								nyk=0;
				   			 	cph=1;
							}
							else if (nyk==-1)
							{
				    				nyk=nyk+block[1];
				    				cph=-1;
							}
							else 
							{	
								cph=0;
							}
							
							
							
							if (nzk >= block[2])
			       				{    
								nzk=0;
				   			 	cp=1;
							}
							else if (nzk==-1)
							{
				    				nzk=nzk+block[2];
				    				cp=-1;
							}
							else 
							{	
								cp=0;
							}
							
							
	 						
							dum=wdr[nxk][nyk][nzk];
							for (jdt=0; jdt<dum; jdt++)
							{
								da=dri[nxk][nyk][nzk][jdt];  
								rpid2=pid[da];         
								xf=ami[da][0]+(float)cpsh*f[0];
								yf=ami[da][1]+(float)cph*f[1];
								zf=ami[da][2]+(float)cp*f[2];
								
								jtt=(cpsh+1)*9+(cph+1)*3+cp+1;
								
								if (da != k) 
								{
					    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));
					    				bdm=dss[rpid1][rpid2];
					    				if (rijs <= dss[rpid1][rpid2])
					    				{
						    				if (flag[k]==0 || flag[da]==0)
						    				{ 
										        if (rijs <= ds[rpid1][rpid2])
											{
												trim=1;
												break;
											}   
											else
											{
												pit[wpi][0]=da;
												pit[wpi][1]=jtt;
												wpi=wpi+1;
											}
										}
										else
										{
											pit[wpi][0]=da;
											pit[wpi][1]=jtt;
											wpi=wpi+1;
										
										}
									}
								}
								
							} 
							
							
							if (trim==1) break;
							
						
						}
					}	
					
					
					
					if (trim == 2)
					{
						
						del1[0]=el[k][0];
						del1[1]=el[k][1];
						del1[2]=el[k][2];
						del1[3]=el[k][3];
						del1[4]=el[k][4];
						del1[5]=el[k][5];
						del1[6]=el[k][6];
						del1[7]=el[k][7];
						del1[8]=el[k][8];
						
						w1=li[k][0];
						w2=li[k][1];
						w3=li[k][2];
						w4=li[k][3];
						w5=li[k][4];
						w6=li[k][5];
						w7=li[k][6];
						w8=li[k][7];
						w9=li[k][8];
	 					for (ip=0; ip<wpi; ip++)
						{
							
							puc=pit[ip][0];
							jtt=pit[ip][1];
							rpid2=pid[puc];	
										
							u1=li[puc][0];
							u2=li[puc][1];
							u3=li[puc][2];
							u4=li[puc][3];
							u5=li[puc][4];
							u6=li[puc][5];
							u7=li[puc][6];
							u8=li[puc][7];
							u9=li[puc][8];
							
							re1=ami[puc][0]-xi;
							re2=ami[puc][1]-yi;
							re3=ami[puc][2]-zi;
							
							re1=re1+(float)tnkx[jtt][0]*f[0];
							re2=re2+(float)tnkx[jtt][1]*f[1];
							re3=re3+(float)tnkx[jtt][2]*f[2];
							
							
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
							
							for (i = 0; i < 3; i++)
								
							{	for (j = 0; j < 3; j++)
								{
								
									abrtm[i][j] = fabs(rtm[i][j]) + EPSILON;
								}
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
							
							
							
							//printf("tri %d\n", tri);
							if (tri==0)
							{
								
								//printf("tri %d\n", tri);
								
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
								lc=0;
								
								dive=0.0;
								bvs=0;
								bvsi=0;
								divea=0.0;
								divei=0.0;
								
								
								if (pps[rpid1][rpid2]==1) { ec=1.0; ep=mr/ds[rpid1][rpid2]; goto label12; }	
								
								while (trig>0.0001)
								{
									
									
									tri=tri+1;
									lc=lc+1;
									
									lemm=1.0-c;
									lem=c;
									lemms=lemm*lemm;
									lems=lem*lem;
									wt=0;
									
									for (i=0; i<3; i++)
									{
										for (j=0; j<3; j++)
										{
											lt[i][j]=del1[wt]*lemm+del2[wt]*lem;
											wt=wt+1;
										}
									}
									
									wt=0;
									
									for (i=0; i<3; i++)
									{
										for (j=0; j<3; j++)
										{
											lts[i][j]=del1[wt]*lemms-del2[wt]*lems;
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
											for (mt=0; mt<4; mt++)    lt[jt][mt]=lt[jt][mt]-((lt[it][mt]/lt[it][it])*dumm);
										
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
				    					//if (ec>1.0) trig=0.000001;
				    						
									ep=ec/mr;
				    					cm=mrp-dive;
									ec=ep*cm*cm;
				    					
					    				if (ec>1.0) trig=0.000001;
								}
								
								
								ep=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
								
								label12:
								if (ec<=1.0)
								{
									
									if (ep<1.0)
									{
								
										if (flag[k]==0 || flag[puc]==0)
										{
											trim=1;
											jet=jet+1;
										}
										else
										{
											displacer[0]=displacer[0]+re1+xi;
											displacer[1]=displacer[1]+re2+yi;
											displacer[2]=displacer[2]+re3+zi;
											wdisplacer=wdisplacer+1;
											
											//printf("Its firing on the place\n");
										}
									}
								
									
					
								}
							
								
							}
							
							
							label16:
							
							
							if (trim==1)
							{
								break;
							}
							
							
						}
						
					}
					
					
					if(wi==0)
					{
					
						if (trim==2)
						{
							if(flag[k]==1)
							{
								if (fabs(displacer[0]+displacer[1]+displacer[2])<0.0000001) { flag[k]=0;}
								else
								{
									displacer[0]=xi-(displacer[0]/wdisplacer);
									displacer[1]=yi-(displacer[1]/wdisplacer);
									displacer[2]=zi-(displacer[2]/wdisplacer);
									bdm=dx*displacer[0]+dy*displacer[1]+dz*displacer[2];
									if(bdm>=0.0) trim=2;
									else trim=1;
								}
							}
							
						}
						
						if (trim==2)
						{
							dum=wdr[l][m][n];
							for (j=0;j<dum;j++)
							{	
								
								if (dri[l][m][n][j]==k)
								{
									for (tric=j; tric<dum; tric++)	dri[l][m][n][tric]=dri[l][m][n][tric+1];
									wdr[l][m][n]=wdr[l][m][n]-1;	
									
									break;
								}
										
							}
							if (vfrap==1)
							{
								dum=wdr[xk][yk][zk];
								dri[xk][yk][zk][dum]=k;
								wdr[xk][yk][zk]=dum+1;
							}
							ami[k][0]=xi;
			    				ami[k][1]=yi;
			    				ami[k][2]=zi;
			    				
							ldcxy[k][0]=ldcxy[k][0]+x1per;
							ldcxy[k][1]=ldcxy[k][1]+y1per;
							ldcz[k]=ldcz[k]+tpart;
							ldcxyz[k][0]=ldcxyz[k][0]+dx;
							ldcxyz[k][1]=ldcxyz[k][1]+dy;
							ldcxyz[k][2]=ldcxyz[k][2]+dz;
							
								
							mum=wnrd[k];
							
							for (jt=0; jt<mum; ++jt)
		    					{
								chum=nrdi[k][jt][0];
								drum=wnrd[chum];
								for (sj=0; sj<drum; ++sj)
								{
									if (nrdi[chum][sj][0]==k)
									{
										for (tric=sj; tric<drum; tric++)
										{
											nrdi[chum][tric][0]=nrdi[chum][(tric+1)][0];
											nrdi[chum][tric][1]=nrdi[chum][(tric+1)][1];
										}
										wnrd[chum]=drum-1;
										break;
									}
								
								}
								
								nrdi[k][jt][0]=6000;
								nrdi[k][jt][1]=13;
								
								
							}
								
							wnrd[k]=0;	
							for (jdt=0; jdt < wpi; jdt++)
							{	
								//printf("raam is ");
								
								
								td=pit[jdt][1];
								nrdi[k][jdt][1]=td;
								wnrd[k]=wnrd[k]+1;
									
								td=pit[jdt][0];
								nrdi[k][jdt][0]=td;
								tum=wnrd[td];
								
								nrdi[td][tum][0]=k;
								nrdi[td][tum][1]=26-pit[jdt][1];
								wnrd[td]=tum+1;
							
							}	
							
							
							
						}
					}
					else
					{
						if (trim==2) widoms=widoms+1;
					}
	
				}
				
				
			//printf("snn750 %d\n",mo);
			}
			
			
			
					
				
		}
		printf("Simulation is going on si %d \n", sitt);
		
		FILE* wu;
	
		wu=fopen("Widoms_ratio.txt", "a");
		fprintf(wu, "%d   %lf  %lf  %lf   %lf\n", sitt+1, widoms, twidoms,  widoms/twidoms, -log(widoms/ twidoms));
	
		fclose(wu);
		
		
		
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitt)/100.0));

		FILE* cfu;
		cfu=fopen(filename,"a");
		
			
		fprintf(cfu,"%lf\n", f[0]);
		fprintf(cfu,"%lf\n", f[1]);
		fprintf(cfu,"%lf\n", f[2]);
		for (i=0; i<N; i++)
		{
			fprintf(cfu,"%lf   %lf   %lf\n", ami[i][0], ami[i][1], ami[i][2]);	
			fprintf(cfu, "  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", ldcxyz[i][0], ldcxyz[i][1], ldcxyz[i][2], ldcxy[i][0], ldcxy[i][1], ldcz[i], lom[i][0], lom[i][1], lom[i][2]);
			
			fprintf(cfu,"%lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf\n", li[i][0], li[i][1], li[i][2], li[i][3], li[i][4], li[i][5], li[i][6], li[i][7], li[i][8]);
			//fprintf(cfu,"%lf  %lf  %lf\n", fau[i][0], fau[i][1], fau[i][2]);
		}
		
		
		fclose(cfu);
		
		
		
		
			
		
		
		//fprintf(dwu,"\n");
		
		if (((float)(sitt+1)/10.0)==((float)floor((float)(sitt+1)/10.0)))
		{
		
		
			
			jt=0;
		
			for (j=0; j<war; j++)
			{
				dif=0.0;
				difu=0.0;
				difus=0.0;
				difusl=0.0;
				jtt=jt+nt[j];
				for (i=jt; i<jtt; i++)
				{	
					dif=dif+(ldcxy[i][0]*ldcxy[i][0]+ldcxy[i][1]*ldcxy[i][1]);
					difu=ldcz[i]*ldcz[i]+difu;
					difus=difus+(ldcxyz[i][0]*ldcxyz[i][0]+ldcxyz[i][1]*ldcxyz[i][1]+ldcxyz[i][2]*ldcxyz[i][2]);
					difusl=difusl+((lom[i][0]*lom[i][0])+(lom[i][1]*lom[i][1])+(lom[i][2]*lom[i][2]));
					
					
				}
				
				dif=dif/(float)nt[j];
				
				difu=difu/(float)nt[j];

				difus=difus/(float)nt[j];
				
				difusl=difusl/(float)nt[j];
				
				dif=dif/((float)sitt-probation);
				
				difu=difu/((float)sitt-probation);

				difus=difus/((float)sitt-probation);
				
				difusl=difusl/((float)sitt-probation);
				
				
				
				printf("medium %d  %lf  %lf  %lf %lf \n", (sitt+1), dif, difu, difus, difusl);
				
				FILE* llwwu;
     				llwwu=fopen("Diffusion_coefficient_el_l.txt","a");
				
				fprintf(llwwu,"%d  %lf  %lf  %lf  %lf  \n", (sitt+1), dif, difu, difus, difusl);
				
				fclose(llwwu);
				
				
				jt=jtt;
				
			}
			
			jt=0;
						
			
			
			
			
			
			
			
			
			
		}
		
		
		
		
		if (((float)(sitt+1)/100.0)==((float)floor((float)(sitt+1)/100.0)))
		{
		
		
			
			
			char filename[25] = {0};

    			sprintf(filename, "%d.txt", (int)floor((float)(sitt+1)/100.0));

			FILE* cfu;
			cfu=fopen(filename,"w");
			
			
			fclose(cfu);
			
			
			
			FILE* ku;
			ku=fopen("Relaxed_system.txt","w");
			
			fprintf(ku,"%d\n", sitt+1);
			
			fprintf(ku,"%lf\n", f[0]);
			fprintf(ku,"%lf\n", f[1]);
			fprintf(ku,"%lf\n", f[2]);
			
			for (i=0; i<N; i++)
			{
				fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
				
				 
				 
				fprintf(ku, "%d\n", wnrd[i]);
				for (j=0; j<wnrd[i]; j++) 
				{
					fprintf(ku,"   %d  %d", nrdi[i][j][0], nrdi[i][j][1]);
				}
				
				fprintf(ku, "\n");
				
				
				
				fprintf(ku,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", li[i][0], li[i][1], li[i][2], li[i][3], li[i][4], li[i][5], li[i][6], li[i][7], li[i][8]);
				
				//fprintf(ku,"%lf  %lf  %lf\n", fau[i][0], fau[i][1], fau[i][2]);
			}
			
			fclose(ku);
			
			
			
			
			//VTK Zone starts!
	
			FILE* fu;
		     	fu=fopen("Cord.vtk","w");
		     	
			fprintf(fu,"# vtk DataFile Version 3.0\n");
			fprintf(fu,"Random data to test tensors\n");
			fprintf(fu,"ASCII\n");
			fprintf(fu,"DATASET POLYDATA\n");
			fprintf(fu,"POINTS %d float\n", N);
			
			for (jt=0; jt<N; jt++)
			{
				for(i=0;i<3;i++)
				{
			  		fprintf(fu,"%lf   ",ami[jt][i]);
				}
				fprintf(fu,"\n");
			}

			fprintf(fu,"\n");
			fprintf(fu,"POINT_DATA %d\n", N);
			fprintf(fu,"\n");
			fprintf(fu,"\n");
			fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
			fprintf(fu,"\n");
			
			for (i=0; i<N; i++)
			{
				rst=r[pid[i]];
				rsmt=rm[pid[i]];
				del[0]=(rsmt*li[i][0]*li[i][0])+(rsmt*li[i][3]*li[i][3])+(rst*li[i][6]*li[i][6]);
				del[1]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
				del[2]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
				del[3]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
				del[4]=(rsmt*li[i][1]*li[i][1])+(rsmt*li[i][4]*li[i][4])+(rst*li[i][7]*li[i][7]);
				del[5]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
				del[6]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
				del[7]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
				del[8]=(rsmt*li[i][2]*li[i][2])+(rsmt*li[i][5]*li[i][5])+(rst*li[i][8]*li[i][8]);
				
				fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
				fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
				fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
				fprintf(fu,"\n");
			}
			
			
			
			fclose(fu);

			//VTK Zone ends!
			
			
			
			
		}
		
		
		
		
		
		

	}
	
	
	

	//SImulaton ends!





	//FINAL CORD Zone starts!

	
	
	//FINAL CORD Zone ends!







	




	
	//The war of Freedom started!
	

	for (i=0;i<N;++i)
	{
		free(li[i]);
	}
		
	free(li);
	
	
	for (i=0;i<N;++i)
	{
		free(el[i]);
	}

	free(el);
	
	for (i=0;i<N;++i)
	{
		free(ami[i]);
	}
		
	
	free(ami);
	

	
	
	
	
	for (k=0; k<N; k++)
	{
		for (j=0; j<UN; j++) free(nrdi[k][j]);
		free(nrdi[k]);
	
	}
		
	free(nrdi);
	free(wnrd);
	
	
		
	for (i=0; i<blockd; i++)
	{
		
		for (j=0; j<blockd; j++)
		{
		
			for (k=0; k<blockd; k++)
			{
				free(dri[i][j][k]);
			}
			free(dri[i][j]);
			free(wdr[i][j]);
		}
		free(dri[i]);
		free(wdr[i]);
		
	}
	free(dri);
	free(wdr);
	

	//Now freedom is achieved!	


	
}




