#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <omp.h>
const int E=20, UN=200, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel, tbox, vtbox;

//void ecfcal(int,double[]);

int  N, trim, check, checkl, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, pnlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, it, nf,  block[3], blockd, jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct, rul, sla, veei, veec, SHA, ptr, ptrr, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt, gv, bvs, bvsi, mpatch, tpatch, count, bpc1, bpc2, vc, pjet, tpjet, sjet, tsjet, tsitt, stsitt, ptt, trims, initial_fp, latrf, bshape_p, boxfx, boxfy, cp, cph, cpsh;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9], f[3], bl1, bl2, bl3, bl[3], f1, f2, f3, ar, per, fr[3], ax1x, ax1y, ax1z, fva[3][3], fvi[3][3], fv[3][3], fvao[3][3], fvio[3][3], fvo[3][3], reserve, bleol, bleop, roott; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda,lembdai, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, fe, fer, tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON,  tbl, pbc,  dif, difu, difus, difl, difusl, adefus, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm, dsss, pree, pre, vol, dlr, press, fer, ene, bf, volo, voln, fr[3], vstep, vpstep, divea, divei, vpress, flim, slredo, slred, dfv, rvstep;  




int main(void)
{
	
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
	
	roott=sqrt(2.0);
	
	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;
	
	
	
	
	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	
	fscanf(bu,"%lf",&bdm);
	initial_fp=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	bshape_p=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	N=(int)round(bdm);
	//if (N==0) printf("raam tere desh men \n");
	// Definition of universal constants starts!
	
	
	
	// Definition of universal constants ends!
	
	
	
	// Definition of global constants starts!
	
	//step=0.01;
	//fscanf(bu,"%lf",&step);
	//fscanf(bu,"%ld",&tsim);
	fscanf(bu,"%lf",&bdm);
	step=bdm;
	fscanf(bu,"%lf",&bdm);
	tsim=(int)round(bdm);
	
	//tsim=100000;
	ivf=0.004;
	
	// Definition of global constants ends       !!!!!!!!!!!!!!!!
	
					
	// Correlation function        !!!!!!!!!!!!!!!!
	
	
	
	
	// Correlation function declaration end         !!!!!!!!!!!!!!!!			
	
	
	
	
  	// Definition of simulation constants starts        !!!!!!!!!!!!!!!!
	

	fscanf(bu,"%lf",&bdm);
	press=bdm;
	
	
	
	
	
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
	
	
	
	
	double epsi;
	
	// Declaration of core bond variable ends         !!!!!!!!!!
	
	
	
	// Evaluation of bond variable starts           !!!!!!!!!
	
	
	
	bn=40;
	bnll=40;
	sst=1;
	
	
	//double epsi1, epsi2, epsi3, omega2, beta1, beta2, beta3, betai;
	
	// Declaration of core bond variable ends!
	
	
	// Evaluation of bond variable starts!
	

	
	fscanf(bu, "%lf", &bdm);
	
	tpatch=(int)bdm;
	
	

	
	
	
	int** pop=malloc(tpatch*sizeof(int*));
	for (i=0; i<tpatch; i++)
	{
		pop[i]=malloc(tpatch*sizeof(int));
	}
	
	
	
	
	jt=0;
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch; j++)
		{
			pop[i][j]=0;
			pop[j][i]=0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				pop[i][j]=(int)round(bdm);
				pop[j][i]=(int)round(bdm);
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	// Patchy potential properties;
	
	double** epsi2=malloc(tpatch*sizeof(double*));
	for(i=0; i<tpatch; i++)
	{
		epsi2[i]=malloc(tpatch*sizeof(double));
	}
	
	double** beta2=malloc(tpatch*sizeof(double*));
	for(i=0; i<tpatch; i++)
	{
		beta2[i]=malloc(tpatch*sizeof(double**));
	
	}
	
	double** ur=malloc(tpatch*sizeof(double*));
	for(i=0; i<tpatch; i++)
	{
		ur[i]=malloc(tpatch*sizeof(double)); 
	
	}
	
	
	
	double* omega2=malloc(tpatch*sizeof(double));
	
	double* om2=malloc(tpatch*sizeof(double));
	
	double* wom=malloc(tpatch*sizeof(double));
	
	double* mom=malloc(tpatch*sizeof(double));
	
	
	
	
	for(i=0; i<tpatch;  i++)
	{
		for(j=i; j<tpatch;  j++)
		{
			epsi2[i][j]=0.0;
			epsi2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch;  j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi2[i][j]=bdm;
				epsi2[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch; j++)
		{
			beta2[i][j]=0.0;
			beta2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta2[i][j]=bdm;
				beta2[j][i]=bdm;
				
				ur[i][j]=log(1.0-(1.0/(1.0+beta2[i][j]))+EPSILON);
				ur[j][i]=ur[i][j];
				
				
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	for (i=0; i<tpatch; i++)
	{
		omega2[i]=0.0;
	}
	
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<tpatch; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			omega2[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	for (i=0; i<tpatch; i++)
	{
		om2[i]=0.0;
	}
	
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<tpatch; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			om2[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	for (i=0; i<tpatch; i++)
	{
		wom[i]=omega2[i]+om2[i];
		mom[i]=omega2[i]-om2[i];
	}
	
	
	
	
	
	int* npatch=malloc(war*sizeof(int));
	
	for (i=0; i<war; i++)
	{
		npatch[i]=0;
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
			npatch[i]=(int)bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	int** patches=malloc(war*sizeof(int*));
	
	for (i=0; i<war; i++)
	{
		patches[i]=malloc(npatch[i]*sizeof(int));
	}
	
	for(i=0; i<war; i++)
	{
		for (j=0; j<npatch[i]; j++)
		{
			patches[i][j]=0;
		}
	}
	
	fscanf(bu, "%lf", &bdm);
	for (i=0; i<war; i++)
	{
		
		for (j=0; j<npatch[i]; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				patches[i][j]=(int)bdm;
			}
		}
		if (bdm==50.0) 
		{
			break;
		}
		
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	

	fclose(bu);
	remove("Base_data.txt");
	
	
	
	
	printf("\n\n\n\nBuddha was born to decore the world \n\n\n\n");

	
	
	
	
	FILE* bbu;
	bbu=fopen("Patch_vectors.txt", "r");
	
	
	double*** orp=malloc(war*sizeof(double**));
	 
	for(i=0; i<war; i++)
	{
		orp[i]=malloc(npatch[i]*sizeof(double*));
		for (j=0; j<npatch[i]; j++)
		{
			orp[i][j]=malloc(3*sizeof(double));
		}
	
	}
	
	
	for (i=0; i<war; i++)
	{
		//printf("Patch orientation for the particle type %d is given by-   \n", i );
		for (j=0; j<npatch[i]; j++)
		{
			//printf("Patch vector %d:  ", j );
			for(k=0; k<3; k++)
			{
				//printf("John!!!!!!! %d\n", k);
				fscanf(bbu, "%lf  ", &orp[i][j][k]);
			
			}
			
			//printf("\n");
		}
	}
	
	fclose(bbu);
	
	
	
	
	//printf("Buddha was born to decore the world \n");
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		for (j=0; j<npatch[i]; j++)
		{
			bdm=sqrt(orp[i][j][0]*orp[i][j][0]+orp[i][j][1]*orp[i][j][1]+orp[i][j][2]*orp[i][j][2]);
			orp[i][j][0]=orp[i][j][0]/bdm;
			orp[i][j][1]=orp[i][j][1]/bdm;
			orp[i][j][2]=orp[i][j][2]/bdm;
		}
	}
	
	
	
	
	
	
	
	
	mpatch=0;
	for (i=0; i<war; i++)
	{
		if (npatch[i]>mpatch) mpatch=npatch[i];
	}
	
	
	double**  tpaxe=malloc(mpatch*sizeof(double*));
	for (i=0; i<mpatch; i++)
	{
		tpaxe[i]=malloc(3*sizeof(double));
	}
	// Declaired!
	
	//printf("Hey here is the mpatch %d\n", mpatch);
	
	double** beta=malloc(tpatch*sizeof(double*));
	for(i=0; i<tpatch; i++)
	{
		beta[i]=malloc(tpatch*sizeof(double)); 
	
	}
	
	
	epsi=0.0;
	
	
	for(i=0; i<tpatch; i++)
	{
		for(j=i; j<tpatch; j++)
		{
			if (pop[i][j]==1)
			{
				if (epsi2[i][j]>epsi) epsi=epsi2[i][j];
			}
		}		
	}	

	
	fvf=0.004;
	ivf=fvf;
	
	
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
	printf("Assigned pressure-    %lf \n", press);
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
	
	
	printf("Number of All type of patches defined is given as (Not necessarily decorated on the particle): %d\n\n\n\n\n", tpatch);
	
	
	printf("Potential parameters associated with the patches are given by:- \n");
	printf("\n");
			 printf("\n");
			 printf("\n");
	
	for (i=0; i<tpatch; i++)
	{
		for (j=0; j<tpatch; j++)
		{
			 printf("Kind of potential between patch type %d and %d is defined as-  \n ", i,j);
			 if (pop[i][j]==0)
			 {
			 	printf("Non interacting patches.\n");
			 	
			 	
			 }
			 else
			 {
			 	printf("Epsi-   %lf\n", epsi2[i][j]);
			 	printf("Beta-   %lf, Temperature-   kBT/u_0 = %lf\n", beta2[i][j], -1.0/ur[i][j]);
			 	printf("Omega lower limit- %lf   %lf,  omega upper limit- %lf   %lf \n", mom[i], mom[j], wom[i], wom[j]);
			 }
			 
			 printf("\n");
			 printf("\n");
		}
	
	}
	
	printf("Selected Epsi to harnesh the particle in the range %lf \n", epsi);
	//printf("Values corresponding to the particular parameter is given by\n");
	
	
	for (i=0; i<war; i++)
	{
		if (npatch[i]>0)printf("Number of patches on type %d particle is: %d\n", i, npatch[i]);
		else printf("No patch defined on type %d particles\n", i);
	}
	
	for (i=0; i<war; i++)
	{
		if (npatch[i]>0)
		{	printf("Patch decoration in the particle type %d is given by-   ", i);
			for (j=0; j<npatch[i]; j++)
			{
				printf("%d  ", patches[i][j]);
			}
		}
		printf("\n\n\n");
	}
	
	for (i=0; i<war; i++)
	{
		if(npatch[i]>0)
		{
			printf("Patch orientation for the particle type %d is given by-   \n", i );
			for (j=0; j<npatch[i]; j++)
			{
				printf("Patch vector %d:  ", j );
				for(k=0; k<3; k++)
				{
					///printf("John!!!!!!! %d\n", k);
					printf("%lf  ", orp[i][j][k]);
				
				}
				
				printf("\n");
			}
		}
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	/*
	printf("Applicable parameters is given as:   \n");
	
	printf("Epsi %lf \n", epsi);
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			if (pop[i][j]>0)
			{
				printf("Applied beta is given between particles %d and %d is given by- %lf   \n", i, j, beta[i][j]);
			
				if (pop[i][j]>2) printf("Applied betai is given between particles %d and %d is given by- %lf   \n", i, j, betai[i][j]);
			}
			else
			{
				printf("Not applicable!!!\n");
			}
		}
		
	}
	// Printing things starts now !!!!!!!
	
	*/
	
	
	
	
	
	
	
	
	
	
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
	
	if (bshape_p==1) ble=ble*1.2;
	
	bleop=ble;
	
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
	
	
	int* pid=malloc(N*sizeof(int));
	
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
	
	
	
	
	
	// Declaration of core variables starts!
	
	int** pit=malloc(UN*sizeof(int*));
	for (i=0; i<UN; i++)
	{
		pit[i]=malloc(2*sizeof(int));
	}
	
	int wpi;
	
	double** li=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		li[i]=malloc(9*sizeof(double));
	}
	
	
	
	double*** paxe=malloc(N*sizeof(double**));
	
	for (i=0; i<N; i++)
	{
		paxe[i]=malloc(npatch[pid[i]]*sizeof(double*));
		for(j=0; j<npatch[pid[i]]; j++)
		{
			paxe[i][j]=malloc(3*sizeof(double));
		}
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
	
	
	double** lfp=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		lfp[i]=malloc(3*sizeof(double));
	}
	
	
	
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
	
	
	int*** tnrdi=malloc(N*sizeof(int**));
	int* wtnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		tnrdi[i]=malloc(UN*sizeof(int*));
		for (j=0; j<UN; j++)
		{
			tnrdi[i][j]=malloc(2*sizeof(int));
		}
		
	}
	
	
	
	
	
	
	
	tbox=blockd*blockd*blockd;
	
	int** dri=malloc(tbox*sizeof(int*));
	int* wdr=malloc(tbox*sizeof(int));
		
	for (vtbox=0; vtbox<tbox; vtbox++)
	{
		dri[vtbox]=malloc(80*sizeof(int));
	}
	
	// Declaration of compresser variables starts!
	
	int vee;
	
	// Declaration of compresser variables ends!
	
	// Declaration of core variables ends!
	
	
	
	
	
	
	// Declaration of bonding variable starts!
	
	
	
	
	int*** bobb=malloc(N*sizeof(int**));
	int* wbobb=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		bobb[k]=malloc(E*sizeof(int*));
		for (i=0; i<E; i++)
		{
			bobb[k][i]=malloc(3*sizeof(int));
		}
	
	}	



	
	
	int*** cbobb=malloc(N*sizeof(int**));
	int* wcbobb=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		cbobb[k]=malloc(E*sizeof(int*));
		for (i=0; i<E; i++)
		{
			cbobb[k][i]=malloc(3*sizeof(int));
		}
	}	
	
	
	int wpbobb;
	int psbobb[E][3];

	int bumbque[E][3];
	
	
	int*** hold=malloc(N*sizeof(int**));
	int* whold=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		hold[k]=malloc(E*sizeof(int*));
		for (i=0; i<E; i++)
		{
			hold[k][i]=malloc(3*sizeof(int));
		}
	}	

	

	

	int  tnkx[27][3];
	
		
		
	

	
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
		wtnrd[i]=0;
		for (j=0; j<UN; j++)
		{
			nrdi[i][j][0]=6000;
			nrdi[i][j][1]=13;
			tnrdi[i][j][0]=6000;
			tnrdi[i][j][1]=13;
		
		}			
		
	}
	
	for (vtbox=0; vtbox<tbox; vtbox++)
	{
		
		
		wdr[vtbox]=0;
		for (l=0; l<80; l++)
		{
			dri[vtbox][l]=6000;
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
	
	
	
		
	
	
	
	
	
	
	
	
	// Initialization of core variables ends!
	
	
	
	
	
	
	
	
	
	// Initialization of bond variables starts!
	
	
	for (i=0; i<N; i++)
	{
		wbobb[i]=0;
		wcbobb[i]=0;
		
		
		for (j=0; j<E; j++)
		{
			bobb[i][j][0]=6000;
			cbobb[i][j][0]=6000;
			
			bobb[i][j][1]=6000;
			cbobb[i][j][1]=6000;
			
			
			bobb[i][j][2]=6000;
			cbobb[i][j][2]=6000;
			
		}			
	}
		
	wpbobb=0;
	for(i=0; i<E; i++)
	{
		psbobb[i][0]=6000;
		psbobb[i][1]=6000;
		psbobb[i][2]=6000;

	}
		
	

	
	
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
	
	
	if (vu==NULL  || initial_fp==0)
	{
		printf("Initial configuration within the simulation box have been started generating: \n \n");
		
		boxfx=block[1]*block[2];   boxfy=block[2];
		
		fv[0][0]=1.0;  fv[0][1]=0.0;  fv[0][2]=0.0;
		fv[1][0]=0.0;  fv[1][1]=1.0;  fv[1][2]=0.0;
		fv[2][0]=0.0;  fv[2][1]=0.0;  fv[2][2]=1.0;
		
		fvi[0][0]=1.0;  fvi[0][1]=0.0;  fvi[0][2]=0.0;
		fvi[1][0]=0.0;  fvi[1][1]=1.0;  fvi[1][2]=0.0;
		fvi[2][0]=0.0;  fvi[2][1]=0.0;  fvi[2][2]=1.0;
		
		fva[0][0]=f[0];  fva[0][1]=0.0;  fva[0][2]=0.0;
		fva[1][0]=0.0;  fva[1][1]=f[1];  fva[1][2]=0.0;
		fva[2][0]=0.0;  fva[2][1]=0.0;  fva[2][2]=f[2];
		
		fvo[0][0]=fv[0][0]; fvo[0][1]=fv[0][1]; fvo[0][2]=fv[0][2];
		fvo[1][0]=fv[1][0]; fvo[1][1]=fv[1][1]; fvo[1][2]=fv[1][2];
		fvo[2][0]=fv[2][0]; fvo[2][1]=fv[2][1]; fvo[2][2]=fv[2][2];
		
		fvio[0][0]=fvi[0][0]; fvio[0][1]=fvi[0][1]; fvio[0][2]=fvi[0][2];
		fvio[1][0]=fvi[1][0]; fvio[1][1]=fvi[1][1]; fvio[1][2]=fvi[1][2];
		fvio[2][0]=fvi[2][0]; fvio[2][1]=fvi[2][1]; fvio[2][2]=fvi[2][2];
		
		fvao[0][0]=fva[0][0]; fvao[0][1]=fva[0][1]; fvao[0][2]=fva[0][2];
		fvao[1][0]=fva[1][0]; fvao[1][1]=fva[1][1]; fvao[1][2]=fva[1][2];
		fvao[2][0]=fva[2][0]; fvao[2][1]=fva[2][1]; fvao[2][2]=fva[2][2];
		
		volo=fabs(fva[0][0]*(fva[1][1]*fva[2][2]-fva[2][1]*fva[1][2])-fva[0][1]*(fva[1][0]*fva[2][2]-fva[2][0]*fva[1][2])+fva[0][2]*(fva[1][0]*fva[2][1]-fva[2][0]*fva[1][1]));
		
		ec=0.0;
		for (i=0; i<3; i++)
		{
			for(j=i+1; j<3; j++)
			{
				ep=acos(fv[i][0]*fv[j][0]+fv[i][1]*fv[j][1]+fv[i][2]*fv[j][2]);
				///printf("%lf  %d  %d \n", ep, i, j);
				
				if (ep>piee/2.0) ep=piee-ep;
				
				ep=piee/2.0-ep;
				if (ep<0.0000001) ep=0.0;
				bdm=bleop/cos(ep);
				if (bdm>ec) ec=bdm;
			}
		
		}
		
		ble=ec;
		bleol=ble;
		
		for(i=0; i<N; i++)	
	    	{	
			
			
			
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
			
			cow=0;
			while (cow==0)	
			{	
				
				xi=((float)(rand() % ix))*bl[0]+(bl[0]/2.0);
		    		yi=((float)(rand() % iy))*bl[1]+(bl[1]/2.0);
		    		zi=((float)(rand() % iz))*bl[2]+(bl[2]/2.0);
		    		xk=(int)floor(xi/bl[0]);
		    		yk=(int)floor(yi/bl[1]);
		    		zk=(int)floor(zi/bl[2]);
		    		//printf("i, iz %d %lf %lf %lf %d \n",i, xi, yi, zi, wdr[xk][yk][zk]);
		    		
		    		vtbox=boxfx*xk+boxfy*yk+zk;
		    		dum=wdr[vtbox];
		    	
				if (dum==0)
				{
					
					dri[vtbox][dum]=i;
					wdr[vtbox]=dum+1;
					ami[i][0]=xi;
					ami[i][1]=yi;
					ami[i][2]=zi;
					cow=1;
					//fprintf(gu,"%lf  %lf  %lf \n", fau[i][0], fau[i][1], fau[i][2]);
				}
			}
			
			lfp[i][0] = ami[i][0]*fv[0][0] + ami[i][1]*fv[1][0] + ami[i][2]*fv[2][0];
		    	lfp[i][1] = ami[i][0]*fv[0][1] + ami[i][1]*fv[1][1] + ami[i][2]*fv[2][1];
		    	lfp[i][2] = ami[i][0]*fv[0][2] + ami[i][1]*fv[1][2] + ami[i][2]*fv[2][2];
			
			//printf("number filled %d\n", i);
		}
		
		
		
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
			
			rpid=pid[k];
			for (i=0; i<npatch[rpid]; i++)
			{
				ax1x=orp[rpid][i][0];
				ax1y=orp[rpid][i][1];
				ax1z=orp[rpid][i][2];
				
				paxe[k][i][0]=ax1x*q11+ax1y*q12+ax1z*q13;
				paxe[k][i][1]=ax1x*q21+ax1y*q22+ax1z*q23;
				paxe[k][i][2]=ax1x*q31+ax1y*q32+ax1z*q33;
			}
			
			
			
			/*
			 ///Randomization of the patch distribution is done by the distribution of patch vector not the the type of patch of patch identity which is limited to finite (finite tyep of bonds)...while orientation can be infinite way... this is just an example of patch randomization
				rpid=pid[k];
			//// Randomization of the patch vector !!!!!!!!!!!!!!!
			it=rand()  % 3;
			for (i=0; i<npatch[rpid]; i++)
			{
				j=it+i;
				if (j>=npatch[rpid]) j=j-npatch[rpid];
				
				//printf("%d %d\n", k, j);
				
				ax1x=orp[rpid][j][0];
				ax1y=orp[rpid][j][1];
				ax1z=orp[rpid][j][2];
				
				paxe[k][i][0]=ax1x*q11+ax1y*q12+ax1z*q13;
				paxe[k][i][1]=ax1x*q21+ax1y*q22+ax1z*q23;
				paxe[k][i][2]=ax1x*q31+ax1y*q32+ax1z*q33;
			}*/
		
			////printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \n", q11, q21, q31, q12, q22, q32, q13, q23, q33);
		}
		
		
			
		
		
		FILE* rpvf;
		rpvf=fopen("Resulted_patch_vector.txt", "w");
		
		for (i=0; i<N; i++)
		{
			for (j=0; j<npatch[pid[i]]; j++)
			{
				fprintf(rpvf, "%lf   %lf   %lf \n", paxe[i][j][0],   paxe[i][j][1] , paxe[i][j][2] );
			
			}
		}
		fclose(rpvf);
		
		cow=0;
		sitts=0;
		quao=0.0;
		quaoi=0.0;
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		fclose(cfu);
		printf("Initial configuration generated successfully!!!!!!! \n \n");
		
		
		
		
		
		
	}
	else
	{
		printf("A pre-prepared system found and intial configuration will be set accordingly. \n");
		printf("Initial configuration is being read from the available file!!! \n \n  ");
		cow=1;
			
		fscanf(vu, "%d", &sitts);
			
		fscanf(vu,"%lf", &f[0]);
		fscanf(vu,"%lf", &f[1]);
		fscanf(vu,"%lf", &f[2]);
		
		fscanf(vu,"%lf", &fv[0][0]); fscanf(vu,"%lf", &fv[0][1]);  fscanf(vu,"%lf", &fv[0][2]);
		fscanf(vu,"%lf", &fv[1][0]); fscanf(vu,"%lf", &fv[1][1]);  fscanf(vu,"%lf", &fv[1][2]);
		fscanf(vu,"%lf", &fv[2][0]); fscanf(vu,"%lf", &fv[2][1]);  fscanf(vu,"%lf", &fv[2][2]);
		
		
		
		
		
		
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
			
			fscanf(vu,"%d", &wcbobb[i]);
			for (j=0; j<wcbobb[i]; j++) 
			{
				fscanf(vu,"%d", &cbobb[i][j][0]);
				fscanf(vu,"%d", &cbobb[i][j][1]);
				fscanf(vu,"%d", &cbobb[i][j][2]);
			}
			
			for (j=0; j<npatch[pid[i]]; j++)
			{
				for (jt=0; jt<3; jt++)  { fscanf(vu, "%lf", &paxe[i][j][jt]); }
			}
			
			
			for(l=0;l<9;l++)
			{
				fscanf(vu,"%lf", &li[i][l]);
				
			}
			
			lfp[i][0] = ami[i][0]*fv[0][0] + ami[i][1]*fv[1][0] + ami[i][2]*fv[2][0];
		    	lfp[i][1] = ami[i][0]*fv[0][1] + ami[i][1]*fv[1][1] + ami[i][2]*fv[2][1];
		    	lfp[i][2] = ami[i][0]*fv[0][2] + ami[i][1]*fv[1][2] + ami[i][2]*fv[2][2];
		
		}
		
		
		ec=0.0;
		for (i=0; i<3; i++)
		{
			for(j=i+1; j<3; j++)
			{
				ep=acos(fv[i][0]*fv[j][0]+fv[i][1]*fv[j][1]+fv[i][2]*fv[j][2]);
				///printf("%lf  %d  %d \n", ep, i, j);
				
				if (ep>piee/2.0) ep=piee-ep;
				
				ep=piee/2.0-ep;
				if (ep<0.0000001) ep=0.0;
				bdm=bleop/cos(ep);
				if (bdm>ec) ec=bdm;
			}
		
		}
		
		ble=ec;
		bleol=ble;
		
		
		
		for ( i=0; i<3; i++)
		{
			
			block[i]=(int)(floor(f[i]/ble));
			bl[i]=f[i]/(float)block[i];
			f[i]=bl[i]*(float)block[i];
			//printf("ptr %d \n", ptr);
		}
		
			
		boxfx=block[1]*block[2];
		boxfy=block[2];
		
		for(i=0; i<N; i++)	
	    	{	
			
			xi=ami[i][0];
			yi=ami[i][1];
			zi=ami[i][2];
	    		xk=(int)floor(xi/bl[0]);
	    		yk=(int)floor(yi/bl[1]);
	    		zk=(int)floor(zi/bl[2]);
	    		
	    		vtbox=boxfx*xk+boxfy*yk+zk;
		    	dum=wdr[vtbox];
	    	
			dri[vtbox][dum]=i;
			wdr[vtbox]=wdr[vtbox]+1;
			
			
		}
		
		
		for (i=0; i<3; i++)
		{
			fva[i][0]=fv[i][0]*f[i];
			fva[i][1]=fv[i][1]*f[i];
			fva[i][2]=fv[i][2]*f[i];
		}
		
		
		dfv = fv[0][0]*(fv[1][1]*fv[2][2]-fv[2][1]*fv[1][2]);
		dfv = dfv-fv[1][0]*(fv[0][1]*fv[2][2]-fv[0][2]*fv[2][1])+fv[2][0]*(fv[0][1]*fv[1][2]-fv[0][2]*fv[1][1]);
		
		fvi[0][0]=(fv[2][2]*fv[1][1]-fv[2][1]*fv[1][2])/dfv;
		fvi[1][0]=-(fv[0][1]*fv[2][2]-fv[2][1]*fv[0][2])/dfv;
		fvi[2][0]=(fv[0][1]*fv[1][2]-fv[1][1]*fv[0][2])/dfv;

		fvi[0][1]=-(fv[1][0]*fv[2][2]-fv[2][0]*fv[1][2])/dfv;
		fvi[1][1]=(fv[0][0]*fv[2][2]-fv[2][0]*fv[0][2])/dfv;
		fvi[2][1]=-(fv[0][0]*fv[1][2]-fv[1][0]*fv[0][2])/dfv;

		fvi[0][2]=(fv[1][0]*fv[2][1]-fv[1][1]*fv[2][0])/dfv;
		fvi[1][2]=-(fv[0][0]*fv[2][1]-fv[2][0]*fv[0][1])/dfv;
		fvi[2][2]=(fv[0][0]*fv[1][1]-fv[1][0]*fv[0][1])/dfv;
		
		
		fvo[0][0]=fv[0][0]; fvo[0][1]=fv[0][1]; fvo[0][2]=fv[0][2];
		fvo[1][0]=fv[1][0]; fvo[1][1]=fv[1][1]; fvo[1][2]=fv[1][2];
		fvo[2][0]=fv[2][0]; fvo[2][1]=fv[2][1]; fvo[2][2]=fv[2][2];
		
		fvio[0][0]=fvi[0][0]; fvio[0][1]=fvi[0][1]; fvio[0][2]=fvi[0][2];
		fvio[1][0]=fvi[1][0]; fvio[1][1]=fvi[1][1]; fvio[1][2]=fvi[1][2];
		fvio[2][0]=fvi[2][0]; fvio[2][1]=fvi[2][1]; fvio[2][2]=fvi[2][2];
		
		fvao[0][0]=fva[0][0]; fvao[0][1]=fva[0][1]; fvao[0][2]=fva[0][2];
		fvao[1][0]=fva[1][0]; fvao[1][1]=fva[1][1]; fvao[1][2]=fva[1][2];
		fvao[2][0]=fva[2][0]; fvao[2][1]=fva[2][1]; fvao[2][2]=fva[2][2];
		
		volo=fabs(fva[0][0]*(fva[1][1]*fva[2][2]-fva[2][1]*fva[1][2])-fva[0][1]*(fva[1][0]*fva[2][2]-fva[2][0]*fva[1][2])+fva[0][2]*(fva[1][0]*fva[2][1]-fva[2][0]*fva[1][1]));
		
		
		
		cow=1;
		
		printf("Initial configuration have been copied from the file and implemented into the system successfully!!!!!!! \n \n  ");
		
			//printf("Go and yourself\n");
		
		fclose(vu);
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		cfu=fopen(filename,"w");
		
		fclose(cfu);
		
	
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
     	
     	

	// Randomization ends here!






	//Preparing for simulation!
     
     
     	
     	FILE* lwu;
     	lwu=fopen("Diffusion_coefficient_el_l.txt","w");
     	fclose (lwu);
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
     
     	
     	
     	flim=3.0*ble+0.1;
     
     	
     	
     	
	kd=0.0;
	istep=1.0/(step*step);
	vpstep=0.01;
	vstep=step;
	rvstep=vstep*pow(3.0, 1.0/2.0);
     	
	lsd=10000;
	iltt=(int)floor((float)tsim/(float)lsd);
	rdd=(float)lsd/10.0;
	kd=0.0;
	istep=1.0/(step*step);
	
	//cow=0;
	titi=0;
	ri=0;
	nlies=2*N;
	pnlies=nlies+1;
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
	
	printf("Data is being printed after number of simulation time-    ");
	printf("%d \n", lsd);
	printf("\n");
	printf("\n");

	printf("Maximum possible data set with the given simulation time, is given by-   ");
	printf("%d\n",iltt);
	printf("\n");
	printf("\n");

	printf("\n");
	printf("\n");
	printf("\n");
	
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
	
	
	vpress=0.05;
	tsitt=10;
    	stsitt=10;
	
	

	
	///printf("%lf  %lf  \n", ble, bleol);
	
	
	
	//Simulation starts here!
	
	for(sitt=sitts; sitt<iltt; sitt++)
	{	
		
		
		
		if(sitt==tsitt)
		{
			
			vpress=vpress*5;
			if (vpress> press)
			{
				vpress=press;		
			}
			else
			{
				tsitt=stsitt+tsitt;
			}
		}
		
		
		sjet=0;
		pjet=0;
		tsjet=0;
		tpjet=0;
		latrf=10;
		
		
		/////printf("Skin 0 \n");
		for(ittt=0; ittt<lsd; ittt++)
		{	
		
			jet=0;
			
			vc=0;
			ptt=N+1;
			
			//printf("Skin 0 \n");
			
			for(mo=0; mo<pnlies; mo++)
			{
				
				
				
            			
            			if (mo<nlies)
            			{
		    			k=rand()  % N;
		    			
		    			boss=drand48();
		    			
		    			//if (cow==0 ) boss=0.4;
		    			rpid1=pid[k];
		    			
		    			
		    		}
		    		else if(mo==nlies)
		    		{
		    		
		    		
		    			if (vc==0)
		    			{
		    				tpjet=tpjet+1;
		    				vc=1;
		    				k=0;
		    				fer=fe;
		    				trims=2;
		    				mo=mo-1;
		    				//rd=drand48();
		    				
		    				
		    				boss=0.4;
		    				
		    				if (bshape_p==0 )
		    				{	
		    					fer=f[0];
		    					voln=exp(log(volo)+(drand48()-0.5)*vpstep);
		    					fe=pow(voln,0.33333333);
		    					fva[0][0]=fe;
		    					fva[1][1]=fe;
		    					fva[2][2]=fe;
		    					f[0]=fe;
		    					f[1]=fe;
		    					f[2]=fe;
		    					fr[0]=fe/fer;
		    					fr[1]=fr[0];
		    					fr[2]=fr[0];
		    				}
		    				else
		    				{
		    					///printf("Buddha was born to decore \n");
		    					jtt=0;
		    					while(jtt==0)
		    					{
			    					jtt=1;
			    					i=rand() % 3;
			    					voln=exp(log(volo)+(drand48()-0.5)*vpstep);
			    					fe=(voln/volo-1.0)*f[i];
			    					j=rand() % 3;
			    					fva[i][j]=fva[i][j]+fe;
			    					
			    					fer=f[i];
								f[i]=sqrt(fva[i][0]*fva[i][0]+fva[i][1]*fva[i][1]+fva[i][2]*fva[i][2]);
		    						
		    						fv[i][0]=fva[i][0]/f[i];
		    						fv[i][1]=fva[i][1]/f[i];
		    						fv[i][2]=fva[i][2]/f[i];
		    						
		    						fr[0]=1.0;
		    						fr[1]=1.0;
		    						fr[2]=1.0;
			    					fr[i]=f[i]/fer;
			    					
			    					dfv = fv[0][0]*(fv[1][1]*fv[2][2]-fv[2][1]*fv[1][2]);
								dfv = dfv-fv[1][0]*(fv[0][1]*fv[2][2]-fv[0][2]*fv[2][1])+fv[2][0]*(fv[0][1]*fv[1][2]-fv[0][2]*fv[1][1]);
			    				
				    				if  ( fabs(dfv) < 0.25  ||  f[i]<flim)
				    				{
				    					f[0]=f[0]/fr[0];
					    				f[1]=f[1]/fr[1];
					    				f[2]=f[2]/fr[2];
					    				
					    				fv[i][0]=fvo[i][0]; fv[i][1]=fvo[i][1]; fv[i][2]=fvo[i][2];
				    					
				    					fva[i][0]=fvao[i][0]; fva[i][1]=fvao[i][1]; fva[i][2]=fvao[i][2];
				    					fe=f[0];
	
									fr[0]=1.0;
					    				fr[1]=1.0;
					    				fr[2]=1.0;
					    				
				    					jtt=0;
				    				}
				    			}
				    			//printf("%lf \n", dfv);
				    			
				    			voln=fabs(fva[0][0]*(fva[1][1]*fva[2][2]-fva[2][1]*fva[1][2])-fva[0][1]*(fva[1][0]*fva[2][2]-fva[2][0]*fva[1][2])+fva[0][2]*(fva[1][0]*fva[2][1]-fva[2][0]*fva[1][1]));
		    				
							fvi[0][0]=(fv[2][2]*fv[1][1]-fv[2][1]*fv[1][2])/dfv;
							fvi[1][0]=-(fv[0][1]*fv[2][2]-fv[2][1]*fv[0][2])/dfv;
							fvi[2][0]=(fv[0][1]*fv[1][2]-fv[1][1]*fv[0][2])/dfv;

							fvi[0][1]=-(fv[1][0]*fv[2][2]-fv[2][0]*fv[1][2])/dfv;
							fvi[1][1]=(fv[0][0]*fv[2][2]-fv[2][0]*fv[0][2])/dfv;
							fvi[2][1]=-(fv[0][0]*fv[1][2]-fv[1][0]*fv[0][2])/dfv;

							fvi[0][2]=(fv[1][0]*fv[2][1]-fv[1][1]*fv[2][0])/dfv;
							fvi[1][2]=-(fv[0][0]*fv[2][1]-fv[2][0]*fv[0][1])/dfv;
							fvi[2][2]=(fv[0][0]*fv[1][1]-fv[1][0]*fv[0][1])/dfv;
							
							//
							
		    				
		    				}
		    				
		    				
		    				///printf("raaaam2  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2], dfv);
		    				
			    			///printf("ptr %d \n", ptr);
		    				for(i=0; i<N; i++)	
					    	{	
					    		xk=(int)floor(ami[i][0]/bl[0]);
					    		yk=(int)floor(ami[i][1]/bl[1]);
					    		zk=(int)floor(ami[i][2]/bl[2]);
					    		
							
							
							vtbox=boxfx*xk+boxfy*yk+zk;
		    					dum=wdr[vtbox];
							
							for (j=0; j<dum; j++)
							{
								dri[vtbox][j]=6000;
							
							}
							wdr[vtbox]=0;
							
							ami[i][0]=ami[i][0]*fr[0];
		    					ami[i][1]=ami[i][1]*fr[1];
		    					ami[i][2]=ami[i][2]*fr[2];
		    					
		    					
		    					
		    					lfp[i][0] = ami[i][0]*fv[0][0] + ami[i][1]*fv[1][0] + ami[i][2]*fv[2][0];
		    					lfp[i][1] = ami[i][0]*fv[0][1] + ami[i][1]*fv[1][1] + ami[i][2]*fv[2][1];
		    					lfp[i][2] = ami[i][0]*fv[0][2] + ami[i][1]*fv[1][2] + ami[i][2]*fv[2][2];
		    					
		    					
		    					
		    					
		    					for (j=0; j<wbobb[i]; j++)	{bobb[i][j][0]=6000;  bobb[i][j][1]=6000;  bobb[i][j][2]=6000;}
							wbobb[i]=0;
							
							for (j=0; j<wtnrd[i]; j++)	{tnrdi[i][j][0]=6000;  tnrdi[i][j][1]=13; }
							wtnrd[i]=0;
							
							for (j=0; j<wcbobb[i]; j++)	{bobb[i][j][0]=cbobb[i][j][0];  bobb[i][j][1]=cbobb[i][j][1]; bobb[i][j][2]=cbobb[i][j][2];}
							wbobb[i]=wcbobb[i];
							
							for (j=0; j<wnrd[i]; j++)	{tnrdi[i][j][0]=nrdi[i][j][0];  tnrdi[i][j][1]=nrdi[i][j][1];}
							wtnrd[i]=wnrd[i];
						}
		    				
		    				
		    				ec=0.0;
						for (i=0; i<3; i++)
						{
							for(j=i+1; j<3; j++)
							{
								ep=acos(fv[i][0]*fv[j][0]+fv[i][1]*fv[j][1]+fv[i][2]*fv[j][2]);
								///printf("%lf  %d  %d \n", ep, i, j);
								
								if (ep<piee/2.0) ep=piee-ep;
								
								ep=ep/2.0;
								if (ep<0.0000001) ep=0.0;
								bdm=bleop/(roott*cos(ep));
								if (bdm>ec) ec=bdm;
							}
						
						}
						
		    				ble=ec;
		    				///printf("at first site %lf  %lf \n", ble, ep);
		    				
		    				for (i=0; i<3; i++)
			    			{
							block[i]=(int)(floor(f[i]/ble));
							bl[i]=f[i]/(float)block[i];
							f[i]=bl[i]*(float)block[i];
							//printf("ptr %d \n", ptr);
						}
						
						boxfx=block[1]*block[2];  boxfy=block[2];
						//printf("budhdha was born to decore\n");
						////printf("Ptr________------ %d   %lf  %lf  %lf  %lf\n", ptr, fr[0], fr[1], fr[2], voln);
						
						
						
						//printf("budhdha was born to decore:    %ld   %ld\n", (long int)(boxfx*block[0]+boxfy*block[1]+block[2]), tbox);
						
						
						for(i=0; i<N; i++)	
					    	{	
							
							///printf("budhdha was born to decore %d\n", i);
					    		xk=(int)floor(ami[i][0]/bl[0]);
					    		yk=(int)floor(ami[i][1]/bl[1]);
					    		zk=(int)floor(ami[i][2]/bl[2]);
					    		
							//	printf("budhdha was born to decore %d  %d   %d  %d   %d\n", i, xk,  yk , zk, blockd);
							//	printf("budhdha was born to decore %d\n", i);						    		
					    		
					    		///if (xk>70 || yk>70 || zk>70 ||  xk<0 ||  yk<0|| zk<0 ) {printf("YOU ARE doomed  %d   %d   %d \n", xk, yk, zk);  }
					    		
					    		vtbox=boxfx*xk+boxfy*yk+zk;
		    					dum=wdr[vtbox];
					    		
							dri[vtbox][dum]=i;
							wdr[vtbox]=dum+1;
							
						}
						
						//printf("ptr hi h   %d\n", ptr);
						
					}
					else if (vc==1)
					{
						
						k=k+1;
			    			mo=mo-1;
			    			
			    			//printf("%d\n", k);
			    			if (k==N && trims==2)
			    			{
			    				vc=2;
			    				mo=mo+1;
			    				ene=0.0;
							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wbobb[i]; j++)
								{
									rpid2=pid[bobb[i][j][0]];
									ene=ene+ur[patches[rpid1][bobb[i][j][1]]][patches[rpid2][bobb[i][j][2]]];
								
								}
								
							}
							
							
							
							
							bdm=0.0;
							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wcbobb[i]; j++)
								{
									rpid2=pid[cbobb[i][j][0]];
									bdm=bdm+ur[patches[rpid1][cbobb[i][j][1]]][patches[rpid2][cbobb[i][j][2]]];
								
								}
								
							}

							
							ene=(bdm-ene)/2.0;
							
							
							bf=exp(-(ene+vpress*(voln-volo)-(float)(N+1)*log(voln/volo)));
							
							
							
							if (drand48()>bf)
							{
			    					trims=0;
			    					
			    					//printf("It was caused here !!!! %d  %lf  %lf  %lf  %lf  %lf   %lf   %lf\n", ittt, volo, voln, voln/volo, ene, bf, voln-volo, (float)N/voln);	
							}
							else
							{
							
								fvo[0][0]=fv[0][0]; fvo[0][1]=fv[0][1]; fvo[0][2]=fv[0][2];
		    						fvo[1][0]=fv[1][0]; fvo[1][1]=fv[1][1]; fvo[1][2]=fv[1][2];
		    						fvo[2][0]=fv[2][0]; fvo[2][1]=fv[2][1]; fvo[2][2]=fv[2][2];
		    				
		    						fvio[0][0]=fvi[0][0]; fvio[0][1]=fvi[0][1]; fvio[0][2]=fvi[0][2];
		    						fvio[1][0]=fvi[1][0]; fvio[1][1]=fvi[1][1]; fvio[1][2]=fvi[1][2];
		    						fvio[2][0]=fvi[2][0]; fvio[2][1]=fvi[2][1]; fvio[2][2]=fvi[2][2];
		    				
		    						fvao[0][0]=fva[0][0]; fvao[0][1]=fva[0][1]; fvao[0][2]=fva[0][2];
		    						fvao[1][0]=fva[1][0]; fvao[1][1]=fva[1][1]; fvao[1][2]=fva[1][2];
		    						fvao[2][0]=fva[2][0]; fvao[2][1]=fva[2][1]; fvao[2][2]=fva[2][2];
		    				
		    						volo=voln;
								bleol=ble;
								pjet=pjet+1;	
					    			k=rand()  % N;
				    				rpid1=pid[k];
					    			boss=drand48();	
							}
			    				
						}
						
						
						if(trims==0)
						{
							///printf("This is the faulty line %d  %d  %lf  %lf  %lf\n",  k, ittt, voln, volo, voln/volo);
							
							/*if (k<N)
							{
								printf("%d\n", wtnrd[k]);
								for (i=0; i<wtnrd[k]; i++)
								{
									v1=ami[k][0]-ami[tnrdi[k][i][0]][0];
									v2=ami[k][1]-ami[tnrdi[k][i][0]][1];
									v3=ami[k][2]-ami[tnrdi[k][i][0]][2];
									bdm=sqrt(v1*v1+v2*v2+v3*v3);
									printf("%d   %lf  %d\n", k, bdm, tnrdi[k][i][0]);
								
								}
							}*/
							
							
							vc=2;
							mo=mo+1;
				    			
				    			f[0]=f[0]/fr[0];
				    			f[1]=f[1]/fr[1];
				    			f[2]=f[2]/fr[2];
				    			
				    			fv[0][0]=fvo[0][0]; fv[0][1]=fvo[0][1]; fv[0][2]=fvo[0][2];
			    				fv[1][0]=fvo[1][0]; fv[1][1]=fvo[1][1]; fv[1][2]=fvo[1][2];
			    				fv[2][0]=fvo[2][0]; fv[2][1]=fvo[2][1]; fv[2][2]=fvo[2][2];
			    				
			    				fvi[0][0]=fvio[0][0]; fvi[0][1]=fvio[0][1]; fvi[0][2]=fvio[0][2];
			    				fvi[1][0]=fvio[1][0]; fvi[1][1]=fvio[1][1]; fvi[1][2]=fvio[1][2];
			    				fvi[2][0]=fvio[2][0]; fvi[2][1]=fvio[2][1]; fvi[2][2]=fvio[2][2];
			    				
			    				fva[0][0]=fvao[0][0]; fva[0][1]=fvao[0][1]; fva[0][2]=fvao[0][2];
			    				fva[1][0]=fvao[1][0]; fva[1][1]=fvao[1][1]; fva[1][2]=fvao[1][2];
			    				fva[2][0]=fvao[2][0]; fva[2][1]=fvao[2][1]; fva[2][2]=fvao[2][2];
				    			
				    			
				    			
				    			for(i=0; i<N; i++)	
						    	{	
						    		xk=(int)floor(ami[i][0]/bl[0]);
						    		yk=(int)floor(ami[i][1]/bl[1]);
						    		zk=(int)floor(ami[i][2]/bl[2]);
								
								vtbox=boxfx*xk+boxfy*yk+zk;
		    						dum=wdr[vtbox];
							
								for (j=0; j<dum; j++) dri[vtbox][j]=6000;
								wdr[vtbox]=0;
								
								
								
								ami[i][0]=ami[i][0]/fr[0];
			    					ami[i][1]=ami[i][1]/fr[1];
			    					ami[i][2]=ami[i][2]/fr[2];
			    					
			    					lfp[i][0] = ami[i][0]*fv[0][0] + ami[i][1]*fv[1][0] + ami[i][2]*fv[2][0];
		    						lfp[i][1] = ami[i][0]*fv[0][1] + ami[i][1]*fv[1][1] + ami[i][2]*fv[2][1];
		    						lfp[i][2] = ami[i][0]*fv[0][2] + ami[i][1]*fv[1][2] + ami[i][2]*fv[2][2];
			    					
							}
				    			ble=bleol;
				    			
			    				for (i=0; i<3; i++)
				    			{
								block[i]=(int)(floor(f[i]/ble));
								bl[i]=f[i]/(float)block[i];
								f[i]=bl[i]*(float)block[i];
								//printf("ptr %d \n", ptr);
							}
							
							boxfx=block[1]*block[2];  boxfy=block[2];
							
							for(i=0; i<N; i++)	
						    	{	
								
						    		xk=(int)floor(ami[i][0]/bl[0]);
						    		yk=(int)floor(ami[i][1]/bl[1]);
						    		zk=(int)floor(ami[i][2]/bl[2]);
						    		
						    		vtbox=boxfx*xk+boxfy*yk+zk;
		    						dum=wdr[vtbox];
								
								dri[vtbox][dum]=i;
								wdr[vtbox]=dum+1;
								
								
								
								for (j=0; j<wcbobb[i]; j++)	{cbobb[i][j][0]=6000;  cbobb[i][j][1]=6000;  cbobb[i][j][2]=6000;}
								wcbobb[i]=0;
							
								for (j=0; j<wnrd[i]; j++)	{nrdi[i][j][0]=6000;  nrdi[i][j][1]=13; }
								wnrd[i]=0;
							
								for (j=0; j<wbobb[i]; j++)	{cbobb[i][j][0]=bobb[i][j][0];  cbobb[i][j][1]=bobb[i][j][1]; cbobb[i][j][2]=bobb[i][j][2];}
								wcbobb[i]=wbobb[i];
							
								for (j=0; j<wtnrd[i]; j++)	{nrdi[i][j][0]=tnrdi[i][j][0];  nrdi[i][j][1]=tnrdi[i][j][1]; }
								wnrd[i]=wtnrd[i];
												
							}
		  
		    					fe=f[0];
		    					voln=volo;
		    					
		    					
				    			
				    			k=rand()  % N;
			    				rpid1=pid[k];
				    			boss=drand48();
						}
						
						boss=0.4;
						
						rpid1=pid[k];
						
						//printf("Buddha was born to decore the world\n");
		    			}
		    		
		    			
		    		}
		    		
            			
            			
            			
            				
            			mxi=ami[k][0];
	    			myi=ami[k][1];
	    			mzi=ami[k][2];
				
	    			l=(int)floor(mxi/bl[0]);
	    			m=(int)floor(myi/bl[1]);
	    			n=(int)floor(mzi/bl[2]);
				//printf("raamtlllllllllllllllll%d\n",k);
				
				wpbobb=0;
				for(i=0; i<E; i++)	
				{
					psbobb[i][0]=6000;
					psbobb[i][1]=6000;
					psbobb[i][2]=6000;
				}
		    		
				//printf("Flag 1\n");
            			if (boss >= 0.5)
				{
					//printf("hey rot!!!\n");
					ra1=drand48();
					theta=acos(1.0-(2.0*ra1));
					//ra2=drand48();
					ra2 = JKISS() / 4294967296.0;
	       				phi=tpiee*ra2;
					
					xf=sin(theta)*cos(phi);
					yf=sin(theta)*sin(phi);
					zf=cos(theta);
					
					
					theta=rvstep*drand48();
					
					
					
					
					rdli[0]=li[k][0];
					rdli[1]=li[k][1];
					rdli[2]=li[k][2];
					rdli[3]=li[k][3];
					rdli[4]=li[k][4];
					rdli[5]=li[k][5];
					rdli[6]=li[k][6];
					rdli[7]=li[k][7];
					rdli[8]=li[k][8];
					
					
					
					q1=cos(theta/2.0);
					q2=xf*sin(theta/2.0);
					q3=yf*sin(theta/2.0);
					q4=zf*sin(theta/2.0);
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
					
					dom0=rdom[0]*q11+rdom[1]*q12+rdom[2]*q13;
					dom1=rdom[0]*q21+rdom[1]*q22+rdom[2]*q23;
					dom2=rdom[0]*q31+rdom[1]*q32+rdom[2]*q33;
					
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
					
					//printf("Hey rot 2\n");
					
					
					rpid1=pid[k];
					for (i=0; i<npatch[rpid1]; i++)
					{	
						tpaxe[i][0]=paxe[k][i][0]*q11+paxe[k][i][1]*q12+paxe[k][i][2]*q13;
						tpaxe[i][1]=paxe[k][i][0]*q21+paxe[k][i][1]*q22+paxe[k][i][2]*q23;
						tpaxe[i][2]=paxe[k][i][0]*q31+paxe[k][i][1]*q32+paxe[k][i][2]*q33;
					}
					
					
					
					//printf("Hey rot 2\n");
					
					
							
						
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
						pdx=ami[puc][0]-mxi+(float)tnkx[jtt][0]*f[0];
						pdy=ami[puc][1]-myi+(float)tnkx[jtt][1]*f[1];
						pdz=ami[puc][2]-mzi+(float)tnkx[jtt][2]*f[2];
						
						re1 = pdx*fv[0][0] + pdy*fv[1][0] + pdz*fv[2][0];
						re2 = pdx*fv[0][1] + pdy*fv[1][1] + pdz*fv[2][1];
						re3 = pdx*fv[0][2] + pdy*fv[1][2] + pdz*fv[2][2];
						
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
							
							dive=0.0;
							bvs=0;
							bvsi=0;
							divea=0.0;
							divei=0.0;
							tri=0;
							
							//printf("Here  is the rot 1\n");
							
							
							
							for(i=0; i<npatch[rpid2]; i++)
							{
								bdm=-(paxe[puc][i][0]*re1+paxe[puc][i][1]*re2+paxe[puc][i][2]*re3)/mrp;
								if ( bdm > mom[patches[rpid2][i]] && bdm < wom[patches[rpid2][i]] )
								{	
									tri=1;
									bpc2=i;
									break;
								}
							}
						//	printf("%d \n", bpc2);
								
							if (tri==1)
							{
								for (j=0; j <npatch[rpid1]; j++)
								{
									if (pop[patches[rpid2][bpc2]][patches[rpid1][j]]==1)
									{
										
										bdm=(tpaxe[j][0]*re1+tpaxe[j][1]*re2+tpaxe[j][2]*re3)/mrp;
										if ( bdm > mom[patches[rpid1][j]]   &&  bdm < wom[patches[rpid1][j]] )
										{
												dive=epsi;
												divea=epsi2[patches[rpid1][j]][patches[rpid2][bpc2]];
												
												bpc1=j;
												
												bvs=1;
												break;
											
												
										}										
									}
								
								}
							}
							
							
							
							//printf("Here  is the rot 2\n");
							
							tri=0;
							
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
							
									check=1;
									jet=jet+1;
									
								}
							
								else
								{	
									ep=ep/mr;
									
									ec=ep*(mrp-divea)*(mrp-divea);
									if (bvs==1 && ec<1.0)
									{
										
										
										psbobb[wpbobb][0]=puc;
										psbobb[wpbobb][1]=bpc1;
										psbobb[wpbobb][2]=bpc2;
										wpbobb=wpbobb+1;
										
										
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
	 				
	 				
	 				//printf("Hey rot 3\n");
	 				
	 					
					
							
					
					if (check==0) 
					{	
						
						
						
						ene=0.0;
						for (it=0; it<wcbobb[k]; it++)
						{
							rpid2=pid[cbobb[k][it][0]];
							ene=ene+ur[patches[rpid1][cbobb[k][it][1]]][patches[rpid2][cbobb[k][it][2]]];
							
						}	
						
						
						
						bdm=0.0;
						for (i=0; i < wpbobb; i++)
						{
							
							rpid2=pid[psbobb[i][0]];
							bdm=bdm+ur[patches[rpid1][psbobb[i][1]]][patches[rpid2][psbobb[i][2]]];
							
						}
						
						
						
						pre=(bdm-ene);
						
						
						if (pre<0.0) bdm=1.1;
						else bdm=exp(-pre);
						
						if (drand48()<bdm)
						check=0;
						else check=1;
					}	
						
	
	
	
	
	
				
			
					if (check==0) 
					{
						
						
						mum=wcbobb[k];
							
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobb[k][jt][0];
							drum=wcbobb[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobb[chum][sj][0]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobb[chum][tric][0]=cbobb[chum][(tric+1)][0];
										cbobb[chum][tric][1]=cbobb[chum][(tric+1)][1];
										cbobb[chum][tric][2]=cbobb[chum][(tric+1)][2];
										
									}
									break;
								}
							
							}
							wcbobb[chum]=drum-1;
							cbobb[k][jt][0]=6000;
						}
							
						wcbobb[k]=0;	
						for (jdt=0; jdt < wpbobb; jdt++)
						{	
							//printf("raam is ");
							td=psbobb[jdt][0];
							cbobb[k][jdt][0]=td;
							cbobb[k][jdt][1]=psbobb[jdt][1];
							cbobb[k][jdt][2]=psbobb[jdt][2];
								tum=wcbobb[td];
								cbobb[td][tum][0]=k;
								cbobb[td][tum][1]=psbobb[jdt][2];
								cbobb[td][tum][2]=psbobb[jdt][1];
							wcbobb[td]=tum+1;
							wcbobb[k]=wcbobb[k]+1;
						}
						
						
						
						
						for (jt=0; jt<9; jt++)
						{
							li[k][jt]=dli[jt];
							el[k][jt]=del[jt];
						}
						
						
						
						for (jt=0; jt<npatch[rpid1]; jt++)
						{
							paxe[k][jt][0]=tpaxe[jt][0];
							paxe[k][jt][1]=tpaxe[jt][1];
							paxe[k][jt][2]=tpaxe[jt][2];
							
						}
						
						
						
						lom[k][0]=lom[k][0]+dom0;
						lom[k][1]=lom[k][1]+dom1;
						lom[k][2]=lom[k][2]+dom2;
						
					}
					//printf("Flag rot\n");  	
				}
                           	
                		else
				{
					trim=2;
					if(vc==1)
					{
						dx=0.0;
						dy=0.0;
						dz=0.0;
					}
					else
					{
						theta=acos(1.0-(2.0*drand48()));
						ra2 = JKISS() / 4294967296.0;
		       				phi=tpiee*ra2;
		       				bdm=vstep*JKISS() / 4294967296.0;
						dx=bdm*sin(theta)*cos(phi);
						dy=bdm*sin(theta)*sin(phi);
						dz=bdm*cos(theta);
						tsjet=tsjet+1;
					}
					
					
					xi=dx+lfp[k][0];
					
					yi=dy+lfp[k][1];
					
					zi=dz+lfp[k][2];
					
					
					dxi=fvi[0][0]*xi+fvi[0][1]*yi+fvi[0][2]*zi;	
					dyi=fvi[1][0]*xi+fvi[1][1]*yi+fvi[1][2]*zi;	
					dzi=fvi[2][0]*xi+fvi[2][1]*yi+fvi[2][2]*zi;
					
				
					if (dxi >= f[0]) dxi=dxi-f[0];
					else if (dxi < 0.0) dxi=dxi+f[0];
					if (dyi >= f[1]) dyi=dyi-f[1];
					else if (dyi < 0.0) dyi=dyi+f[1];
					if (dzi >= f[2]) dzi=dzi-f[2];
					else if (dzi < 0.0) dzi=dzi+f[2];
	        			
	        			xi=fv[0][0]*dxi+fv[1][0]*dyi+fv[2][0]*dzi;	
					yi=fv[0][1]*dxi+fv[1][1]*dyi+fv[2][1]*dzi;	
					zi=fv[0][2]*dxi+fv[1][2]*dyi+fv[2][2]*dzi;
					
	        			xk=(int)floor(dxi/bl[0]);
	        			yk=(int)floor(dyi/bl[1]);
	        			zk=(int)floor(dzi/bl[2]);
					
					
					//printf("raaaam2  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf  \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2]);
					
					
					//printf("here is the problem first    !!!!!\n");	
	    				for (j=0; j<UN; j++)	{pit[j][0]=6000; pit[j][1]=13;}
					wpi=0;	
					
					
					///printf("Flag_trans 2  %d\n", mo);
					
					for(jt=0;jt<27; jt++)
					{
						///printf("Budhdha \n");
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
						
						
						
						vtbox=boxfx*nxk+boxfy*nyk+nzk;
		    				///printf("Budhdha 3333  %ld  %ld   %d   %d   %d   %d   %d  %d\n", vtbox, tbox, nxk, nyk, nzk, xk, yk, zk);
		    				
		    			
		    				
		    				dum=wdr[vtbox];
 						///printf("Budhdha 3333\n");
						
						for (jdt=0; jdt<dum; jdt++)
						{
							da=dri[vtbox][jdt];  
							
							///if (da>5000 ) printf("Buddha will resurrect from the corpses!!!!!!!!  %ld\n", vtbox);
							rpid2=pid[da];         
							pdx=ami[da][0]+(float)cpsh*f[0];
							pdy=ami[da][1]+(float)cph*f[1];
							pdz=ami[da][2]+(float)cp*f[2];
							
							xf = pdx*fv[0][0] + pdy*fv[1][0] + pdz*fv[2][0];
							yf = pdx*fv[0][1] + pdy*fv[1][1] + pdz*fv[2][1];
							zf = pdx*fv[0][2] + pdy*fv[1][2] + pdz*fv[2][2];
							
							jtt=(cpsh+1)*9+(cph+1)*3+cp+1;
							
							if (da != k) 
							{
				    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));
				    				bdm=dss[rpid1][rpid2];
				    				if (rijs <= dss[rpid1][rpid2])
				    				{
					    				
					    					pit[wpi][0]=da;
										pit[wpi][1]=jtt;
										wpi=wpi+1;
								                if (rijs <= ds[rpid1][rpid2])
										{
											if (vc==1) reserve=rijs;
											trim=1;
											break;
										}   
										/*else
										{
											pit[wpi][0]=da;
											pit[wpi][1]=jtt;
											wpi=wpi+1;
										}
										*/
										
									
								}
							}
							
						} 
						
						///printf("Budhdha 2222\n");
						if (trim==1) break;
						
					
					}
						
					///printf("Flag_trans 33\n");
					
					
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
							
							pdx=ami[puc][0]-dxi+(float)tnkx[jtt][0]*f[0];
							pdy=ami[puc][1]-dyi+(float)tnkx[jtt][1]*f[1];
							pdz=ami[puc][2]-dzi+(float)tnkx[jtt][2]*f[2];
							
							
							
							re1 = pdx*fv[0][0] + pdy*fv[1][0] + pdz*fv[2][0];
							re2 = pdx*fv[0][1] + pdy*fv[1][1] + pdz*fv[2][1];
							re3 = pdx*fv[0][2] + pdy*fv[1][2] + pdz*fv[2][2];
							
							
							
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
							
								lc=0;
								
								dive=0.0;
								bvs=0;
								bvsi=0;
								divea=0.0;
								divei=0.0;
								tri=0;
								
								
								for(i=0; i < npatch[rpid2]; i++)
								{
									
									bdm=-(paxe[puc][i][0]*re1+paxe[puc][i][1]*re2+paxe[puc][i][2]*re3)/mrp;
									if ( bdm > mom[patches[rpid2][i]] &&  bdm < wom[patches[rpid2][i]] )
									{	
										tri=1;
										bpc2=i;
										break;
									}
								}
									
								//printf("Trans sub1  %d   %d   %d \n", patches[rpid2][jt], rpid2, jt);
								if (tri==1)
								{
									for (j = 0; j < npatch[rpid1]; j++)
									{
										if (pop[patches[rpid2][bpc2]][patches[rpid1][j]]==1)
										{
											bdm=(paxe[k][j][0]*re1+paxe[k][j][1]*re2+paxe[k][j][2]*re3)/mrp;
											
											if (  bdm > mom[patches[rpid1][j]]  &&  bdm < wom[patches[rpid1][j]] )
											{
													
													
													dive=epsi;
													divea=epsi2[patches[rpid1][j]][patches[rpid2][bpc2]];
												
													bpc1=j;
												
													bvs=1;
													break;
													
											}										
										}
									
									}
								}
								
								
								
								//printf("Trans sub2\n");
								tri=0;
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
								
										trim=1;
										
									}
								
									else
									{	
										ep=ep/mr;
										
										ec=ep*(mrp-divea)*(mrp-divea);
										if (bvs==1 && ec<1.0)
										{
											
											
											psbobb[wpbobb][0]=puc;
											psbobb[wpbobb][1]=bpc1;
											psbobb[wpbobb][2]=bpc2;
											wpbobb=wpbobb+1;
												
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
					
					//printf("Flag_trans4\n");
					
					
					if (trim==2 && vc!=1)
					{
					
						
						ene=0.0;
						for (it=0; it<wcbobb[k]; it++)
						{
							rpid2=pid[cbobb[k][it][0]];
							ene=ene+ur[patches[rpid1][cbobb[k][it][1]]][patches[rpid2][cbobb[k][it][2]]];
							
						}	
						
						
						
						bdm=0.0;
						for (i=0; i < wpbobb; i++)
						{
							
							rpid2=pid[psbobb[i][0]];
							bdm=bdm+ur[patches[rpid1][psbobb[i][1]]][patches[rpid2][psbobb[i][2]]];
							//printf("Bonds are being formed %lf   %d   %d\n", bdm, patches[rpid1][psbobb[i][1]], patches[rpid2][psbobb[i][2]] );
							 
						}
						
						pre=(bdm-ene);
						
						if (pre<0.0) bdm=1.1; 
						else bdm=exp(-pre);
						
						
						if (drand48()<bdm)
						trim=2;
						else trim=1;
						
					}
					
					//printf("Flag trans   55\n");
					if (trim==2)
					{
						if (vc!=1) sjet=sjet+1;
						
						//printf("Flag trans   555\n");
					
						vtbox=boxfx*l+boxfy*m+n;
		    				dum=wdr[vtbox];

						for (j=0;j<dum;j++)
						{	
							
							if (dri[vtbox][j]==k)
							{
								for (tric=j; tric<dum; tric++)	dri[vtbox][tric]=dri[vtbox][tric+1];
									
								
								break;
							}
									
						}
						wdr[vtbox]=wdr[vtbox]-1;
						
						//printf("Flag trans   666\n");
						vtbox=boxfx*xk+boxfy*yk+zk;
		    				dum=wdr[vtbox];
		
						dri[vtbox][dum]=k;
						wdr[vtbox]=dum+1;
						
						lfp[k][0]=xi;
		    				ami[k][0]=dxi;
		    				lfp[k][1]=yi;
		    				ami[k][1]=dyi;
		    				lfp[k][2]=zi;
		    				ami[k][2]=dzi;
		    				
		    				//printf("Flag trans   666\n");
					
		    				
		    				
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
						
						
						
							
						
						
						
						mum=wcbobb[k];
							
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobb[k][jt][0];
							drum=wcbobb[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobb[chum][sj][0]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobb[chum][tric][0]=cbobb[chum][(tric+1)][0];
										cbobb[chum][tric][1]=cbobb[chum][(tric+1)][1];
										cbobb[chum][tric][2]=cbobb[chum][(tric+1)][2];
									}
									break;
								}
							
							}
							wcbobb[chum]=drum-1;
							cbobb[k][jt][0]=6000;
						}
							
						wcbobb[k]=0;	
						for (jdt=0; jdt < wpbobb; jdt++)
						{	
							//printf("raam is ");
							td=psbobb[jdt][0];
							cbobb[k][jdt][0]=td;
							cbobb[k][jdt][1]=psbobb[jdt][1];
							cbobb[k][jdt][2]=psbobb[jdt][2];
								tum=wcbobb[td];
								cbobb[td][tum][0]=k;
								cbobb[td][tum][1]=psbobb[jdt][2];
								cbobb[td][tum][2]=psbobb[jdt][1];
							wcbobb[td]=tum+1;
							wcbobb[k]=wcbobb[k]+1;
						}
							
						
					}
					
					else if(vc==1 && trim==1)
					{
					
						trims=0;  ///printf("\n\n\n Buddha %d %d  %d  %d \n\n\n", ittt, k, pit[wpi-1][0], trim);
						/*for (i=0; i<wpi; i++)
						{
							v1=lfp[pit[i][0]][0] - lfp[k][0] + f[0]*tnkx[pit[i][1]][0];
							v2=lfp[pit[i][0]][1] - lfp[k][1] + f[1]*tnkx[pit[i][1]][1];
							v3=lfp[pit[i][0]][2] - lfp[k][2] + f[2]*tnkx[pit[i][1]][2];
							bdm=(v1*v1+v2*v2+v3*v3);
							printf(" Hey there - %d   %d   %lf  %d  %lf  %lf    \n", k, pit[i][0], bdm, pit[i][1], dss[0][0], reserve);
							
							//printf("raaaam2  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf  \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2]);
							//printf("raaaam2  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf  \n", fv[0][0], fv[0][1], fv[0][2], fv[1][0], fv[1][1], fv[1][2], fv[2][0], fv[2][1], fv[2][2]);
						}*/
					}
					
					
					
					//printf("Flag trans 3\n");
					
					
	
				}
				
				
			//printf("snn750 %d\n",mo);
			}
			
			
			
			if (ittt==latrf && bshape_p==1)
			{
				///printf("Skin 1\n");
			        
				for(i=0; i<N; i++)	
			    	{	
			    		xk=(int)floor(ami[i][0]/bl[0]);
			    		yk=(int)floor(ami[i][1]/bl[1]);
			    		zk=(int)floor(ami[i][2]/bl[2]);
			    		
			    		vtbox=boxfx*xk+boxfy*yk+zk;
		    			dum=wdr[vtbox];
					for (j=0; j<dum; j++)
					{
						dri[vtbox][j]=6000;
					
					}
					wdr[vtbox]=0;
				}
				
				///////   Lattice reduction area !!!!!!!!!!!!!!
				
				w1=fva[0][1]*fva[1][2]-fva[1][1]*fva[0][2];
				w2=fva[1][0]*fva[0][2]-fva[0][0]*fva[1][2];
				w3=fva[0][0]*fva[1][1]-fva[1][0]*fva[0][1];
				
				w4=fva[1][1]*fva[2][2]-fva[2][1]*fva[1][2];
				w5=fva[2][0]*fva[1][2]-fva[1][0]*fva[2][2];
				w6=fva[1][0]*fva[2][1]-fva[2][0]*fva[1][1];
				
				w7=fva[2][1]*fva[0][2]-fva[0][1]*fva[2][2];
				w8=fva[0][0]*fva[2][2]-fva[2][0]*fva[0][2];
				w9=fva[2][0]*fva[0][1]-fva[0][0]*fva[2][1];
				
				ree[0]=w1*w1+w2*w2+w3*w3;
				ree[1]=w4*w4+w5*w5+w6*w6;
				ree[2]=w7*w7+w8*w8+w9*w9;
				
				ree[0]=sqrt(ree[0]);
				ree[1]=sqrt(ree[1]);
				ree[2]=sqrt(ree[2]);
				
				slredo=ree[0]+ree[1]+ree[2];
			
				tri=1;
			    	while (tri==1)
			    	{
					tri=0;
					
					for(i=0; i<3; i++)
					{
						
						for(j=i+1; j<3; j++)
						{
							l=i;
							m=j;
							for(it=0; it<2; it++)
							{
								pre=1.0;
								for(jt=0; jt<2; jt++)
								{
									rdli[0]=fva[0][0]; rdli[1]=fva[0][1]; rdli[2]=fva[0][2];			
			    						rdli[3]=fva[1][0]; rdli[4]=fva[1][1]; rdli[5]=fva[1][2];
			    						rdli[6]=fva[2][0]; rdli[7]=fva[2][1]; rdli[8]=fva[2][2];
									kt=m*3;
									
									rdli[kt]=fva[m][0]+pre*fva[l][0];
									rdli[kt+1]=fva[m][1]+pre*fva[l][1];
									rdli[kt+2]=fva[m][2]+pre*fva[l][2];
									
									
								
									w1=rdli[1]*rdli[5]-rdli[4]*rdli[2];
									w2=rdli[3]*rdli[2]-rdli[0]*rdli[5];
									w3=rdli[0]*rdli[4]-rdli[3]*rdli[1];
									w4=rdli[4]*rdli[8]-rdli[7]*rdli[5];
									w5=rdli[6]*rdli[5]-rdli[3]*rdli[8];
									w6=rdli[3]*rdli[7]-rdli[6]*rdli[4];
									w7=rdli[7]*rdli[2]-rdli[1]*rdli[8];
									w8=rdli[0]*rdli[8]-rdli[6]*rdli[2];
									w9=rdli[6]*rdli[1]-rdli[0]*rdli[7];
									
									ree[0]=w1*w1+w2*w2+w3*w3;
									ree[1]=w4*w4+w5*w5+w6*w6;
									ree[2]=w7*w7+w8*w8+w9*w9;
									
									ree[0]=sqrt(ree[0]);
									ree[1]=sqrt(ree[1]);
									ree[2]=sqrt(ree[2]);
									
									slred=ree[0]+ree[1]+ree[2];
								
					    				
									if(slred<slredo)
									{
										
										u1=rdli[0]; u2=rdli[1]; u3=rdli[2];			
					    					u4=rdli[3]; u5=rdli[4]; u6=rdli[5];
					    					u7=rdli[6]; u8=rdli[7]; u9=rdli[8];
					 
					    					slredo=slred;
					    					
					    					tri=1;
					    					
									}
									pre=-1;
								}
								l=j;
								m=i;
							}
						}
					}
					
					if(tri==1)
					{
						fva[0][0]=u1;  fva[0][1]=u2;  fva[0][2]=u3;
						fva[1][0]=u4;  fva[1][1]=u5;  fva[1][2]=u6;
						fva[2][0]=u7;  fva[2][1]=u8;  fva[2][2]=u9;
						
						f[0]=sqrt(fva[0][0]*fva[0][0]+fva[0][1]*fva[0][1]+fva[0][2]*fva[0][2]);
						fv[0][0]=fva[0][0]/f[0]; fv[0][1]=fva[0][1]/f[0]; fv[0][2]=fva[0][2]/f[0];
						
						f[1]=sqrt(fva[1][0]*fva[1][0]+fva[1][1]*fva[1][1]+fva[1][2]*fva[1][2]);
						fv[1][0]=fva[1][0]/f[1]; fv[1][1]=fva[1][1]/f[1]; fv[1][2]=fva[1][2]/f[1]; 
						
						f[2]=sqrt(fva[2][0]*fva[2][0]+fva[2][1]*fva[2][1]+fva[2][2]*fva[2][2]);
						fv[2][0]=fva[2][0]/f[2]; fv[2][1]=fva[2][1]/f[2]; fv[2][2]=fva[2][2]/f[2]; 
						
						
						
						dfv = fv[0][0]*(fv[1][1]*fv[2][2]-fv[2][1]*fv[1][2])-fv[1][0]*(fv[0][1]*fv[2][2]-fv[0][2]*fv[2][1])+fv[2][0]*(fv[0][1]*fv[1][2]-fv[0][2]*fv[1][1]);
						
						fvi[0][0]=(fv[2][2]*fv[1][1]-fv[2][1]*fv[1][2])/dfv;
						fvi[1][0]=-(fv[0][1]*fv[2][2]-fv[2][1]*fv[0][2])/dfv;
						fvi[2][0]=(fv[0][1]*fv[1][2]-fv[1][1]*fv[0][2])/dfv;

						fvi[0][1]=-(fv[1][0]*fv[2][2]-fv[2][0]*fv[1][2])/dfv;
						fvi[1][1]=(fv[0][0]*fv[2][2]-fv[2][0]*fv[0][2])/dfv;
						fvi[2][1]=-(fv[0][0]*fv[1][2]-fv[1][0]*fv[0][2])/dfv;

						fvi[0][2]=(fv[1][0]*fv[2][1]-fv[1][1]*fv[2][0])/dfv;
						fvi[1][2]=-(fv[0][0]*fv[2][1]-fv[2][0]*fv[0][1])/dfv;
						fvi[2][2]=(fv[0][0]*fv[1][1]-fv[1][0]*fv[0][1])/dfv;

					    	///printf("raaaam3  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2], dfv);
					    	///printf("raaaam3  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fv[0][0], fv[0][1], fv[0][2], fv[1][0], fv[1][1], fv[1][2], fv[2][0], fv[2][1], fv[2][2], dfv);
						for (it=0; it<N; it++)
						{
							
							xi = lfp[it][0];
							yi = lfp[it][1];
							zi = lfp[it][2];
						
							dxi = fvi[0][0]*xi+fvi[0][1]*yi+fvi[0][2]*zi;		
							dyi = fvi[1][0]*xi+fvi[1][1]*yi+fvi[1][2]*zi;
							dzi = fvi[2][0]*xi+fvi[2][1]*yi+fvi[2][2]*zi;
							
							xi=floor(dxi/f[0]);
							dxi=dxi-xi*f[0];
							yi=floor(dyi/f[1]);
							dyi=dyi-yi*f[1];	
							zi=floor(dzi/f[2]);
							dzi=dzi-zi*f[2];
							
							ami[it][0]=dxi;
							ami[it][1]=dyi;
							ami[it][2]=dzi;
							
							lfp[it][0] = fv[0][0]*dxi+fv[1][0]*dyi+fv[2][0]*dzi;	
							lfp[it][1] = fv[0][1]*dxi+fv[1][1]*dyi+fv[2][1]*dzi;	
							lfp[it][2] = fv[0][2]*dxi+fv[1][2]*dyi+fv[2][2]*dzi;
							
							
						}
						
						
					}
					
					
				
				}
				
				
				
				/////// Box alingnement towards the xy-plane !!!!!!!!!!!!!
				
				
				
				
				
				if (fabs(fv[0][0])>0.99999) 
				{
					dx=0.0;
					dy=0.0;
					dz=1.0;
					rdth=acos(fv[0][0]);
					if(fv[0][0]>=1.0) rdth=0.0;
					else if (fv[0][0]<=-1.0) rdth=piee;
				}
				else
				{
					dx=0.0;
					dy=-fv[0][2]/sqrt(fv[0][2]*fv[0][2]+fv[0][1]*fv[0][1]);
					dz=fv[0][1]/sqrt(fv[0][2]*fv[0][2]+fv[0][1]*fv[0][1]);
					rdth=acos(fv[0][0]);
				}
				///printf("%lf   %lf   %lf   %lf\n", dx, dy, dz, rdth);
				
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
				
				dli[0]=fv[0][0]*q11+fv[0][1]*q12+fv[0][2]*q13;
				dli[1]=fv[0][0]*q21+fv[0][1]*q22+fv[0][2]*q23;
				dli[2]=fv[0][0]*q31+fv[0][1]*q32+fv[0][2]*q33;
				
				dli[3]=fv[1][0]*q11+fv[1][1]*q12+fv[1][2]*q13;
				dli[4]=fv[1][0]*q21+fv[1][1]*q22+fv[1][2]*q23;
				dli[5]=fv[1][0]*q31+fv[1][1]*q32+fv[1][2]*q33;
				
				dli[6]=fv[2][0]*q11+fv[2][1]*q12+fv[2][2]*q13;
				dli[7]=fv[2][0]*q21+fv[2][1]*q22+fv[2][2]*q23;
				dli[8]=fv[2][0]*q31+fv[2][1]*q32+fv[2][2]*q33;
				
				fv[0][0]=dli[0];  fv[0][1]=dli[1];   fv[0][2]=dli[2];
				fv[1][0]=dli[3];  fv[1][1]=dli[4];   fv[1][2]=dli[5];
				fv[2][0]=dli[6];  fv[2][1]=dli[7];   fv[2][2]=dli[8];
				
				
				
				
				
				bdm=sqrt(fv[1][1]*fv[1][1]+fv[1][2]*fv[1][2]);
				
				
							
				rdth = acos(fv[1][1]/bdm);
				if (fv[1][2]>0.0) dx=fv[0][0];
				else dx=-fv[0][0];			
				dy=0.0;
				dz=0.0;
				
				
				
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
				
				dli[3]=fv[1][0]*q11+fv[1][1]*q12+fv[1][2]*q13;
				dli[4]=fv[1][0]*q21+fv[1][1]*q22+fv[1][2]*q23;
				dli[5]=fv[1][0]*q31+fv[1][1]*q32+fv[1][2]*q33;
				
				dli[6]=fv[2][0]*q11+fv[2][1]*q12+fv[2][2]*q13;
				dli[7]=fv[2][0]*q21+fv[2][1]*q22+fv[2][2]*q23;
				dli[8]=fv[2][0]*q31+fv[2][1]*q32+fv[2][2]*q33;
				
				fv[1][0]=dli[3];  fv[1][1]=dli[4];   fv[1][2]=dli[5];
				fv[2][0]=dli[6];  fv[2][1]=dli[7];   fv[2][2]=dli[8];
				
				fva[0][0]=f[0]*fv[0][0];  fva[0][1]=f[0]*fv[0][1];  fva[0][2]=f[0]*fv[0][2];
				fva[1][0]=f[1]*fv[1][0];  fva[1][1]=f[1]*fv[1][1];  fva[1][2]=f[1]*fv[1][2];
				fva[2][0]=f[2]*fv[2][0];  fva[2][1]=f[2]*fv[2][1];  fva[2][2]=f[2]*fv[2][2];
				
				dfv = fv[0][0]*(fv[1][1]*fv[2][2]-fv[2][1]*fv[1][2]) - fv[1][0]*(fv[0][1]*fv[2][2]-fv[0][2]*fv[2][1]) + fv[2][0]*(fv[0][1]*fv[1][2]-fv[0][2]*fv[1][1]);
				//printf("raaaam4  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2], dfv);		
				
				
				fvi[0][0] = (fv[2][2]*fv[1][1]-fv[2][1]*fv[1][2])/dfv;
				fvi[1][0] = -(fv[0][1]*fv[2][2]-fv[2][1]*fv[0][2])/dfv;
				fvi[2][0] = (fv[0][1]*fv[1][2]-fv[1][1]*fv[0][2])/dfv;

				fvi[0][1] = -(fv[1][0]*fv[2][2]-fv[2][0]*fv[1][2])/dfv;
				fvi[1][1] = (fv[0][0]*fv[2][2]-fv[2][0]*fv[0][2])/dfv;
				fvi[2][1] = -(fv[0][0]*fv[1][2]-fv[1][0]*fv[0][2])/dfv;

				fvi[0][2] = (fv[1][0]*fv[2][1]-fv[1][1]*fv[2][0])/dfv;
				fvi[1][2] = -(fv[0][0]*fv[2][1]-fv[2][0]*fv[0][1])/dfv;
				fvi[2][2] = (fv[0][0]*fv[1][1]-fv[1][0]*fv[0][1])/dfv;
				
				//printf("raaaam5  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2], dfv);
				
				
				for (it=0; it<N; it++)
				{
					
					lfp[it][0] = fv[0][0]*ami[it][0]+fv[1][0]*ami[it][1]+fv[2][0]*ami[it][2];	
					lfp[it][1] = fv[0][1]*ami[it][0]+fv[1][1]*ami[it][1]+fv[2][1]*ami[it][2];	
					lfp[it][2] = fv[0][2]*ami[it][0]+fv[1][2]*ami[it][1]+fv[2][2]*ami[it][2];
					
				}
				
				
				
				//////////////////////////           /////////////////////
				
				fvo[0][0]=fv[0][0]; fvo[0][1]=fv[0][1]; fvo[0][2]=fv[0][2];			
				fvo[1][0]=fv[1][0]; fvo[1][1]=fv[1][1]; fvo[1][2]=fv[1][2];
				fvo[2][0]=fv[2][0]; fvo[2][1]=fv[2][1]; fvo[2][2]=fv[2][2];
				
				fvio[0][0]=fvi[0][0]; fvio[0][1]=fvi[0][1]; fvio[0][2]=fvi[0][2];			
				fvio[1][0]=fvi[1][0]; fvio[1][1]=fvi[1][1]; fvio[1][2]=fvi[1][2];
				fvio[2][0]=fvi[2][0]; fvio[2][1]=fvi[2][1]; fvio[2][2]=fvi[2][2];
				
				fvao[0][0]=fva[0][0]; fvao[0][1]=fva[0][1]; fvao[0][2]=fva[0][2];			
				fvao[1][0]=fva[1][0]; fvao[1][1]=fva[1][1]; fvao[1][2]=fva[1][2];
				fvao[2][0]=fva[2][0]; fvao[2][1]=fva[2][1]; fvao[2][2]=fva[2][2];
				
				ec=0.0;
				for (i=0; i<3; i++)
				{
					for(j=i+1; j<3; j++)
					{
						ep=acos(fv[i][0]*fv[j][0]+fv[i][1]*fv[j][1]+fv[i][2]*fv[j][2]);
						
						
						
						if (ep<piee/2.0) ep=piee-ep;
								
						ep=ep/2.0;
						if (ep<0.0000001) ep=0.0;
						bdm=bleop/(roott*cos(ep));
						
						
						if (bdm>ec) ec=bdm;
					}
				
				}
				ble=ec;
				bleol=ble;
				
				
				for (i=0; i<3; i++)
	    			{
					block[i]=(int)(floor(f[i]/ble));
					bl[i]=f[i]/(float)block[i];
					f[i]=bl[i]*(float)block[i];
					
					//printf("ptr %d \n", ptr);
				}
				boxfx=block[1]*block[2];  boxfy=block[2];
				
				for(i=0; i<N; i++)	
			    	{	
					
			    		xk=(int)floor(ami[i][0]/bl[0]);
			    		yk=(int)floor(ami[i][1]/bl[1]);
			    		zk=(int)floor(ami[i][2]/bl[2]);
					
					vtbox=boxfx*xk+boxfy*yk+zk;
		    			dum=wdr[vtbox];
					
					dri[vtbox][dum]=i;
					wdr[vtbox]=dum+1;
					
					for (j=0; j<wnrd[i]; j++)
					{
						cp=0;
						cph=0;
						cpsh=0;
						
						if (ami[nrdi[i][j][0]][0]-ami[i][0]<-bl[0]) cpsh=1;
						else if(ami[nrdi[i][j][0]][0]-ami[i][0]>bl[0]) cpsh=-1;
						
						if (ami[nrdi[i][j][0]][1]-ami[i][1]<-bl[1]) cph=1;
						else if(ami[nrdi[i][j][0]][1]-ami[i][1]>bl[1]) cph=-1;
						
						if (ami[nrdi[i][j][0]][2]-ami[i][2]<-bl[2]) cp=1;
						else if(ami[nrdi[i][j][0]][2]-ami[i][2]>bl[2]) cp=-1;
						
						jtt=(cpsh+1)*9+(cph+1)*3+cp+1;
						nrdi[i][j][1]=jtt;
					
					}
					
					
				}
				latrf=latrf+10;	
				
				///printf("skin 22222222222222  %d\n", ittt);
			}
		
			
				
				
		}
		printf("Simulation is going on si %d\n", sitt+1);
		
		printf("Simulation properties are: vstep %lf and pstep %lf \n", vstep, vpstep);
		
		dfv = fv[0][0]*(fv[1][1]*fv[2][2]-fv[2][1]*fv[1][2]);
		dfv =fabs( dfv-fv[1][0]*(fv[0][1]*fv[2][2]-fv[0][2]*fv[2][1])+fv[2][0]*(fv[0][1]*fv[1][2]-fv[0][2]*fv[1][1]));
		printf("Thermodynamic parameters are-   Temp: %lf,   Current pressure:  %lf,   Target pressure:  %lf, Number fraction:  %lf\n", -1.0/ur[0][0], vpress, press, (float)N/volo);
	
	
		printf("Box properties Edge:  %lf, Unit shell volume:  %lf, Unit block length:  %lf\n\n\n", f[2], dfv, ble);
		///printf("raaaam3  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fv[0][0], fv[0][1], fv[0][2], fv[1][0], fv[1][1], fv[1][2], fv[2][0], fv[2][1], fv[2][2], dfv);
		///printf("raaaam3  %lf %lf %lf    %lf %lf %lf    %lf %lf %lf   %lf \n", fvi[0][0], fvi[0][1], fvi[0][2], fvi[1][0], fvi[1][1], fvi[1][2], fvi[2][0], fvi[2][1], fvi[2][2], dfv);
		bdm=(float)sjet/(float)tsjet;
		
		if (bdm<0.3) vstep=vstep*0.95;
		else if (bdm > 0.4)  vstep=vstep*1.05;
		
		
	
		bdm=(float)pjet/(float)tpjet;
		
		if (bdm<0.3) vpstep=vpstep*0.95;
		else if (bdm > 0.4)  vpstep=vpstep*1.05;
		
		
		
		if(vstep>0.1) vstep=0.1;
		else if (vstep<0.001) vstep=0.001;
		
		if(vpstep>0.1) vpstep=0.1;
		else if (vpstep<0.001) vpstep=0.001;
		
		rvstep=pow(3.0, 1.0/2.0)*vstep;
		
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitt)/100.0));

		FILE* cfu;
		cfu=fopen(filename,"a");
		
			
		fprintf(cfu,"%lf\n", f[0]);
		fprintf(cfu,"%lf\n", f[1]);
		fprintf(cfu,"%lf\n", f[2]);
		
		fprintf(cfu,"%lf   %lf    %lf\n", fv[0][0], fv[0][1], fv[0][2]);
		fprintf(cfu,"%lf   %lf    %lf\n", fv[1][0], fv[1][1], fv[1][2]);
		fprintf(cfu,"%lf   %lf    %lf\n", fv[2][0], fv[2][1], fv[2][2]);
		
		for (i=0; i<N; i++)
		{
			fprintf(cfu,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
			
			fprintf(cfu, "%d\n", wcbobb[i]);
			for (j=0; j<wcbobb[i]; j++) 
			{
				fprintf(cfu,"   %d", cbobb[i][j][0]);
				//printf("It is being done!!\n");
			}
			fprintf(cfu, "\n");
			
			
			rpid=pid[i];
			
			fprintf(cfu,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",  li[i][0], li[i][1], li[i][2],  li[i][3], li[i][4], li[i][5],  li[i][6], li[i][7], li[i][8]);
		}
		
		
		fclose(cfu);
		
		
		
		
			
		
		
		//fprintf(dwu,"\n");
		
		
		
		
		
		
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
			
			fprintf(ku,"%lf   %lf    %lf\n", fv[0][0], fv[0][1], fv[0][2]);
			fprintf(ku,"%lf   %lf    %lf\n", fv[1][0], fv[1][1], fv[1][2]);
			fprintf(ku,"%lf   %lf    %lf\n", fv[2][0], fv[2][1], fv[2][2]);
			
			
			for (i=0; i<N; i++)
			{
				fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
				
				fprintf(ku, "%d\n", wnrd[i]);
				for (j=0; j<wnrd[i]; j++) 
				{
					fprintf(ku,"   %d    %d", nrdi[i][j][0], nrdi[i][j][1]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku, "%d\n", wcbobb[i]);
				for (j=0; j<wcbobb[i]; j++) 
				{
					fprintf(ku,"   %d     %d     %d", cbobb[i][j][0],  cbobb[i][j][1],  cbobb[i][j][2]);
				}
				fprintf(ku, "\n");
				
				for (j=0; j<npatch[pid[i]]; j++)
				{
					fprintf(ku, "%lf   %lf   %lf \n", paxe[i][j][0],   paxe[i][j][1] , paxe[i][j][2] );
				
				}
			
				
				fprintf(ku,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \n", li[i][0], li[i][1], li[i][2], li[i][3], li[i][4], li[i][5], li[i][6], li[i][7], li[i][8]);
				
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
			  		fprintf(fu,"%lf   ",lfp[jt][i]);
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
			
			
			
			
			
			
			for (it=0; it<tpatch; it++)
			{
				
				for (jt=it; jt<tpatch; jt++)
				{
					if (pop[it][jt]==1)
					{
						char filenameb[25] = {0};

		    				sprintf(filenameb, "Cord_bond_type_%d_%d.csv", it, jt);
						
						
						FILE* yu;
				     		yu=fopen(filenameb,"w");
				     	
						fprintf(yu,"xcord, ycord, zcord, icord, jcord, kcord \n");
						
						for (i=0; i<N; i++)
						{
							rpid1=pid[i];
							for (j=0; j<wcbobb[i]; j++)
							{
								rpid2=pid[cbobb[i][j][0]];
								v1=ami[i][0]-ami[cbobb[i][j][0]][0];
								v2=ami[i][1]-ami[cbobb[i][j][0]][1];
								v3=ami[i][2]-ami[cbobb[i][j][0]][2];
								
								if (v1<-f[0]/2.0) v1=v1+f[0];
								else if(v1>f[0]/2.0) v1=v1-f[0];
								if (v2<-f[1]/2.0) v2=v2+f[1];
								else if(v2>f[1]/2.0) v2=v2-f[1];
								if (v3<-f[2]/2.0) v3=v3+f[2];
								else if(v3>f[2]/2.0) v3=v3-f[2];
								
								
								re1 = v1*fv[0][0] + v2*fv[1][0] + v3*fv[2][0];
								re2 = v1*fv[0][1] + v2*fv[1][1] + v3*fv[2][1];
								re3 = v1*fv[0][2] + v2*fv[1][2] + v3*fv[2][2];
								v1=re1;
								v2=re2;
								v3=re3;
								
								ec=sqrt(v1*v1+v2*v2+v3*v3);
								v1=v1/ec;
								v2=v2/ec;
								v3=v3/ec;
								
								dxi=lfp[cbobb[i][j][0]][0]+v1*ec*0.5;
								dyi=lfp[cbobb[i][j][0]][1]+v2*ec*0.5;
								dzi=lfp[cbobb[i][j][0]][2]+v3*ec*0.5;
								
								if(patches[rpid1][cbobb[i][j][1]]==it || patches[rpid2][cbobb[i][j][2]]==it)
								{
									if (patches[rpid2][cbobb[i][j][2]]==jt || patches[rpid1][cbobb[i][j][1]]==jt)
									{
										if(i>cbobb[i][j][0])	fprintf(yu, "%lf, %lf, %lf, %lf, %lf, %lf \n", dxi, dyi, dzi, v1, v2, v3);
									}
								}
							}
							
						}
						fclose(yu);
					}
				}
				
				
			}
			
			
			FILE* yuu;
	     		yuu=fopen("Cord_bonds.csv","w");
	     	
			fprintf(yuu,"xcord, ycord, zcord, icord, jcord, kcord \n");
			
			for (i=0; i<N; i++)
			{
				
				for (j=0; j<wcbobb[i]; j++)
				{
					
					v1=ami[i][0]-ami[cbobb[i][j][0]][0];
					v2=ami[i][1]-ami[cbobb[i][j][0]][1];
					v3=ami[i][2]-ami[cbobb[i][j][0]][2];
					
					if (v1<-f[0]/2.0) v1=v1+f[0];
					else if(v1>f[0]/2.0) v1=v1-f[0];
					if (v2<-f[1]/2.0) v2=v2+f[1];
					else if(v2>f[1]/2.0) v2=v2-f[1];
					if (v3<-f[2]/2.0) v3=v3+f[2];
					else if(v3>f[2]/2.0) v3=v3-f[2];
					
					re1 = v1*fv[0][0] + v2*fv[1][0] + v3*fv[2][0];
					re2 = v1*fv[0][1] + v2*fv[1][1] + v3*fv[2][1];
					re3 = v1*fv[0][2] + v2*fv[1][2] + v3*fv[2][2];
					
					v1=re1;
					v2=re2;
					v3=re3;
					
					ec=sqrt(v1*v1+v2*v2+v3*v3);
					v1=v1/ec;
					v2=v2/ec;
					v3=v3/ec;
					
					dxi=lfp[cbobb[i][j][0]][0]+v1*ec*0.5;
					dyi=lfp[cbobb[i][j][0]][1]+v2*ec*0.5;
					dzi=lfp[cbobb[i][j][0]][2]+v3*ec*0.5;
					
					
					
					if(i>cbobb[i][j][0])	fprintf(yuu, "%lf, %lf, %lf, %lf, %lf, %lf \n", dxi, dyi, dzi, v1, v2, v3);
				}
				
			}
			fclose(yuu);
			
			
			
			
			
			
			
			
			
			
			for (i=0; i<tpatch; i++)
			{
				char filenamev[25] = {0};

    				sprintf(filenamev, "Cord_patch_vector_%d.csv", i);

				FILE* vfu;
				vfu=fopen(filenamev,"w");
				
				fprintf(vfu,"xcord, ycord, zcord, icord, jcord, kcord \n");
				
				for(j=0; j<N; j++)
				{
					jt=pid[j];
					for (k=0; k<npatch[jt]; k++)
					{
						if(patches[jt][k]==i)
						{
							fprintf(vfu, "%lf, %lf, %lf, %lf, %lf, %lf \n", lfp[j][0], lfp[j][1], lfp[j][2], paxe[j][k][0], paxe[j][k][1], paxe[j][k][2]);
						}
					}
				}
				fclose(vfu);
				
			
			}
			
			FILE* svfu;
			svfu=fopen("Cord_symmetry_axes.csv", "w");
			fprintf(svfu,"xcord, ycord, zcord, icord, jcord, kcord \n");
			for (i=0; i<N; i++)
			{
				rpid=pid[i];
				v1=0.0;
				v2=0.0;
				v3=0.0;
				for (j=0; j<npatch[rpid]; j++)
				{
					v1=v1+paxe[i][j][0];
					v2=v2+paxe[i][j][1];
					v3=v3+paxe[i][j][2];
				}
				bdm=sqrt(v1*v1+v2*v2+v3*v3);
				v1=v1/bdm;
				v2=v2/bdm;
				v3=v3/bdm;
				
				fprintf(svfu, "%lf, %lf, %lf, %lf, %lf, %lf \n", lfp[i][0], lfp[i][1], lfp[i][2], v1, v2, v3);
			
			}
			
			fclose(svfu);
			
			
			

			//VTK Zone ends!
			
			
			
			
		}
		
		
		
		
		

	}
	
	
	
	
	
	
	
	
	
	




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
		free(nrdi[k]);
	
	}
		
	free(nrdi);
	free(wnrd);
	
	
		
	for (vtbox=0; vtbox<tbox; vtbox++)
	{
		free(dri[vtbox]);	
	}
	free(dri);
	free(wdr);
	

	//Now freedom is achieved!	


	
}






