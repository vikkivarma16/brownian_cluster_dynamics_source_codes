#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
const int E=20, UN=200, na=60;

long int tsim, ci, dpuc, dpuc1,  vel, sjet, pjet, tsjet, tpjet;

//void ecfcal(int,double[]);

int  N, trim, check, checkl, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, nf, buf_block[3], buf_blockd, block_en[2][3], blockd_en[2], jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct,  veec, SHA, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt, gv, thi, nc, sdri, pboo, inf, vc, trims, ptt, kn, isitt, bvs, bvsi, ren_id, ren_id1, ren_id2, n_en1, n_en2, flag_v, flag_n, flag_d, flag_npi, flag_npit, count, n_check;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, cp, cph, cpsh, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9],  ar, per,  buf_f[3],  buf_bl[3], f_en[2][3], bl_en[2][3], ncv1, ncv2, ennch, g_area; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda,lembdai, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, buf_fe, fe_en[2], tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON,  tbl, pbc,  dif, difu, difus, difl, difusl, adefus, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm, dsss, pree, pre, vol, dlr, press, buf_fer, fer_en[2], ene, bf, en1_aro, en1_arn, en2_aro, en2_arn, fr_en[2][3], vstep, vpress, vpstep, divea, divei, dive, rsphere_en[2], cm_en[2][3];  



int main(void)
{

	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	fscanf(bu,"%lf",&bdm);
	N=2*(int)round(bdm);
		
     	srand48(time(NULL));
	srand(time(NULL));
	
	//if (N==0) printf("raam tere desh men \n");
	// Definition of universal constants starts!
	
	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;
	pboo=1;
	
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
	ivf=0.17;
	
	// Definition of global constants ends       !!!!!!!!!!!!!!!!
	
	
  	// Definition of simulation constants starts        !!!!!!!!!!!!!!!!
	
	
	
	fscanf(bu,"%lf",&bdm);
	fvf=bdm;
	
	
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
	
	
	
	
	// Definition of compression parameter ends            !!!!!!!!!!
	
	

	// Declaration of core bond variable starts          !!!!!!!!!!
	
	double epsi;
	
	// Declaration of core bond variable ends         !!!!!!!!!!
	
	
	
	// Evaluation of bond variable starts           !!!!!!!!!
	
	
	
	bn=40;
	bnll=40;
	sst=1;
	
	
	//double epsi1, epsi2, epsi3, omega2, beta1, beta2, beta3, betai;
	
	// Declaration of core bond variable ends!
	
	
	// Evaluation of bond variable starts!
	

	//pop =  0 ;              // 0: hard core;   1: HC+Jenus;    2: HC+Patchy;    3: HC+Isotropic;     4: HC+jenus+Isotropic;   5: HC+patchy+Isotropic; 

	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			pop[i][j]=0;
			pop[j][i]=0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
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



	
	
	// Potential properties define here!
	
	
	
	
	
	
	
	
	
	double* poten=malloc((((war*(war-1))/2)+war)*sizeof(double*));
	
	// Jenus potential properties;
	
	double** epsi1=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi1[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta1=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta1[i]=malloc(war*sizeof(double)); 
	
	}
	
	
	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi1[i][j]=0.0;
			epsi1[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi1[i][j]=bdm;
				epsi1[j][i]=bdm;
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
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta1[i][j]=0.0;
			beta1[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta1[i][j]=bdm;
				beta1[j][i]=bdm;
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
	
	double** epsi2=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi2[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta2=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta2[i]=malloc(war*sizeof(double)); 
	
	}
	
	double* omega2=malloc(war*sizeof(double));
	double* wom=malloc(war*sizeof(double));
	double* mom=malloc(war*sizeof(double));
	double* maom=malloc(war*sizeof(double));
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi2[i][j]=0.0;
			epsi2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
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
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta2[i][j]=0.0;
			beta2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
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
	
	
	
	
	for (i=0; i<war; i++)
	{
		omega2[i]=0.0;
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
			omega2[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	for (i=0; i<war; i++)
	{
		wom[i]=0.0;
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
			wom[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	for (i=0; i<war; i++)
	{
		mom[i]=omega2[i]-wom[i];
	
	}
	for (i=0; i<war; i++)
	{
		maom[i]=omega2[i]+wom[i];
	
	}
	
	
	
	// Isotropic potential properties;
	
	double** epsi3=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi3[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta3=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta3[i]=malloc(war*sizeof(double)); 
	
	}
	

	
	
	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi3[i][j]=0.0;
			epsi3[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi3[i][j]=bdm;
				epsi3[j][i]=bdm;
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
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta3[i][j]=0.0;
			beta3[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta3[i][j]=bdm;
				beta3[j][i]=bdm;
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
	
	
	
	
	
			
		
	
	
	// Thermodynamic protocol active region ends !!!!!!!!!!!!!!!!!!!!
	

	fclose(bu);
	remove("Base_data.txt");



	
	
	// Declaired!
	
	double** beta=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** betai=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		betai[i]=malloc(war*sizeof(double)); 
	
	}
	
	
	double** ur=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		ur[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** uri=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		uri[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** urh=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		urh[i]=malloc(war*sizeof(double)); 
	
	}
	
	
	epsi=0.0;
	for(i=0; i<war; i++)
	{
		for(j=0; j<war; j++)
		{
			
	
			if (pop[i][j]==0)
			{
				epsi=0.0;
			}
			
			else if (pop[i][j]==1)
			{
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi1[i][j]>epsi) epsi=epsi1[i][j];
						beta[i][j]=beta1[i][j];
						beta[j][i]=beta1[i][j];
						
					}
				}
			}
			else if (pop[i][j]==2)
			{
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi2[i][j]>epsi) epsi=epsi2[i][j];
						beta[i][j]=beta2[i][j];
						beta[j][i]=beta2[i][j];
						
					}
				}
			}
			else if (pop[i][j]==3)
			{
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
				
				
			}
			else if (pop[i][j]==4)
			{
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi1[i][j]>epsi) epsi=epsi1[i][j];
						beta[i][j]=beta1[i][j];
						beta[j][i]=beta1[i][j];
						
					}
				}
				
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
				
				
			}	
			else if (pop[i][j]==5)
			{
				
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi2[i][j]>epsi) epsi=epsi2[i][j];
						beta[i][j]=beta2[i][j];
						beta[j][i]=beta2[i][j];
						
					}
				}
				
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
			}
			
		}
	}
	
	for(i=0; i<war; i++)
	{
		for(j=0; j<war; j++)
		{
			ur[i][j]=log(1.0-(1.0/(1.0+beta[i][j]))+EPSILON);
			uri[i][j]=log(1.0-(1.0/(1.0+betai[i][j]))+EPSILON);
			urh[i][j]=1.0;
		}	
	}
	
	printf("________________________________________fffffffffffffffffff___________________________%lf  %d  \n\n\n\n\n", ur[0][0], pop[0][0]);
	fvf=ivf;
	
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
	
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			 printf("Kind of potential between family of particles %d and %d is defined as-   ", i,j);
			 if (pop[i][j]==0)
			 {
			 	printf("Only hard core.\n");
			 	printf("Epsi-   Not applicable.\n");
			 	printf("Beta-   Not applicable.\n");
			 	printf("Omega-   Not applicable.\n");
			 	
			 }
			 else if (pop[i][j]==1)
			 {
			 	printf("Hard core + Jenus.\n");
			 	printf("Epsi-   %lf \n", epsi1[i][j]);
			 	printf("Beta-   %lf \n", beta1[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	 
			 
			  else if (pop[i][j]==2)
			 {
			 	printf("Hard core + Patchy.\n");
			 	printf("Epsi-   %lf \n", epsi2[i][j]);
			 	printf("Beta-   %lf \n", beta2[i][j]);
			 	printf("Omega for particle %d, for the patchy potential defined between type of particle %d and %d-   %lf \n", i, i, j, omega2[i]);
			 }	
			 
			  else if (pop[i][j]==3)
			 {
			 	printf("Hard core + Jenus.\n");
			 	printf("Epsi-   %lf \n", epsi3[i][j]);
			 	printf("Beta-   %lf \n", beta3[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	
			  else if (pop[i][j]==4)
			 {
			 	printf("Hard core + Jenus + Isotropic.\n");
			 	printf("Epsi and Epsi-isotropic-   %lf  %lf \n", epsi1[i][j], epsi3[i][j]);
			 	printf("Beta and Beta-isotropic-   %lf  %lf \n", beta1[i][j], beta3[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	
			 
			 else if (pop[i][j]==5)
			 {
			 	printf("Hard core + Patchy + Isotropic.\n");
			 	printf("Epsi and Epsi-isotropic-   %lf  %lf \n", epsi2[i][j], epsi3[i][j]);
			 	printf("Beta-\n");
			 	printf("Beta and Beta-isotropic-   %lf  %lf \n", beta2[i][j], beta3[i][j]);
			 	printf("Omega for particle %d, for the patchy potential defined between type of particle %d and %d-   %lf \n", i, i, j, omega2[i]);
			 }	
			 printf("\n");
			 printf("\n");
			 printf("\n");
		 
		}
	
	}
	
	printf("Selected Epsi to harnesh the particle in the range %lf \n", epsi);
	//printf("Values corresponding to the particular parameter is given by\n");
	
	
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	
	
	dlr=0.001; // It is defined for the pressure calculation width, where we have to take all the pairs lying within the range !!!!!!!!!!!!!!!! 
	
	
	
	
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
			
			e[i][0]=dlr+rm[i]+(epsi/2.0);
			e[i][1]=dlr+rm[i]+(epsi/2.0);
			e[i][2]=dlr+r[i]+(epsi/2.0);
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
				
				e[i][0]=dlr+rm[i]+(epsi/2.0);
				e[i][1]=dlr+rm[i]+(epsi/2.0);
				e[i][2]=dlr+r[i]+(epsi/2.0);
			}
			else
			{
				d[i]=sar[i]*ar[i];
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=dlr+rm[i]+(epsi/2.0);
				e[i][1]=dlr+rm[i]+(epsi/2.0);
				e[i][2]=dlr+r[i]+(epsi/2.0);
			
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
	//printf("%lf\n", ds[1][0]);
	buf_fe=0.0;
	ivf=fvf;
	
	for (i=0; i<war; i++)
	{
		buf_fe=buf_fe+((piee*r[i]*r[i]*(float)nt[i])/(ar[i]*ar[i]*fvf));
	}
	buf_fe=buf_fe/(4.0*piee);
	rsphere_en[0]=pow(buf_fe,(1.0/2.0));
	rsphere_en[1]=rsphere_en[0];
	buf_fe=4*rsphere_en[0];
	buf_fer=buf_fe;
	buf_f[0]=buf_fe;
	buf_f[1]=buf_f[0];
	buf_f[2]=buf_f[0];
	
	buf_fer=buf_f[0];
	buf_fe=buf_f[0];
	
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
	
	ble=ble+epsi+dlr;
	
	
	
	
	
	
	buf_bl[0]=ble;
	buf_bl[1]=ble;
	buf_bl[2]=ble;
	
	buf_block[0]=(int)floor(buf_f[0]/buf_bl[0]);
	buf_bl[0]=buf_f[0]/(float)buf_block[0];
	buf_f[0]=buf_bl[0]*(float)buf_block[0];
	
	buf_block[1]=(int)floor(buf_f[1]/buf_bl[1]);
	buf_bl[1]=buf_f[1]/(float)buf_block[1];
	buf_f[1]=buf_bl[1]*(float)buf_block[1];
	
	buf_block[2]=(int)floor(buf_f[2]/buf_bl[2]);
	buf_bl[2]=buf_f[2]/(float)buf_block[2];
	buf_f[2]=buf_bl[2]*(float)buf_block[2];
	
	buf_blockd=buf_block[0]+20;
	
	
	
	
	
	
	
	
	
	
	
	
	printf("                                ******:::::::   Calculated parameters are given by   :::::::******\n");
	printf("\n");
	printf("Starting edge is given by-  %lf\n", buf_f[0]);
	printf("\n");
	printf("\n");
	printf("Aimed edge is given by-  %lf\n", buf_fe);
	printf("\n");
	printf("\n");
	printf("Radius of sphere -  %lf\n", rsphere_en[0]);
	printf("\n");
	printf("\n");
	printf("Starting number of shell at each edge is given by-  %d\n", buf_block[0]);
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
	
	
	
	
	double** ami=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		ami[i]=malloc(3*sizeof(double));
	}
	int* pid=malloc(N*sizeof(int));
	
	
	
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
	
	
	
	int* ensemble1=malloc(N*sizeof(int));
	int* ensemble2=malloc(N*sizeof(int));
	
	
	
	double ldcxyz[N][3],  lom[N][3];
	
	
	
	int** nrdi=malloc(N*sizeof(int*));
	int* wnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		nrdi[i]=malloc(UN*sizeof(int));
	}
	
	int** tnrdi=malloc(N*sizeof(int*));
	int* wtnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		tnrdi[i]=malloc(UN*sizeof(int));
	}
	
	
	int**** dri=malloc(buf_blockd*sizeof(int***));
	int*** wdr=malloc(buf_blockd*sizeof(int**));
		
	for (j=0; j<buf_blockd; j++)
	{
		dri[j]=malloc(buf_blockd*sizeof(int**));
		wdr[j]=malloc(buf_blockd*sizeof(int*));
		for (i=0; i<buf_blockd; i++)
		{
			dri[j][i]=malloc(buf_blockd*sizeof(int*));
			wdr[j][i]=malloc(buf_blockd*sizeof(int));
	
			for (k=0; k<buf_blockd; k++)
			{
				dri[j][i][k]=malloc(30*sizeof(int));
			}
		}	
	}
	
	
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

	int** bobbi=malloc(N*sizeof(int*));
	int* wbobbi=malloc(N*sizeof(int));
	
	
	int** cbobbi=malloc(N*sizeof(int*));
	int* wcbobbi=malloc(N*sizeof(int));
	for (k=0; k<N; k++)
	{
		bobbi[k]=malloc(E*sizeof(int));
		cbobbi[k]=malloc(E*sizeof(int));
	
	}	
		
	
	
	int wpbobbi; 
	int psbobbi[E];

	int  tnkx[27][3];
		
		
		
		
	int* ensemble_id=malloc(N*sizeof(int));
	
	jt=0;
	
	for (i=0; i<N; i++)
	{
		ensemble1[i]=6000;
		ensemble2[i]=6000;
	
	}
	
	n_en1=0;
	n_en2=0;
	for(i=0; i<war; i++)
	{
		jdt=nt[i];
		for(j=jt; j<jt+jdt; j++)
		{
			if(j<0.5*(jt+jdt))
			{
				ensemble_id[j]=0;
				ensemble1[n_en1]=j;
				n_en1=n_en1+1;
			}
			else 
			{
				ensemble_id[j]=1;
				ensemble2[n_en2]=j;
				n_en2=n_en2+1;
			}
		}
		jt=jt+jdt;
	}
	
	
	for(i=0; i<N; i++)
	{
		//printf("Ensemble identity:  %d  %d\n", i,  ensemble_id[i]);
	}
	printf("%d  %d  \n", n_en1, n_en2);

	
	// Declaration of bonding variable ends! 
	
	
	
	
	
	
	
	// Initialization of core variables start!
	
	
	adefus=0.0;
	adefusl=0.0;
	for (i=0; i<N; i++)
	{
		
		
		
		ldcxyz[i][0]=0.0;
		ldcxyz[i][1]=0.0;
		ldcxyz[i][2]=0.0;
		
		lom[i][0]=0.0;
		lom[i][1]=0.0;
		lom[i][2]=0.0;
		
		wtnrd[i]=0;
		wnrd[i]=0;
		for (j=0; j<UN; j++)
		{
			nrdi[i][j]=6000;
			tnrdi[i][j]=6000;
		
		}	
		
	}
	
	for (i=0; i<buf_blockd; i++)
	{
		
		for (j=0; j<buf_blockd; j++)
		{
			for (k=0; k<buf_blockd; k++)
			{
				wdr[i][j][k]=0;
				for (l=0; l<30; l++)
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
	
	
	
	// Initialization of bond variables starts!
	
	for (i=0; i<N; i++)
	{
		wbobb[i]=0;
		wcbobb[i]=0;
		for (j=0; j<E; j++)
		{
			bobb[i][j]=6000;
			bobbi[i][j]=6000;
			cbobb[i][j]=6000;
			cbobbi[i][j]=6000;
			
		}
		
		wbobbi[i]=0;
		wcbobbi[i]=0;
						
	}
		
	
	wpbobbi=0;
	wpbobb=0;
	for(i=0; i<E; i++)
	{
		psbobb[i]=6000;
		psbobbi[i]=6000;
	}
	
	// Initialization of the bond variable ends here !!!!!!!!!!!!!!!!
	
	
	
	
	
	
	
	
	
	
	for (i=0; i<N; i++)
	{
		jt=0;
		jtt=0;
		for (j=0; j<war; j++)
		{
			jtt=jtt+(int)round(per[j]*(float)N);
			if (i>jtt)    jt=jt+1;		
		}
		pid[i]=jt;
		//printf("%d   %d\n", i , pid[i]);
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
	
	

	
	
	
	
	
	
	// Setting initial coordinates and direction starts!
	
	
	
	FILE* cfu;
	
	FILE* vu;
	vu=fopen("posB-p0.4002-s0.0100-L020-t0.000000.dat","r");
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	if (vu==NULL)
	{
		printf("Initial configuration within the simulation box have been started generating: \n \n");
		
		for ( i=0; i<3; i++)
		{
			f_en[0][i]=buf_fe;
			block_en[0][i]=(int)(floor(f_en[0][i]/ble));
			bl_en[0][i]=f_en[0][i]/(float)block_en[0][i];
			f_en[0][i]=bl_en[0][i]*(float)block_en[0][i];
			//printf("ptr %d \n", ptr);
		}
		
		
		for ( i=0; i<3; i++)
		{
			f_en[1][i]=buf_fe;
			block_en[1][i]=(int)(floor(f_en[1][i]/ble));
			bl_en[1][i]=f_en[1][i]/(float)block_en[1][i];
			f_en[1][i]=bl_en[1][i]*(float)block_en[1][i];
		}
		
		fe_en[0]=buf_fe;
		fer_en[0]=buf_fer;
		
		fe_en[1]=buf_fe;
		fer_en[1]=buf_fer;
		
		cm_en[0][0]=0.5*fe_en[0];
		cm_en[0][1]=0.5*fe_en[0];
		cm_en[0][2]=0.5*fe_en[0];
		cm_en[1][0]=0.5*fe_en[1];
		cm_en[1][1]=0.5*fe_en[1];
		cm_en[1][2]=0.5*fe_en[1];
		
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
			ren_id1=ensemble_id[i];
			
			ra1=drand48();
	        	theta=acos(1.0-(2.0*ra1));
			ra2=drand48();
	       		phi=tpiee*ra2;
			
			pdx=sin(theta)*cos(phi);
			pdy=sin(theta)*sin(phi);
			pdz=cos(theta);
			
			xi=rsphere_en[ren_id1]*pdx+cm_en[ren_id1][0];
	    		yi=rsphere_en[ren_id1]*pdy+cm_en[ren_id1][1];
	    		zi=rsphere_en[ren_id1]*pdz+cm_en[ren_id1][2];
			
			xk=(int)floor(xi/bl_en[ren_id1][0]);
	    		yk=(int)floor(yi/bl_en[ren_id1][1]);
	    		zk=(int)floor(zi/bl_en[ren_id1][2]);
			
	    		dum=wdr[xk][yk][zk];
			dri[xk][yk][zk][dum]=i;
			wdr[xk][yk][zk]=dum+1;
			ami[i][0]=xi;
			ami[i][1]=yi;
			ami[i][2]=zi;
							
		}
		
		cow=0;
		
		//cow=0;
		sitts=0;
		quao=0.0;
		quaoi=0.0;
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		printf("Initial configuration generated successfully!!!!!!! \n \n");
		fclose(cfu);
	}
	else
	{
		printf("A pre-prepared system found and intial configuration will be set accordingly. \n");
		printf("Initial configuration is being read from the available file!!! \n \n  ");
		
		cow=1;
		
		fscanf(vu, "%d", &sitts);
		fscanf(vu, "%lf", &quao);
		fscanf(vu, "%lf", &quaoi);
			
		fscanf(vu,"%lf", &f_en[0][0]);
		fscanf(vu,"%lf", &f_en[0][1]);
		fscanf(vu,"%lf", &f_en[0][2]);
		
		
		fscanf(vu,"%lf", &f_en[1][0]);
		fscanf(vu,"%lf", &f_en[1][1]);
		fscanf(vu,"%lf", &f_en[1][2]);
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		
		fclose(cfu);
		
		for (i=0; i<N; i++)
		{
			
			for(l=0;l<3;l++)
			{
				fscanf(vu,"%lf", &ami[i][l]);
			}
			
			fscanf(vu,"%d", &wnrd[i]);
			for (j=0; j<wnrd[i]; j++) 
			{
				fscanf(vu,"%d", &nrdi[i][j]);
			}
			
			fscanf(vu,"%d", &wbobb[i]);
			for (j=0; j<wbobb[i]; j++) 
			{
				fscanf(vu,"%d", &bobb[i][j]);
			}
			
			fscanf(vu,"%d", &wbobbi[i]);
			for (j=0; j<wbobbi[i]; j++) 
			{
				fscanf(vu,"%d", &bobbi[i][j]);
			}
			
			for(l=6;l<9;l++)
			{
				fscanf(vu,"%lf", &li[i][l]);
				
			}
			
		}
		
		
		
		for ( i=0; i<3; i++)
		{
			
			block_en[0][i]=(int)(floor(f_en[0][i]/ble));
			bl_en[0][i]=f_en[0][i]/(float)block_en[0][i];
			f_en[0][i]=bl_en[0][i]*(float)block_en[0][i];
			//printf("ptr %d \n", ptr);
		}
		
		
		for ( i=0; i<3; i++)
		{
			
			block_en[1][i]=(int)(floor(f_en[1][i]/ble));
			bl_en[1][i]=f_en[1][i]/(float)block_en[1][i];
			f_en[1][i]=bl_en[1][i]*(float)block_en[1][i];
			//printf("ptr %d \n", ptr);
		}
		
		
		for ( i=0; i<3; i++)
		{
			
			buf_block[i]=(int)(floor(f_en[0][i]/ble));
			buf_bl[i]=f_en[0][i]/(float)buf_block[i];
			buf_f[i]=buf_bl[i]*(float)buf_block[i];
			//printf("ptr %d \n", ptr);
		}
		
		for(i=0; i<N; i++)	
	    	{	
	    		
	    		ren_id1=ensemble_id[i];
	    		
	    		xk=(int)floor(ami[i][0]/bl_en[ren_id1][0]);
	    		yk=(int)floor(ami[i][1]/bl_en[ren_id1][1]);
	    		zk=(int)floor(ami[i][2]/bl_en[ren_id1][2]);
	    	
			dum=wdr[xk][yk][zk];
			dri[xk][yk][zk][dum]=i;
			wdr[xk][yk][zk]=dum+1;
			
		}
		cow=1;
		
		printf("Initial configuration have been copied from the file and implemented into the system successfully!!!!!!! \n \n  ");
		
			//printf("Go and yourself\n");
		
		fclose(vu);
	}
	
	//printf("raam\n");
	
	en1_arn=4*piee*rsphere_en[0]*rsphere_en[0];
	en2_arn=4*piee*rsphere_en[1]*rsphere_en[1];
	
	g_area=en1_arn+en2_arn;
	
	
		
     	
     	// Thermodynamic protocol active region starts  !!!!!!!
     	
     	int*** pres=malloc(N*sizeof(int**));
     	for (i=0; i<N; i++)
     	{
     		pres[i]=malloc(4*sizeof(int*)); /// as there are four kind of potential     !!!!!!!!!!!!!!
     		for(j=0; j<4; j++)
     		{
     			pres[i][j]=malloc(E*sizeof(int));
     		}
     		
     	}
     	
     	
     	
     	int** wpres=malloc(N*sizeof(int*));
     	for (i=0; i<N; i++)
     	{
     		wpres[i]=malloc(4*sizeof(int));
     	}
     	
     	for (i=0; i<N; i++)
     	{
     		for(k=0; k<4; k++)
     		{
     		
	     		for (j=0; j<E; j++)
	     		{
	     			pres[i][k][j]=6000;
	     		}
	     		wpres[i][k]=0; 
	     	
	     	}
     	
     	}	
     	 
     	int ppres[4][E], wppres[4]; 
     	 


     	
    
     	
     	
	FILE* pu;
	pu=fopen("Pressure_histrogram.txt","w");
	
	fclose(pu);
	
	
	
     	printf("%lf\n dssssss", dss[0][0]);
     	
	
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
	
	
	
	
	
	
	
	FILE * hu;
	hu=fopen("Initial_conf.vtk","w");

	fprintf(hu,"# vtk DataFile Version 3.0\n");
	fprintf(hu,"Random data to test tensors\n");
	fprintf(hu,"ASCII\n");
	fprintf(hu,"DATASET POLYDATA\n");
	fprintf(hu,"POINTS %d float\n", N);
	
	for (jt=0; jt<N; jt++)
	{
		for(i=0;i<3;i++)
		{
	  		fprintf(hu,"%lf   ",ami[jt][i]);
		}
		fprintf(hu,"\n");
	}

	fprintf(hu,"\n");
	fprintf(hu,"POINT_DATA %d\n", N);
	fprintf(hu,"\n");
	fprintf(hu,"\n");
	fprintf(hu,"TENSORS non-spherical_ellipsoid float\n");
	fprintf(hu,"\n");
	
	for (i=0; i<N; i++)
	{
		rpid1=pid[i];
		del[0]=(rm[rpid1]*li[i][0]*li[i][0])+(rm[rpid1]*li[i][3]*li[i][3])+(r[rpid1]*li[i][6]*li[i][6]);
		del[1]=(rm[rpid1]*li[i][0]*li[i][1])+(rm[rpid1]*li[i][3]*li[i][4])+(r[rpid1]*li[i][6]*li[i][7]);
		del[2]=(rm[rpid1]*li[i][0]*li[i][2])+(rm[rpid1]*li[i][3]*li[i][5])+(r[rpid1]*li[i][6]*li[i][8]);
		del[3]=(rm[rpid1]*li[i][0]*li[i][1])+(rm[rpid1]*li[i][3]*li[i][4])+(r[rpid1]*li[i][6]*li[i][7]);
		del[4]=(rm[rpid1]*li[i][1]*li[i][1])+(rm[rpid1]*li[i][4]*li[i][4])+(r[rpid1]*li[i][7]*li[i][7]);
		del[5]=(rm[rpid1]*li[i][1]*li[i][2])+(rm[rpid1]*li[i][4]*li[i][5])+(r[rpid1]*li[i][7]*li[i][8]);
		del[6]=(rm[rpid1]*li[i][0]*li[i][2])+(rm[rpid1]*li[i][3]*li[i][5])+(r[rpid1]*li[i][6]*li[i][8]);
		del[7]=(rm[rpid1]*li[i][1]*li[i][2])+(rm[rpid1]*li[i][4]*li[i][5])+(r[rpid1]*li[i][7]*li[i][8]);
		del[8]=(rm[rpid1]*li[i][2]*li[i][2])+(rm[rpid1]*li[i][5]*li[i][5])+(r[rpid1]*li[i][8]*li[i][8]);
		
		fprintf(hu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
		fprintf(hu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
		fprintf(hu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
		fprintf(hu,"\n");
	}
	fclose(hu);

	
	
	
	
	
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
     
     
     	
     	FILE* mmwwu;
     	mmwwu=fopen("Diffusion_coefficient_el_l.txt","w");
     	fclose (mmwwu);
     	
     
     	
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
     
     
     
     
     
     
     
     	//press=-ur[0][0]*press;
	
	
	
	printf("Calculation is being done for the pressure:    %lf  \n\n\n", press);
	
     
     	
     	
     	
     	
     	
     	
	lsd=(int)floor((1.0/(step*step)));
	iltt=(int)floor((float)tsim/(float)lsd);
	rdd=(float)lsd/10.0;
	
	
	kd=0.0;
	istep=1.0/(step*step);
	vpstep=0.05;
	vstep=0.02;

	titi=0;
	ri=0;
	nlies=2*N;
	
	isitt=1;
	
	
	
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
	
	
	printf("Volume change Start !!!!!  %lf   %lf  %d\n", fer_en[0], fe_en[0], block_en[0][0]);
	
	flag_d=nlies;
	//exit(0);
	nlies=nlies+2;
	
	flag_npit=(int)(0.01*N*0.5);
	
	
	
	cow=0;
	
	
	int* flag=malloc(N*sizeof(int));
	
	int** tpop=malloc(war*sizeof(int*));
	for (i=0; i<war; i++)
	{
		tpop[i]=malloc(war*sizeof(int));
	}
	
	for (i=0; i<war; i++)
	{
		for(j=0; j<war; j++)
		{
			tpop[i][j]=pop[i][j];
		
		}
	
	}
	if ( cow==0)
	{
		
		for (i=0; i<war; i++)
		{
			for(j=0; j<war; j++)
			{
				pop[i][j]=0;
			
			}
		
		}
		
		
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
	double* displacer=malloc(3*sizeof(double));
	int wdisplacer;
	double rotator;
	
	
	printf("Here is the radius 1: %lf and 2: %lf \n", rsphere_en[0], rsphere_en[1]);
	
	printf("Here is the center of mass 1: (%lf   %lf   %lf ) and 2:  (%lf   %lf   %lf)\n", cm_en[0][0], cm_en[0][1], cm_en[0][2], cm_en[1][0], cm_en[1][1], cm_en[1][2]);
	
	printf("Here is the center of mass 2: (%lf   %lf   %lf ) and 2:  (%lf   %lf   %lf)\n", f_en[0][0], f_en[0][1], f_en[0][2], f_en[1][0], f_en[1][1], f_en[1][2]);
	
	
	for(sitt=sitts; sitt<iltt; sitt++)
	{	
		
		
		sjet=0;
		pjet=0;
		tsjet=0;
		tpjet=0;
		if (sitt==sitts+isitt)
		{
		
			for (i=0; i<N; i++)
			{
				ldcxyz[i][0]=0.0;
				ldcxyz[i][1]=0.0;
				ldcxyz[i][2]=0.0;
				
				lom[i][0]=0.0;
				lom[i][1]=0.0;
				lom[i][2]=0.0;
			}
			
			
			for (i=0; i<war; i++)
			{
				for(j=0; j<war; j++)
				{
					pop[i][j]=tpop[i][j];
				
				}
			
			}
			
			vpstep=0.05;
		}
		
		
		flag_n=2;
		
		for(ittt=0; ittt<lsd; ittt++)
		{	
	
			jet=0;
			vc=0;
			ptt=N;
			for (i=0; i<war; i++) trac[i]=0;
		     	nc=0;
		     	flag_v=0;
		     	
		     	flag_npi=0;
		     	
			if (cow==0)
			{
				ittt=0;
				cow=1;
				count=0;
				nlies=flag_d+2;
				for (i=0; i<N; i++)
				{
					if (flag[i]==1)		{cow=0; count=count+1;}
				}
				//printf("You got %d number of particles left behind the line.\n", count);
				
				
				if(cow==1)
				{
					
					FILE* ku;
     					ku=fopen("Compressed_system.txt","w");
     					fprintf(ku,"%d\n", sitt+1);
			
					fprintf(ku,"%lf\n", quao);
					
					fprintf(ku,"%lf\n", quaoi);
					
					
					fprintf(ku,"%lf   %lf   %lf \n", f_en[0][0], f_en[0][1], f_en[0][2]);
					fprintf(ku,"%lf   %lf   %lf \n", f_en[1][0], f_en[1][1], f_en[1][2]);
				
					for (i=0; i<N; i++)
					{
						fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
						
						fprintf(ku, "%d\n", wnrd[i]);
						for (j=0; j<wnrd[i]; j++) 
						{
							fprintf(ku,"   %d", nrdi[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku, "%d\n", wbobb[i]);
						for (j=0; j<wbobb[i]; j++) 
						{
							fprintf(ku,"   %d", bobb[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku, "%d\n", wbobbi[i]);
						for (j=0; j<wbobbi[i]; j++) 
						{
							fprintf(ku,"   %d", bobbi[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
					}
					
					
					
					fclose(ku);
					printf("Compression happend!!! Now the actual simulation will start!!! \n");
					
					FILE* fu;
				     	fu=fopen("Initial_cord_ensemble_1.vtk","w");
					fprintf(fu,"# vtk DataFile Version 3.0\n");
					fprintf(fu,"Random data to test tensors\n");
					fprintf(fu,"ASCII\n");
					fprintf(fu,"DATASET POLYDATA\n");
					fprintf(fu,"POINTS %d float\n", (int)(n_en1));
					
					
					FILE* fu1;
				     	fu1=fopen("Initial_cord_ensemble_2.vtk","w");
					fprintf(fu1,"# vtk DataFile Version 3.0\n");
					fprintf(fu1,"Random data to test tensors\n");
					fprintf(fu1,"ASCII\n");
					fprintf(fu1,"DATASET POLYDATA\n");
					fprintf(fu1,"POINTS %d float\n", (int)(n_en2));
					
					
					for (jt=0; jt<N; jt++)
					{
						if (ensemble_id[jt]==0)
						{
							
						  	fprintf(fu,"%lf    %lf    %lf  \n",ami[jt][0], ami[jt][1], ami[jt][2]);
						}
						else
						{
							fprintf(fu1,"%lf    %lf    %lf  \n",ami[jt][0]+f_en[0][0]+3.0, ami[jt][1], ami[jt][2]);	
						}
					}

					fprintf(fu,"\n");
					fprintf(fu,"POINT_DATA %d\n", (int)(n_en1));
					fprintf(fu,"\n");
					fprintf(fu,"\n");
					fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
					fprintf(fu,"\n");
					
					fprintf(fu1,"\n");
					fprintf(fu1,"POINT_DATA %d\n", (int)(n_en2));
					fprintf(fu1,"\n");
					fprintf(fu1,"\n");
					fprintf(fu1,"TENSORS non-spherical_ellipsoid float\n");
					fprintf(fu1,"\n");
					
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
						
						if (ensemble_id[i]==0)
						{
							fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
							fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
							fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
							fprintf(fu,"\n");
							
						}
						else
						{
							fprintf(fu1,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
							fprintf(fu1,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
							fprintf(fu1,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
							fprintf(fu1,"\n");
						
						
						}
						
					}
					
					
					fclose(fu1);
					fclose(fu);
					
				}
				else
				{
					nlies=nlies-2;
				}
				//printf("You got %d number of particles left behind the line.   %d\n", count, nlies);	 
			}
		
		
			//nlies=flag_d;
		
			n_check=0;
		
			for(mo=0; mo<nlies; mo++)
			{
            			
            			//printf("Buddha was born to decore the world!!!  %d\n", mo );
            			if (mo<flag_d)
            			{
		    			ntt=ptt;
		    			
		    			k=rand()  % ntt;
		    			
		    			boss=drand48();
		    			
		    			//if (cow==0 ) boss=0.4;
		    			rpid1=pid[k];
		    			
		    			
		    		}
		    		else if(mo==flag_d)
		    		{
		    		
		    		
		    			if (vc==0)
		    			{
		    				tpjet=tpjet+1;
		    				mo=mo-1;
		    				vc=1;
		    				k=0;
		    				rpid1=pid[k];
		    				boss=0.4;
		    				
		    				fer_en[0]=rsphere_en[0];
		    				fer_en[1]=rsphere_en[1];
		    				
		    				trims=2;
		    				//rd=drand48();
		    				//printf("Volume change Start !!!!!  %lf   %lf \n", fer_en[0], fe_en[0]);
		    				
		    				
		    				
		    				//printf("Volume change Start before !!!!!  %lf  %lf  \n", rsphere_en[0], rsphere_en[1]);
		    				
		    				en1_aro=4*piee*rsphere_en[0]*rsphere_en[0];
		    				en2_aro=4*piee*rsphere_en[1]*rsphere_en[1];
		    				bdm=log(en1_aro/en2_aro)+(drand48()-0.5)*vpstep;
		    				en1_arn=g_area*exp(bdm)/(1.0+exp(bdm));
		    				
		    				rsphere_en[0]=pow((en1_arn/(4.0*piee)),0.5);
		    				fe_en[0]=rsphere_en[0];
		    				fr_en[0][0]=fe_en[0]/fer_en[0];
		    				
			    			en2_arn=g_area-en1_arn;
			    			
			    			rsphere_en[1]=pow((en2_arn/(4.0*piee)),0.5);
		    				fe_en[1]=rsphere_en[1];
		    				fr_en[1][0]=fe_en[1]/fer_en[1];
		    				
		    				
			    			
			    			///printf("Volume change Start  after !!!!!  %lf    %lf   %lf    %lf \n", rsphere_en[0], rsphere_en[1], en1_arn, en1_aro);
			    			
			    			
		    				for(i=0; i<N; i++)	
					    	{	
					    		
					    		ren_id=ensemble_id[i];
					    		xk=(int)floor(ami[i][0]/bl_en[ren_id][0]);
					    		yk=(int)floor(ami[i][1]/bl_en[ren_id][1]);
					    		zk=(int)floor(ami[i][2]/bl_en[ren_id][2]);
					    		
							dum=wdr[xk][yk][zk];
							for (j=0; j<dum; j++)
							{
								dri[xk][yk][zk][j]=6000;
							
							}
							wdr[xk][yk][zk]=0;
							
							ami[i][0]=(ami[i][0]-cm_en[ren_id][0])*fr_en[ren_id][0];
		    					ami[i][1]=(ami[i][1]-cm_en[ren_id][1])*fr_en[ren_id][0];
		    					ami[i][2]=(ami[i][2]-cm_en[ren_id][2])*fr_en[ren_id][0];
		    					ami[i][0]=ami[i][0]+cm_en[ren_id][0];
			    				ami[i][1]=ami[i][1]+cm_en[ren_id][1];
			    				ami[i][2]=ami[i][2]+cm_en[ren_id][2];
		    					
		    					for (j=0; j<wbobb[i]; j++)	bobb[i][j]=6000;
							wbobb[i]=0;
							
							for (j=0; j<wbobbi[i]; j++)	bobbi[i][j]=6000;
							wbobbi[i]=0;
							
							for (j=0; j<wtnrd[i]; j++)	tnrdi[i][j]=6000;
							wtnrd[i]=0;
							
							for (j=0; j<wcbobb[i]; j++)	bobb[i][j]=cbobb[i][j];
							wbobb[i]=wcbobb[i];
							
							for (j=0; j<wcbobbi[i]; j++)	bobbi[i][j]=cbobbi[i][j];
							wbobbi[i]=wcbobbi[i];
							
							for (j=0; j<wnrd[i]; j++)	tnrdi[i][j]=nrdi[i][j];
							wtnrd[i]=wnrd[i];
		    					
							
						}
		    				
		    				
		    				
						
						///printf("Volume change Start rewsrefs !!!!!\n");
						
						for(i=0; i<N; i++)	
					    	{	
							//printf("%d\n", ren_id);
							ren_id=ensemble_id[i];
							
							///printf("%d\n", ren_id);
							
							
					    		xk=(int)floor(ami[i][0]/bl_en[ren_id][0]);
					    		yk=(int)floor(ami[i][1]/bl_en[ren_id][1]);
					    		zk=(int)floor(ami[i][2]/bl_en[ren_id][2]);
					    		//printf("sdhfdf   %d  %d   %d  %lf   %lf   \n", xk, yk ,zk, bl_en[ren_id][0], f_en[ren_id][0]);
							dum=wdr[xk][yk][zk];
							
							
							dri[xk][yk][zk][dum]=i;
							wdr[xk][yk][zk]=dum+1;
							//if (xk>block[0] || yk>block[1] || zk>block[2]) printf("YOU ARE doomed%d\n",i);
							
						}
					}
					else if (vc==1)
					{
						k=k+1;
			    			
			    		//	printf("Buddha was born to decore the world %d\n", k );
			    			mo=mo-1;
			    			boss=0.4;
						rpid1=pid[k];
			    			if (k==N && trims==2)
			    			{
			    				vc=2;
			    				ene=0.0;
			    				
							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wbobb[i]; j++)
								{
									rpid2=pid[bobb[i][j]];
									ene=ene+ur[rpid1][rpid2];
								
								}
								
							}
							
							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wbobbi[i]; j++)
								{
									rpid2=pid[bobbi[i][j]];
									ene=ene+uri[rpid1][rpid2];
								
								}
								
							}
							bdm=0.0;
							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wcbobb[i]; j++)
								{
									rpid2=pid[cbobb[i][j]];
									bdm=bdm+ur[rpid1][rpid2];
								
								}
								
							}

							for(i=0; i<N; i++)
							{
								rpid1=pid[i];
								for (j=0; j<wcbobbi[i]; j++)
								{
									rpid2=pid[cbobbi[i][j]];
									bdm=bdm+uri[rpid1][rpid2];
								
								}
								
							}
							ene=(bdm-ene)/2.0;
							bf=exp(-(ene-(float)(n_en1)*log(en1_arn/en1_aro)-(float)(n_en2)*log(en2_arn/en2_aro)));
							
							//printf("%lf  %d %d\n", ene, ittt, flag_n);
							
							if (drand48()>bf)
							{
			    					trims=0;	
							}
							else
							{
								pjet=pjet+1;	
								
								ntt=N;
					    			
					    			k=rand()  % ntt;
					    			
				    				rpid1=pid[k];
					    			
					    			boss=drand48();
							}
			    				
						}
						
						
						if(trims==0)
						{
							vc=2;
						
							fr_en[0][0]=1.0/fr_en[0][0];
				    			rsphere_en[0]=rsphere_en[0]*fr_en[0][0];
				    				
				    			fr_en[1][0]=1.0/fr_en[1][0];
				    			rsphere_en[1]=rsphere_en[1]*fr_en[1][0];
				    			
				    			//printf("Volume change failed!!!!!  %lf  %lf  \n", rsphere_en[0], rsphere_en[1]);		    		
				    			
				    			for(i=0; i<N; i++)	
						    	{	
						    		ren_id=ensemble_id[i];
						    		xk=(int)floor(ami[i][0]/bl_en[ren_id][0]);
						    		yk=(int)floor(ami[i][1]/bl_en[ren_id][1]);
						    		zk=(int)floor(ami[i][2]/bl_en[ren_id][2]);
								dum=wdr[xk][yk][zk];
								for (j=0; j<dum; j++) dri[xk][yk][zk][j]=6000;
								wdr[xk][yk][zk]=0;
								ami[i][0]=(ami[i][0]-cm_en[ren_id][0])*fr_en[ren_id][0];
			    					ami[i][1]=(ami[i][1]-cm_en[ren_id][1])*fr_en[ren_id][0];
			    					ami[i][2]=(ami[i][2]-cm_en[ren_id][2])*fr_en[ren_id][0];
			    					ami[i][0]=ami[i][0]+cm_en[ren_id][0];
			    					ami[i][1]=ami[i][1]+cm_en[ren_id][1];
			    					ami[i][2]=ami[i][2]+cm_en[ren_id][2];
				    				
							}
			    				
							
							
							
							for(i=0; i<N; i++)	
						    	{	
								
						    		
						    		ren_id=ensemble_id[i];
						    		xk=(int)floor(ami[i][0]/bl_en[ren_id][0]);
						    		yk=(int)floor(ami[i][1]/bl_en[ren_id][1]);
						    		zk=(int)floor(ami[i][2]/bl_en[ren_id][2]);
								dum=wdr[xk][yk][zk];
								dri[xk][yk][zk][dum]=i;
								wdr[xk][yk][zk]=dum+1;
								
								for (j=0; j<wnrd[i]; j++) nrdi[i][j]=6000;
								wnrd[i]=0;
								
								for (j=0; j<wcbobb[i]; j++)  cbobb[i][j]=6000;
								wcbobb[i]=0;
								
								for (j=0; j<wcbobbi[i]; j++)	cbobbi[i][j]=6000;
								wcbobbi[i]=0;
								
								for (j=0; j<wtnrd[i]; j++) nrdi[i][j]=tnrdi[i][j];
								wnrd[i]=wtnrd[i];
								
								for (j=0; j<wbobb[i]; j++)  cbobb[i][j]=bobb[i][j];
								wcbobb[i]=wbobb[i];
								
								for (j=0; j<wbobbi[i]; j++)	cbobbi[i][j]=bobbi[i][j];
								wcbobbi[i]=wbobbi[i];
												
							}
		    					fe_en[0]=rsphere_en[0];	
		    					fe_en[1]=rsphere_en[1];
		    					en1_arn=en1_aro;
		    					en2_arn=en2_aro;
							pjet=pjet+1;		
							
				    			
				    			ntt=N;
				    			
				    			k=rand()  % ntt;
				    			
			    				rpid1=pid[k];
				    			boss=drand48();
				    			
				    			
						}	
					
		    			
		    			}
		    		
		    		
		    		}
		    		else 
		    		{
		    			
		    			if (ittt==flag_n)
		    			{
			    			
			    			if (flag_npi<flag_npit) {mo=mo-1;} 
			    			else {nc=2; flag_n=flag_n+1;}
			    			
			    			
			    			if (n_check>nlies) {nc=2; flag_n=flag_n+1; mo=mo+1; }
			    			
		    				nc=1;
		    				
		    				bdm=drand48();
		    				if (n_en1==0) bdm=0.6;
		    				else if (n_en2==0) bdm=0.4;
		    				
		    				if (bdm<0.5)
		    				{
		    					k=ensemble1[rand() % n_en1];
		    					ncv1=(float)n_en1-1.0;
		    					ncv2=(float)n_en2+1.0;
		    					ren_id=1;
		    					ennch=log((en1_arn*ncv2)/(en2_arn*(float)n_en1));
		    				}
		    				else
		    				{
		    					k=ensemble2[rand() % n_en2];
		    					ncv1=(float)n_en1+1.0;
		    					ncv2=(float)n_en2-1.0;
		    					ren_id=0;
		    					ennch=log((en2_arn*ncv1)/((float)n_en2*en1_arn));
		    				}
		    				
		    				boss=0.4;
		    				n_check=n_check+1;
		    				
		    				//printf("here it is consistently fired %d\n", nc);
		    				rpid1=pid[k];
		    			}
		    			else 
		    			{
		    				k=rand() %N;
		    			}
		    	
		    		}
		    		
            			//printf("Buddha was born out of something  !!!!!!!!!!!!!!!!!!  %d   %d\n", mo, k);
				 			
            			mxi=ami[k][0];
	    			myi=ami[k][1];
	    			mzi=ami[k][2];
	    			
	    			ren_id1=ensemble_id[k];

	    			l=(int)floor(mxi/bl_en[ren_id1][0]);
	    			m=(int)floor(myi/bl_en[ren_id1][1]);
	    			n=(int)floor(mzi/bl_en[ren_id1][2]);
	    			
				wpbobbi=0;
				wpbobb=0;
				for(i=0; i<E; i++)
				{
					psbobb[i]=6000;
					psbobbi[i]=6000;
				}
		    		
				
				for (i=0; i<4; i++)
				{
					for (j=0; j<E; j++)
					{
						ppres[i][j]=6000;
					}
					wppres[i]=0;
				}
				
				//printf("Buddha was born out of something  !!!!!!!!!!!!!!   %lf\n", boss);
				
				
				
				displacer[0]=0.0;
				displacer[1]=0.0;
				displacer[2]=0.0;
				wdisplacer=0;
				rotator=0.0;
				
				
				
            			if (boss >= 0.5)
				{
					
					//printf("Buddha was born in rotation !!! \n");
					
					ra1=drand48();
					theta=acos(1.0-(2.0*ra1));
					//ra2=drand48();
					ra2 = JKISS() / 4294967296.0;
	       				phi=2.0*piee*ra2;
				
					xf=sin(theta)*cos(phi);
					yf=sin(theta)*sin(phi);
					zf=cos(theta);
					
					theta=2.0*vstep*drand48();
					//if (theta<0.0) theta=tpiee-theta;  
					
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
					
					
					//printf("raamtt4444\n");
					check=0;
					
					//check=0;
					
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
						
						puc=nrdi[k][ip];
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
						if (re1 < -pbc) re1=re1+f_en[ren_id1][0];
						else if (re1 > pbc) re1=re1-f_en[ren_id1][0];
						if (re2 < -pbc) re2=re2+f_en[ren_id1][1];
						else if (re2 > pbc) re2=re2-f_en[ren_id1][1];
						if (re3 < -pbc) re3=re3+f_en[ren_id1][2];
						else if (re3 > pbc) re3=re3-f_en[ren_id1][2];
						
						
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
							
							ri=0;
							dive=0.0;
							bvs=0;
							bvsi=0;
							divea=0.0;
							divei=0.0;
							if (pop[rpid1][rpid2]==1)
								
							{
								
								zt1=(re1*w7+re2*w8+re3*w9)/mrp;
								zt2=(re1*u7+re2*u8+re3*u9)/mrp;
								if ( zt1 > 0.0 && zt2 < 0.0 )
								{
									bvs=1;
									dive=epsi1[rpid1][rpid2];
									divea=dive;
									
								}
									
								
							
							}
							
							
							else if (pop[rpid1][rpid2]==2)
							
							{
								zt1=(re1*w7+re2*w8+re3*w9)/mrp;
								zt2=-(re1*u7+re2*u8+re3*u9)/mrp;
								if ( zt1 > mom[rpid1] && zt1 < maom[rpid1] && zt2 > mom[rpid2] && zt2 < maom[rpid2] )
								{
							
									bvs=1;
									dive=epsi2[rpid1][rpid2];
									divea=dive;
									
								}
							}
							
							
							else if (pop[rpid1][rpid2]==3)
							{
								
								bvsi=1;
								dive=epsi3[rpid1][rpid2];
								divei=dive;
							}
							
							else if (pop[rpid1][rpid2]==4)
							
							{
							
								
								zt1=(re1*w7+re2*w8+re3*w9)/mrp;
								zt2=(re1*u7+re2*u8+re3*u9)/mrp;
								if ( zt1 > 0.0 && zt2 < 0.0 )
								{
									bvs=1;
									divea=epsi1[rpid1][rpid2];
									
									
								}
									
								
								bvsi=1;
								dive=epsi;
								divei=epsi3[rpid1][rpid2];
								
							
							}
							
							
							else if (pop[rpid1][rpid2]==5)
							
							{
																			
								zt1=(fabs(re1*w7+re2*w8+re3*w9))/mrp;
								zt2=(fabs(re1*u7+re2*u8+re3*u9))/mrp;
								if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
								{
									bvs=1;
									divea=epsi2[rpid1][rpid2];
									
									
								}
								
								bvsi=1;
								dive=epsi;
								divei=epsi3[rpid1][rpid2];
								
							
							}
							
							
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
										
									}
									else
									{
										rotator=rotator+fabs(w7*u7+w8*u8+w9*u9)-fabs(li[k][6]*u7+li[k][7]*u8+li[k][8]*u9);
										
										
									}
									
									
								}
							
								else
								{	
									ep=ep/mr;
									
									ec=ep*(mrp-divea)*(mrp-divea);
									if (bvs==1 && ec<1.0)
									{
										psbobb[wpbobb]=puc;
										wpbobb=wpbobb+1;	
										
										
									}
									
									ec=ep*(mrp-divei)*(mrp-divei);
									if (bvsi==1 && ec<1.0)
									{
									
										psbobbi[wpbobbi]=puc;
										wpbobbi=wpbobbi+1;	
										
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
							if (rotator<0.00000001)	check=1;
						}
						
						
						ene=0.0;
						for (it=0; it<wcbobb[k]; it++)
						{
							rpid2=pid[cbobb[k][it]];
							ene=ene+ur[rpid1][rpid2];
							
						}	
						
						for (it=0; it<wcbobbi[k]; it++)
						{
							rpid2=pid[cbobbi[k][it]];
							ene=ene+uri[rpid1][rpid2];
							
							
							
						}	
						
						bdm=0.0;
						for (i=0; i < wpbobb; i++)
						{
							
							rpid2=pid[psbobb[i]];
							bdm=bdm+ur[rpid1][rpid2];
							
						}
						
						for (i=0; i < wpbobbi; i++)
						{
							
							rpid2=pid[psbobbi[i]];
							bdm=bdm+uri[rpid1][rpid2];
							
						}
						
						pre=(bdm-ene);
						
						
						if (pre<0.0) bdm=1.1;
						else bdm=exp(-pre);
						
						if (drand48()>bdm) {check=1;}
						
					}	
						
					if (check==0) 
					{
						
						mum=wcbobb[k];
							
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobb[k][jt];
							drum=wcbobb[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobb[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobb[chum][tric]=cbobb[chum][(tric+1)];
									}
									break;
								}
							
							}
							wcbobb[chum]=drum-1;
							cbobb[k][jt]=6000;
						}
							
						wcbobb[k]=0;	
						for (jdt=0; jdt < wpbobb; jdt++)
						{	
							//printf("raam is ");
							td=psbobb[jdt];
							cbobb[k][jdt]=td;
								tum=wcbobb[td];
								cbobb[td][tum]=k;
							wcbobb[td]=tum+1;
							wcbobb[k]=wcbobb[k]+1;
						}
							
						mum=wcbobbi[k];
						
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobbi[k][jt];
							drum=wcbobbi[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobbi[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobbi[chum][tric]=cbobbi[chum][(tric+1)];
									}
									break;
								}
							
							}
							wcbobbi[chum]=drum-1;
							cbobbi[k][jt]=6000;
						}
							
						wcbobbi[k]=0;	
						for (jdt=0; jdt < wpbobbi; jdt++)
						{	
							//printf("raam is ");
							td=psbobbi[jdt];
							cbobbi[k][jdt]=td;
								tum=wcbobbi[td];
								cbobbi[td][tum]=k;
							wcbobbi[td]=tum+1;
							wcbobbi[k]=wcbobbi[k]+1;
						}
						
						
						for (jt=0; jt<9; jt++)
						{
							li[k][jt]=dli[jt];
							el[k][jt]=del[jt];
							
							
						}
						
						
						lom[k][0]=lom[k][0]+dom0;
						lom[k][1]=lom[k][1]+dom1;
						lom[k][2]=lom[k][2]+dom2;
						
					}
					//printf("raamtttt\n");  	
				}
                           	
                		else
				{
					//printf("Buddha was born in translation !!! \n");          			
					trim=2;
					
	        			
					
					
					
					//bdm=(mxi-cm_en[ren_id1][0])*(mxi-cm_en[ren_id1][0])+(myi-cm_en[ren_id1][1])*(myi-cm_en[ren_id1][1])+(mzi-cm_en[ren_id1][2])*(mzi-cm_en[ren_id1][2]);
					
					//pre=(xi-cm_en[ren_id1][0])*(xi-cm_en[ren_id1][0])+(yi-cm_en[ren_id1][1])*(yi-cm_en[ren_id1][1])+(zi-cm_en[ren_id1][2])*(zi-cm_en[ren_id1][2]);
					//printf("%d  %d  %d  %lf   %lf   %lf\n", vc, nc, sqrt(bdm), sqrt(pre), rsphere_en[ren_id1]);
					
					
					
					if (nc==1)
					{
						ren_id1=ren_id;
						
						dx=0.0;
						dy=0.0;
						dz=0.0;
						ra1= JKISS() / 4294967296.0;
						theta=acos(1.0-(2.0*ra1));
						ra2=drand48();
				       		phi=tpiee*ra2;
						
						pdx=sin(theta)*cos(phi);
						pdy=sin(theta)*sin(phi);
						pdz=cos(theta);
						
						xi=rsphere_en[ren_id1]*pdx+cm_en[ren_id1][0];
				    		yi=rsphere_en[ren_id1]*pdy+cm_en[ren_id1][1];
				    		zi=rsphere_en[ren_id1]*pdz+cm_en[ren_id1][2];
			
					}
					else if (vc==1)
					{
						dx=0.0;
						dy=0.0;
						dz=0.0;
						
						xi=dx+ami[k][0];	
						yi=dy+ami[k][1];
						zi=dz+ami[k][2];
					}
					else
					{
						rdth=vstep;
						if (cow==0) rdth=0.01;
						else tsjet=tsjet+1;
						
						rdth=rdth*drand48();
						ra2 = JKISS() / 4294967296.0;
		       				rdphi=tpiee*ra2;
		       				
		       				xf=sin(rdth)*cos(rdphi);
						yf=sin(rdth)*sin(rdphi);
						zf=cos(rdth);
					//	printf("Here is the communication:   %lf   %lf   %lf\n", xf, yf, zf);
						xi=(mxi-cm_en[ren_id1][0])/rsphere_en[ren_id1];
						yi=(myi-cm_en[ren_id1][1])/rsphere_en[ren_id1];
						zi=(mzi-cm_en[ren_id1][2])/rsphere_en[ren_id1];
					//	printf("Here is the communication:   %lf   %lf   %lf\n", xi, yi, zi);
						theta=acos(zi);
						if (fabs(zi) > 0.99999)
						{
							dx=1.0;
							dy=0.0;
							dz=0.0;
						}
						else
						{
							dx=-yi*(1.0/pow(((xi*xi)+(yi*yi)),0.5));
							dy=xi*(1.0/pow(((xi*xi)+(yi*yi)),0.5));
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
						
						xi=(xf*q11+yf*q12+zf*q13)*rsphere_en[ren_id1]+cm_en[ren_id1][0];
						yi=(xf*q21+yf*q22+zf*q23)*rsphere_en[ren_id1]+cm_en[ren_id1][1];
						zi=(xf*q31+yf*q32+zf*q33)*rsphere_en[ren_id1]+cm_en[ren_id1][2];
							
						dx=xi-mxi;
						dy=yi-myi;
						dz=zi-mzi;
								
					}
					
					//bdm=(mxi-cm_en[ren_id1][0])*(mxi-cm_en[ren_id1][0])+(myi-cm_en[ren_id1][1])*(myi-cm_en[ren_id1][1])+(mzi-cm_en[ren_id1][2])*(mzi-cm_en[ren_id1][2]);
					
					//pre=(xi-cm_en[ren_id1][0])*(xi-cm_en[ren_id1][0])+(yi-cm_en[ren_id1][1])*(yi-cm_en[ren_id1][1])+(zi-cm_en[ren_id1][2])*(zi-cm_en[ren_id1][2]);
					
					
					//printf("%d  %d  %lf   %lf   %lf\n", vc, nc, sqrt(bdm), sqrt(pre), rsphere_en[ren_id1]);
					
	        			if (xi >= f_en[ren_id1][0]) {xi=xi-f_en[ren_id1][0];}
	        			else if (xi < 0.0) {xi=xi+f_en[ren_id1][0];}
					if (yi >= f_en[ren_id1][1]) {yi=yi-f_en[ren_id1][1];}
	        			else if (yi < 0.0) {yi=yi+f_en[ren_id1][1];}
					if (zi >= f_en[ren_id1][2]) {zi=zi-f_en[ren_id1][2];}
	        			else if (zi < 0.0) {zi=zi+f_en[ren_id1][2];}
	        			
	        		///	printf("snn7rrrrrrrrrrrrrrrr50 %d\n",mo);
	        			
	        			
	        			xk=(int)floor(xi/bl_en[ren_id1][0]);
	        			yk=(int)floor(yi/bl_en[ren_id1][1]);
	        			zk=(int)floor(zi/bl_en[ren_id1][2]);
				//	printf("raaaam2 %d %d %d %d %lf %lf %lf  %d %d %d\n",  nc,  xk,   yk,   zk,   xi,   yi,  zi, block_en[0][0], block_en[0][1], block_en[0][2] );
						
	    				for (j=0; j<UN; j++)	pit[j]=6000;
					wpi=0;	
					
					for(jt=0;jt<27; jt++)
					{
						
						nxk=xk+tnkx[jt][0];
						nyk=yk+tnkx[jt][1];
						nzk=zk+tnkx[jt][2];
						
						
						if (nxk >= block_en[ren_id1][0])
						{
			    				nxk=0;
			   	 			cpsh=f_en[ren_id1][0];
						}
						
						else if (nxk==-1)
						{
			    				nxk=nxk+block_en[ren_id1][0];
			    				cpsh=-f_en[ren_id1][0];
						}
						else 
						{	
							cpsh=0.0;
						}
						
						
						if (nyk >= block_en[ren_id1][1])
		       			{    
							nyk=0;
			   			 	cph=f_en[ren_id1][1];
						}
						else if (nyk==-1)
						{
			    				nyk=nyk+block_en[ren_id1][1];
			    				cph=-f_en[ren_id1][1];
						}
						else 
						{	
							cph=0.0;
						}
						
						
						
						if (nzk >= block_en[ren_id1][2])
		       			{    
							nzk=0;
			   			 	cp=f_en[ren_id1][2];
						}
						else if (nzk==-1)
						{
			    				nzk=nzk+block_en[ren_id1][2];
			    				cp=-f_en[ren_id1][2];
						}
						else 
						{	
							cp=0.0;
						}
						
						
 						
						dum=wdr[nxk][nyk][nzk];
						for (jdt=0; jdt<dum; jdt++)
						{
							da=dri[nxk][nyk][nzk][jdt];  
							rpid2=pid[da];         
							xf=ami[da][0]+cpsh;
							yf=ami[da][1]+cph;
							zf=ami[da][2]+cp;
							if (da != k && ren_id1==ensemble_id[da])
							{
				    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi)); 
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
											pit[wpi]=da;
											wpi=wpi+1;
										}
						   				
										
									}
									else
									{
										pit[wpi]=da;
										wpi=wpi+1;
									
									}
								}
							}
						} 
						
						
						if (trim==1) break;
						
					
					}
						
				//	printf("sn111111111111111n750 %d\n",mo);
					
					
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
							
							puc=pit[ip];
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
							if (re1 < -pbc) re1=re1+f_en[ren_id1][0];
							else if (re1 > pbc) re1=re1-f_en[ren_id1][0];
							if (re2 < -pbc) re2=re2+f_en[ren_id1][1];
							else if (re2 > pbc) re2=re2-f_en[ren_id1][1];
							if (re3 < -pbc) re3=re3+f_en[ren_id1][2];
							else if (re3 > pbc) re3=re3-f_en[ren_id1][2];
							
							
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
								if (pop[rpid1][rpid2]==1)
									
								{
									
									zt1=(re1*w7+re2*w8+re3*w9)/mrp;
									zt2=(re1*u7+re2*u8+re3*u9)/mrp;
									if ( zt1 > 0.0 && zt2 < 0.0 )
									{
										bvs=1;
										dive=epsi1[rpid1][rpid2];
										divea=dive;
										
									}
										
									
								
								}
								
								
								else if (pop[rpid1][rpid2]==2)
								{
									
									
									zt1=(re1*w7+re2*w8+re3*w9)/mrp;
									zt2=-(re1*u7+re2*u8+re3*u9)/mrp;
									if ( zt1 > mom[rpid1] && zt1 < maom[rpid1] && zt2 > mom[rpid2] && zt2 < maom[rpid2] )
									{
									
										bvs=1;
										dive=epsi2[rpid1][rpid2];
										divea=dive;
										
									}
									
								}
								
								
								else if (pop[rpid1][rpid2]==3)
								{
									
									bvsi=1;
									dive=epsi3[rpid1][rpid2];
									divei=dive;
								}
								
								else if (pop[rpid1][rpid2]==4)
								
								{
								
									
									zt1=(re1*w7+re2*w8+re3*w9)/mrp;
									zt2=(re1*u7+re2*u8+re3*u9)/mrp;
									if ( zt1 > 0.0 && zt2 < 0.0 )
									{
										bvs=1;
										divea=epsi1[rpid1][rpid2];
										
										
									}
										
									
									bvsi=1;
									dive=epsi;
									divei=epsi3[rpid1][rpid2];
									
								
								}
								
								
								else if (pop[rpid1][rpid2]==5)
								
								{
																				
									zt1=(fabs(re1*w7+re2*w8+re3*w9))/mrp;
									zt2=(fabs(re1*u7+re2*u8+re3*u9))/mrp;
									if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
									{
										bvs=1;
										divea=epsi2[rpid1][rpid2];
										
										
									}
									
									bvsi=1;
									dive=epsi;
									divei=epsi3[rpid1][rpid2];
									
								
								}
								
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
								
									else
									{	
										ep=ep/mr;
										
										
										
										
										ec=ep*(mrp-divea)*(mrp-divea);
										if (bvs==1 && ec<1.0)
										{
											psbobb[wpbobb]=puc;
											wpbobb=wpbobb+1;	
											
											
										}
										
										ec=ep*(mrp-divei)*(mrp-divei);
										if (bvsi==1 && ec<1.0)
										{
										
											psbobbi[wpbobbi]=puc;
											wpbobbi=wpbobbi+1;	
											
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
					
				//	printf("snn750ffffffffffffff %d\n",mo);
					
					
					if (trim==2 && vc!=1)
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
						
					
						ene=0.0;
						for (it=0; it<wcbobb[k]; it++)
						{
							rpid2=pid[cbobb[k][it]];
							ene=ene+ur[rpid1][rpid2];
						}	
						
						for (it=0; it<wcbobbi[k]; it++)
						{
							rpid2=pid[cbobbi[k][it]];
							ene=ene+uri[rpid1][rpid2];
						}	
						
						//printf
						
						bdm=0.0;
						for (i=0; i < wpbobb; i++)
						{
							
							rpid2=pid[psbobb[i]];
							bdm=bdm+ur[rpid1][rpid2];
							
						}
						
						for (i=0; i < wpbobbi; i++)
						{
							
							rpid2=pid[psbobbi[i]];
							bdm=bdm+uri[rpid1][rpid2];
							
						}
						
						
						pre=(bdm-ene);
						
						if (nc==1)
						{
							//printf("During Particle change  %d %lf  \n", sitt,   pre);
							pre=pre+ennch;
						}
							
						if (pre<0.0) bdm=1.1;
						else bdm=exp(-pre);
						
						//printf("%d %lf \n", trim, bdm);
						
						
						if (drand48()>bdm) {trim=1;}
						
						//printf("%d \n", trim);
						
					}
					
					
				//	printf("snn750dfgdf %d\n",mo);
					
					if (trim==2)
					{
						if (vc != 1 && nc != 1) sjet=sjet+1;
						if (nc==1)
						{
							if (ensemble_id[k]==0)
							{
								for (i=0; i<n_en1; i++)
								{
									if(ensemble1[i]==k)
									{
										for (j=i; j<n_en1; j++)
										ensemble1[j]=ensemble1[j+1];
										n_en1=n_en1-1;
										break;
									}
								
								}
								ensemble_id[k]=1;
								ensemble2[n_en2]=k;
								n_en2=n_en2+1;
							}
							else
							{
								for (i=0; i<n_en2; i++)
								{
									if(ensemble2[i]==k)
									{
										for (j=i; j<n_en2; j++)
										ensemble2[j]=ensemble2[j+1];
										n_en2=n_en2-1;
										break;
									}
								
								}
								ensemble_id[k]=0;
								ensemble1[n_en1]=k;
								n_en1=n_en1+1;
							}
							flag_npi=flag_npi+1;
							
							//printf("It is running but not the others,\n");
							
						}
						
						dum=wdr[l][m][n];
						for (j=0;j<dum;j++)
						{	
							
							if (dri[l][m][n][j]==k)
							{
								for (tric=j; tric<dum; tric++)	dri[l][m][n][tric]=dri[l][m][n][tric+1];
								break;
							}
									
						}
						wdr[l][m][n]=wdr[l][m][n]-1;
						dum=wdr[xk][yk][zk];
						dri[xk][yk][zk][dum]=k;
						wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
						
						ami[k][0]=xi;
		    				
		    				ami[k][1]=yi;
		    				
		    				ami[k][2]=zi;
		    				
						ldcxyz[k][0]=ldcxyz[k][0]+dx;
						ldcxyz[k][1]=ldcxyz[k][1]+dy;
						ldcxyz[k][2]=ldcxyz[k][2]+dz;
						
						mum=wnrd[k];
						
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=nrdi[k][jt];
							drum=wnrd[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (nrdi[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)	nrdi[chum][tric]=nrdi[chum][(tric+1)];
									
									break;
								}
							
							}
							wnrd[chum]=drum-1;
							nrdi[k][jt]=6000;
						}
							
						wnrd[k]=0;	
						for (jdt=0; jdt < wpi; jdt++)
						{	
							//printf("raam is ");
							td=pit[jdt];
							nrdi[k][jdt]=td;
								tum=wnrd[td];
								nrdi[td][tum]=k;
							wnrd[td]=tum+1;
							wnrd[k]=wnrd[k]+1;
						}	
						
						
						
						mum=wcbobb[k];
							
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobb[k][jt];
							drum=wcbobb[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobb[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobb[chum][tric]=cbobb[chum][(tric+1)];
									}
									break;
								}
							
							}
							wcbobb[chum]=drum-1;
							cbobb[k][jt]=6000;
						}
							
						wcbobb[k]=0;	
						for (jdt=0; jdt < wpbobb; jdt++)
						{	
							//printf("raam is ");
							td=psbobb[jdt];
							cbobb[k][jdt]=td;
								tum=wcbobb[td];
								cbobb[td][tum]=k;
							wcbobb[td]=tum+1;
							wcbobb[k]=wcbobb[k]+1;
						}
						
						
						
						mum=wcbobbi[k];
						
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=cbobbi[k][jt];
							drum=wcbobbi[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (cbobbi[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										cbobbi[chum][tric]=cbobbi[chum][(tric+1)];
									}
									break;
								}
							
							}
							wcbobbi[chum]=drum-1;
							cbobbi[k][jt]=6000;
						}
							
						wcbobbi[k]=0;	
						for (jdt=0; jdt < wpbobbi; jdt++)
						{	
							//printf("raam is ");
							td=psbobbi[jdt];
							cbobbi[k][jdt]=td;
								tum=wcbobbi[td];
								cbobbi[td][tum]=k;
							wcbobbi[td]=tum+1;
							wcbobbi[k]=wcbobbi[k]+1;
						}
						
						
						
						for (i=0; i<4; i++)
						{
							mum=wpres[k][i];
								
							for (jt=0; jt<mum; ++jt)
		    					{
								chum=pres[k][i][jt];
								drum=wpres[chum][i];
								for (sj=0; sj<drum; ++sj)
								{
									if (pres[chum][i][sj]==k)
									{
										for (tric=sj; tric<drum; tric++)
										{
											pres[chum][i][tric]=pres[chum][i][(tric+1)];
										}
										break;
									}
								
								}
								wpres[chum][i]=drum-1;
								pres[k][i][jt]=6000;
							}
								
							wpres[k][i]=0;	
							for (jdt=0; jdt < wppres[i]; jdt++)
							{	
								//printf("raam is ");
								td=ppres[i][jdt];
								pres[k][i][jdt]=td;
									tum=wpres[td][i];
									pres[td][i][tum]=k;
								wpres[td][i]=tum+1;
								wpres[k][i]=wpres[k][i]+1;
							}
						
						}

						
					}
					else
					{
						
						trims=0; //printf("raam %d %d \n", ittt, k);
					}
				//	printf("snn750 %d\n",mo);
				
	
				}
				
				
			//printf("snn750 %d\n",mo);
			}
			
			
			
			//printf("Buddhdha was born at the number %d\n", ittt);
		//printf("%lf  %lf\n", (f_en[0][1]*f_en[0][1]*f_en[0][1])+(f_en[1][1]*f_en[1][1]*f_en[1][1]), g_vol);
	//	if (cow==1)
	//	{	
	//	printf("Simulation is going on itt %d  %d  %lf  %lf  %lf    %lf  Rho = %lf  \n", ittt, n_en1, fer_en[0], f_en[0][0], f_en[0][1], f_en[0][2], (float)n_en1/(f_en[0][0]*f_en[0][1]*f_en[0][2]));
		
	//	printf("Simulation is going on itt %d  %d  %lf  %lf  %lf    %lf  Rho = %lf  \n\n", ittt, n_en2, fer_en[1], f_en[1][0], f_en[1][1], f_en[1][2], (float)n_en2/(f_en[1][0]*f_en[1][1]*f_en[1][2]));
		
	//	}	
						
		}
		printf("Simulation is going on itt %d  %d     %lf  Rho = %lf  \n", sitt, n_en1, rsphere_en[0], (float)n_en1/(4.0*piee*rsphere_en[0]*rsphere_en[0]));
		
		printf("Simulation is going on itt %d  %d     %lf  Rho = %lf  \n\n", sitt, n_en2, rsphere_en[1], (float)n_en2/(4.0*piee*rsphere_en[1]*rsphere_en[1]));
		
		bdm=(float)sjet/(float)tsjet;
		
		if (bdm<0.3) vstep=vstep*0.95;
		else if (bdm > 0.4)  vstep=vstep*1.05;
		
		
	
		bdm=(float)pjet/(float)tpjet;
		
		if (bdm<0.3) vpstep=vpstep*0.95;
		else if (bdm > 0.4)  vpstep=vpstep*1.05;
		
		
		
		if(vstep>0.4) vstep=0.4;
		else if (vstep<0.001) vstep=0.001;
		
		if(vpstep>0.4) vpstep=0.4;
		else if (vpstep<0.001) vpstep=0.001;
		
			
		mxi=0.0;
		myi=0.0;
		mzi=0.0;
		for (k=0; k<4; k++)
		{
			for(i=0; i<N; i++)
			{
				
				ren_id1=ensemble_id[i];
				if (ren_id1==0)
				{
				
					for(j=0; j<wpres[i][k]; j++)
					{
						
						da=pres[i][k][j];
						
						//printf("%d      %d     %d\n", i, da, k);
						xi=ami[i][0]-ami[da][0];
						if (xi>pbc) xi=xi-f_en[ren_id1][0];
						else if(xi<-pbc) xi=xi+f_en[ren_id1][0];
						
						yi=ami[i][1]-ami[da][1];
						if (yi>pbc) yi=yi-f_en[ren_id1][1];
						else if(yi<-pbc) yi=yi+f_en[ren_id1][1];
						
						zi=ami[i][2]-ami[da][2];
						if (zi>pbc) zi=zi-f_en[ren_id1][2];
						else if(zi<-pbc) zi=zi+f_en[ren_id1][2];
						
						rpid1=pid[i];
						rpid2=pid[da];
						
						mr=sqrt(xi*xi+yi*yi+zi*zi);
						if(k==0) 
						{
							ep=urh[rpid1][rpid2];
						}
						else if(k==3)
						{
							ep=uri[rpid1][rpid2];
						}
						else
						{
							ep=ur[rpid1][rpid2];
						}
						
						if (da>i)
						{
							mxi=mxi+((xi*xi*ep)/(mr*dlr));
							myi=myi+((yi*yi*ep)/(mr*dlr));
							mzi=mzi+((zi*zi*ep)/(mr*dlr));
						}
						//printf("%lf\n", ep);
						
					}
				}
					
			}
		}
		ren_id1=0;
		ec=((float)n_en1/(f_en[ren_id1][0]*f_en[ren_id1][1]*f_en[ren_id1][2]))+(mxi+myi+mzi)/(3.0*f_en[ren_id1][0]*f_en[ren_id1][1]*f_en[ren_id1][2]);
		//ec=ec+((float)N/(f_en[ren_id1][0]*f[1]*f[2]));
		FILE* pu;
		pu=fopen("Pressure_histrogram.txt","a");
		fprintf(pu,"%lf\n", ec);
		
		
		
		
		mxi=0.0;
		myi=0.0;
		mzi=0.0;
		for (k=0; k<4; k++)
		{
			for(i=0; i<N; i++)
			{
				
				ren_id1=ensemble_id[i];
				if (ren_id1==1)
				{
				
					for(j=0; j<wpres[i][k]; j++)
					{
						
						da=pres[i][k][j];
						
						//printf("%d      %d     %d\n", i, da, k);
						xi=ami[i][0]-ami[da][0];
						if (xi>pbc) xi=xi-f_en[ren_id1][0];
						else if(xi<-pbc) xi=xi+f_en[ren_id1][0];
						
						yi=ami[i][1]-ami[da][1];
						if (yi>pbc) yi=yi-f_en[ren_id1][1];
						else if(yi<-pbc) yi=yi+f_en[ren_id1][1];
						
						zi=ami[i][2]-ami[da][2];
						if (zi>pbc) zi=zi-f_en[ren_id1][2];
						else if(zi<-pbc) zi=zi+f_en[ren_id1][2];
						
						rpid1=pid[i];
						rpid2=pid[da];
						
						mr=sqrt(xi*xi+yi*yi+zi*zi);
						if(k==0) 
						{
							ep=urh[rpid1][rpid2];
						}
						else if(k==3)
						{
							ep=uri[rpid1][rpid2];
						}
						else
						{
							ep=ur[rpid1][rpid2];
						}
						
						if (da>i)
						{
							mxi=mxi+((xi*xi*ep)/(mr*dlr));
							myi=myi+((yi*yi*ep)/(mr*dlr));
							mzi=mzi+((zi*zi*ep)/(mr*dlr));
						}
						//printf("%lf\n", ep);
						
					}
				}
					
			}
		}
		ren_id1=1;
		ec=((float)n_en2/(f_en[ren_id1][0]*f_en[ren_id1][1]*f_en[ren_id1][2]))+(mxi+myi+mzi)/(3.0*f_en[ren_id1][0]*f_en[ren_id1][1]*f_en[ren_id1][2]);
		//ec=ec+((float)N/(f[0]*f[1]*f[2]));
		
		fprintf(pu,"%lf\n", ec);
		fclose(pu);
		
		
		
		
		
		
		
		char filename[25] = {0};

    		sprintf(filename, "%d.txt", (int)floor((float)(sitt+1)/100.0));
		
		FILE* cfu;
		cfu=fopen(filename,"a");
			
		fprintf(cfu,"%lf   ", f_en[0][0]);
		fprintf(cfu,"%lf   ", f_en[0][1]);
		fprintf(cfu,"%lf\n", f_en[0][2]);
		
		fprintf(cfu,"%lf   ", f_en[1][0]);
		fprintf(cfu,"%lf   ", f_en[1][1]);
		fprintf(cfu,"%lf\n", f_en[1][2]);
		
		fprintf(cfu,"\n");
		for (i=0; i<N; i++)
		{
			fprintf(cfu,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
			
			fprintf(cfu, "%d\n", wcbobb[i]);
			for (j=0; j<wcbobb[i]; j++) 
			{
				fprintf(cfu,"   %d", cbobb[i][j]);
			}
			fprintf(cfu, "\n");
			
			
			fprintf(cfu, "%d\n", wcbobbi[i]);
			for (j=0; j<wcbobbi[i]; j++) 
			{
				fprintf(cfu,"   %d", cbobbi[i][j]);
			}
			fprintf(cfu, "\n");
			
			
			fprintf(cfu,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
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
					
					difus=difus+(ldcxyz[i][0]*ldcxyz[i][0]+ldcxyz[i][1]*ldcxyz[i][1]+ldcxyz[i][2]*ldcxyz[i][2]);
					difusl=difusl+((lom[i][0]*lom[i][0])+(lom[i][1]*lom[i][1])+(lom[i][2]*lom[i][2]));
				}
				
				
				difus=difus/(float)nt[j];
				
				difusl=(difusl*rdtheta[j]*rdtheta[j])/(float)nt[j];
				
				
				difus=difus/(float)(sitt-isitt);
				
				difusl=difusl/(float)(sitt-isitt);
				
				
				
				printf("Large:  %d   %lf %lf %d  \n", (sitt+1), difus, difusl, jet);
				
				FILE* mmwwu;
     				mmwwu=fopen("Diffusion_coefficient_el_l.txt","a");
				
				fprintf(mmwwu,"%d   %lf  %lf  \n", (sitt+1), difus, difusl);
				
				fclose(mmwwu);
				
				
				jt=jtt;
				
			}
			
		}
		
		FILE* fu;
	     	fu=fopen("Cord_ensemble_1.vtk","w");
		fprintf(fu,"# vtk DataFile Version 3.0\n");
		fprintf(fu,"Random data to test tensors\n");
		fprintf(fu,"ASCII\n");
		fprintf(fu,"DATASET POLYDATA\n");
		fprintf(fu,"POINTS %d float\n", (int)(n_en1));
		
		
		FILE* fu1;
	     	fu1=fopen("Cord_ensemble_2.vtk","w");
		fprintf(fu1,"# vtk DataFile Version 3.0\n");
		fprintf(fu1,"Random data to test tensors\n");
		fprintf(fu1,"ASCII\n");
		fprintf(fu1,"DATASET POLYDATA\n");
		fprintf(fu1,"POINTS %d float\n", (int)(n_en2));
		
		
		for (jt=0; jt<N; jt++)
		{
			if (ensemble_id[jt]==0)
			{
				
			  	fprintf(fu,"%lf    %lf    %lf  \n",ami[jt][0], ami[jt][1], ami[jt][2]);
			}
			else
			{
				fprintf(fu1,"%lf    %lf    %lf  \n",ami[jt][0]+f_en[0][0]+3.0, ami[jt][1], ami[jt][2]);	
			}
		}

		fprintf(fu,"\n");
		fprintf(fu,"POINT_DATA %d\n", (int)(n_en1));
		fprintf(fu,"\n");
		fprintf(fu,"\n");
		fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
		fprintf(fu,"\n");
		
		fprintf(fu1,"\n");
		fprintf(fu1,"POINT_DATA %d\n", (int)(n_en2));
		fprintf(fu1,"\n");
		fprintf(fu1,"\n");
		fprintf(fu1,"TENSORS non-spherical_ellipsoid float\n");
		fprintf(fu1,"\n");
		
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
			
			if (ensemble_id[i]==0)
			{
				fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
				fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
				fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
				fprintf(fu,"\n");
				
			}
			else
			{
				fprintf(fu1,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
				fprintf(fu1,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
				fprintf(fu1,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
				fprintf(fu1,"\n");
			
			
			}
			
		}
		
		
		fclose(fu1);
		fclose(fu);
		
		
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
			
			fprintf(ku,"%lf\n", quao);
			
			fprintf(ku,"%lf\n", quaoi);
			
				
			fprintf(ku,"%lf   ", f_en[0][0]);
			fprintf(ku,"%lf   ", f_en[0][1]);
			fprintf(ku,"%lf\n", f_en[0][2]);
			
			fprintf(ku,"%lf   ", f_en[1][0]);
			fprintf(ku,"%lf   ", f_en[1][1]);
			fprintf(ku,"%lf\n", f_en[1][2]);
			
			for (i=0; i<N; i++)
			{
				
				
				fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
				
				fprintf(ku, "%d\n", ensemble_id[i]);
				
				fprintf(ku, "%d\n", wnrd[i]);
				for (j=0; j<wnrd[i]; j++) 
				{
					fprintf(ku,"   %d", nrdi[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku, "%d\n", wcbobb[i]);
				for (j=0; j<wcbobb[i]; j++) 
				{
					fprintf(ku,"   %d", cbobb[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku, "%d\n", wcbobbi[i]);
				for (j=0; j<wcbobbi[i]; j++) 
				{
					fprintf(ku,"   %d", cbobbi[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
				
			}
			
			fclose(ku);
			
			
			
			
			
			//VTK Zone starts!
	
			FILE* fu;
		     	fu=fopen("Cord_ensemble_1.vtk","w");
			fprintf(fu,"# vtk DataFile Version 3.0\n");
			fprintf(fu,"Random data to test tensors\n");
			fprintf(fu,"ASCII\n");
			fprintf(fu,"DATASET POLYDATA\n");
			fprintf(fu,"POINTS %d float\n", (int)(n_en1));
			
			
			FILE* fu1;
		     	fu1=fopen("Cord_ensemble_2.vtk","w");
			fprintf(fu1,"# vtk DataFile Version 3.0\n");
			fprintf(fu1,"Random data to test tensors\n");
			fprintf(fu1,"ASCII\n");
			fprintf(fu1,"DATASET POLYDATA\n");
			fprintf(fu1,"POINTS %d float\n", (int)(n_en2));
			
			
			for (jt=0; jt<N; jt++)
			{
				if (ensemble_id[jt]==0)
				{
					
				  	fprintf(fu,"%lf    %lf    %lf  \n",ami[jt][0], ami[jt][1], ami[jt][2]);
				}
				else
				{
					fprintf(fu1,"%lf    %lf    %lf  \n",ami[jt][0]+f_en[0][0]+3.0, ami[jt][1], ami[jt][2]);	
				}
			}

			fprintf(fu,"\n");
			fprintf(fu,"POINT_DATA %d\n", (int)(n_en1));
			fprintf(fu,"\n");
			fprintf(fu,"\n");
			fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
			fprintf(fu,"\n");
			
			fprintf(fu1,"\n");
			fprintf(fu1,"POINT_DATA %d\n", (int)(n_en2));
			fprintf(fu1,"\n");
			fprintf(fu1,"\n");
			fprintf(fu1,"TENSORS non-spherical_ellipsoid float\n");
			fprintf(fu1,"\n");
			
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
				
				if (ensemble_id[i]==0)
				{
					fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
					fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
					fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
					fprintf(fu,"\n");
					
				}
				else
				{
					fprintf(fu1,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
					fprintf(fu1,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
					fprintf(fu1,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
					fprintf(fu1,"\n");
				
				
				}
				
			}
			
			
			fclose(fu1);
			fclose(fu);
			

			//VTK Zone ends!

			
			
			
			
			
		}
		
		if (((float)(sitt+1)/1000.0)==((float)floor((float)(sitt+1)/1000.0)))
		{
		
		
			qua=0.0;
	
			for (i=0; i<N; i++)
			{
				qua=qua+(float)wbobb[i];
			
			}
			qua=qua/(float)N;
			
			quai=0.0;
	
			for (i=0; i<N; i++)
			{
				quai=quai+(float)wbobbi[i];
			
			}
			quai=quai/(float)N;

			
				if (fabs(qua-quao) < 0.05 && fabs(quai-quaoi) < 0.05) { printf("please don't run leave me alone\n"); break;}
				quao=qua;
				quaoi=quai;
		
		
		
		}
		
		
		
		

	}
	
	
	
	
	
	
	
	
	
	
	fclose (mmwwu);
	
	fclose (cfu);

	//SImulaton ends!


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
	
	
		
	for (i=0; i<buf_blockd; i++)
	{
		
		for (j=0; j<buf_blockd; j++)
		{
		
			for (k=0; k<buf_blockd; k++)
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





