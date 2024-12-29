#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <omp.h>
const int E=20, UN=600, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel;

//void ecfcal(int,double[]);

int  N, trim, check, checkl, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, it, nf,  block[3], blockd, jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct, rul, sla, veei, veec, SHA, ptr, ptrr, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt, gv, bvs, bvsi, cloop, icloop, tphy1, tphy2, ftaverage, time_interval, target_time, starter_flag, n_sec, n_k_sec,  n_a_sec, vfrap, probation, initia_p, blockob[3], xkob, ykob, zkob, nxkob, nykob, nzkob, pidres, wi, wnlies;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9], f[3], bl1, bl2, bl3, bl[3], f1, f2, f3, ar, per, loopd, blob[3], bleob,sss, qqq, rho, widoms, twidoms, cp, cph, cpsh; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, fe, tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON, eb, teb, tbl, pbc,  dif, difu, difus, difl, difusl, adefus, dive, divea, divei, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm, robs, obsd, ofact, pbcp, pbco, radial_width, radial_range, k_width, k_range, coef_1, coef_2, coef_3, director_ore[3], angle_width, angle_range, pfvf;  

double stretched_fitting(double*, double*, double*, int, double**);

int main(void)
{

	srand48(time(NULL));
	srand(time(NULL));
	
	
	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;
	

	static unsigned long long  x=123456789,y=987654321,z=43219876,cc=6543217; 
	unsigned int JKISS()
	{ 
		unsigned long long t;
		x=314527869*x+1234567; 
		y^=y <<5;y^=y>>7;y^=y<<22;
		t = 4294584393ULL*z+cc; cc = t>>32; z= t;
		return x+y+z; 
	}	
	
	
		
     	srand48(time(NULL));
	srand(time(NULL));
	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	
	fscanf(bu,"%lf",&bdm);
	initia_p=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	vfrap=(int)round(bdm);
	
	
	fscanf (bu,"%lf", &bdm);
	robs=bdm;
	printf("Here is the bdm in the same sense which could be distinguished by the same kind of act %lf \n", bdm);
	
	
	
	fscanf (bu,"%lf", &bdm);
	obsd=bdm;
	
	
	

	fscanf(bu,"%lf",&bdm);
	step=bdm;
	fscanf(bu,"%lf",&bdm);
	tsim=(int)round(bdm);
	
	fscanf(bu,"%lf",&bdm);
	probation=(int)round(bdm);
	
	
	//tsim=100000;
	ivf=0.005;
	
	// Definition of global constants ends       !!!!!!!!!!!!!!!!
	
					
	// Correlation function        !!!!!!!!!!!!!!!!
	
	
	
	
	// Correlation function declaration end         !!!!!!!!!!!!!!!!			
	
	
	
	
  	// Definition of simulation constants starts        !!!!!!!!!!!!!!!!
	
	
	
	
	//fvf= 0.02 ;
	fscanf(bu,"%lf",&bdm);
	fvf=bdm;
	
	///printf("hishshshisisisias %lf \n",fvf);
	
	int war;
	
	///fscanf(bu,"%d",&war);
	fscanf(bu,"%lf",&bdm);
	war=(int)ceil(bdm);
	war=war+1;
	

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
	
	
	
	
	for (i=1; i<war; i++)
	{
		sco[i]=0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=1; i<war; i++)
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
	
	
	//sco[0]=sco[1];
	//printf("Here it is me who is going to reduce %d  %d\n", sco[0], sco[1]);
	
	
	for (i=1; i<war; i++)
	{
		ar[i]=1.0;
	}
	fscanf(bu,"%lf",&bdm);
	
	
	for (i=1; i<war; i++)
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
	
	//printf("Who am i %lf  %lf \n", ar[0], ar[1]);
	ar[0]=3.0;
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		sar[i]=1.0;
	}
	fscanf(bu,"%lf",&bdm);
	
	////printf("bdmbdm bmd bmd bmdkb fg %lf ", bdm);
	for (i=1; i<war; i++)
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
	for (i=1; i<war-1; i++)
	{
		red=red+per[i];
	}
	
	per[war-1]=1.0-red;
	
	
	
	
	
	fscanf(bu,"%lf",&bdm);
	tphy1=(int) bdm;
	
	fscanf(bu,"%lf",&bdm);
	tphy2=(int) bdm;
	
	fscanf(bu,"%lf",&bdm);
	ftaverage=(int) bdm;
	
	
	fscanf(bu,"%lf",&bdm);
	time_interval=(int) bdm;
	
	printf("Here is the problem %d  %d  %d  %d\n", tphy1, tphy2, ftaverage, time_interval);
	
	
	
	
	
	
	
	
	
	
	//fscanf (bu,"%lf", &bdm);

	
	//printf("Here is the bdm in the same sense which could be distinguished by the same kind of act %lf  %lf %lf\n", bdm, sar[0], sar[1]);
	
	
	
	
	//Last fraction        !!!!!!!

	///printf("Here is the one part which could not be done in the same one %lf \n", robs);
	
	
	// Shape parameters 0 and < 1.0 ; Oblate and 1 and > 1.0 ; Prolate          !!!!!!!!!!!!!! 
	
	SHA = 1 ;
	
	
	
	
	// Definition of simulation constants ends        !!!!!!!!!!!
	
	
	
	// Definition of compression parameter starts          !!!!!!!!!!
	
	veec=200;
	vees=0.005;    // Compression steps!
	vel=(int)(100.0*(fvf/0.5));
	//vel=10;
	
	
	// Definition of compression parameter ends            !!!!!!!!!!
	
	

	// Declaration of core bond variable starts          !!!!!!!!!!
	
	double epsi;
	
	// Declaration of core bond variable ends         !!!!!!!!!!
	
	
	
	// Evaluation of bond variable starts           !!!!!!!!!
	
	epsi=EPSILON;
	
	bn=40;
	bnll=40;
	sst=1;
	
	
	//double epsi1, epsi2, epsi3, omega2, beta1, beta2, beta3, betai;
	
	// Declaration of core bond variable ends!
	
	
	// Evaluation of bond variable starts!
	

	
	
	

	fclose(bu);
	remove("Base_data.txt");



	
	
	
	
	
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
	
	for (i=0; i<war; i++)
	{
		if (sco[i]==0)
		{
			
			d[i]=sar[i]*pow(ar[i],(2.0/3.0));
			if(i==0)
			{
				d[i]=2.0*robs*ar[0]; 
				
			} 
			r[i]=d[i]/2.0;
			rs[i]=r[i]*r[i];
			rsm[i]=rs[i]/(ar[i]*ar[i]);
			rm[i]=r[i]/ar[i];
			
			if(i==0)
			{
				rsm[i]=robs*robs;
				rm[i]=robs;
			} 
			
			
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
			if (i==0)
			{
			
				ds[i][j]=robs;
				dss[i][j]=robs;
			
			
			}
			else
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
			}
			if (j==0)
			{
				
				ds[i][j]=ds[i][j]+robs;
				dss[i][j]=dss[i][j]+robs;
			
			}
			else
			{
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
			}
			
			dss[i][j]=dss[i][j]+epsi;
			
			dss[i][j]=dss[i][j]*dss[i][j];
			ds[i][j]=ds[i][j]*ds[i][j];
			
			
		}
		
		
	}
	
	
	ble=1.0;
	for (i=1; i<war; i++)
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
	
	
	bleob=0.0;
	for (i=1; i<war; i++)
	{
		if (bleob < sqrt(dss[0][i])) bleob=sqrt(dss[0][i]);
		
		/////////printf("humgamamammmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm %lf \n", sqrt(dss[0][i]));
	}
	
	
	
	
	
	
	
	
	
	
	

	nt[0]=0;
	it=0;
	i=0;
	
	while(it<500)
	{
		i=i+1;
		nt[0]=i*i;
		fe=pow(nt[0]*piee*robs*robs/obsd, 0.5);
		it=(int)(((fe*fe*fe)-piee*fe*robs*robs*(float)nt[0])*fvf*6.0/piee);
		
	}	
	
	f[0]=fe;
	f[1]=fe;
	f[2]=fe;
	
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
	
	if (it>600)
	{
		tri=0;
		while(tri==0)
		{
			if(it<600 || block[1]<4)
			{
				tri=1;
			}
			else
			{
				block[1]=block[1]-1;
				f[1]=bl[1]*block[1];
				it=(int)(((f[0]*f[1]*f[2])-piee*f[1]*robs*robs*(float)nt[0])*fvf*6.0/piee);
			}
		}
	}
	
	for (i=1; i<war  ; i++)
	{
		nt[i]=(int)round((float)it*per[i]);
		//printf("%d  \n", nt[i]);
	}
	N=0;
	for (i=0; i<war; i++) N=N+nt[i];
	
	
	
	
	
	
	
	
	
	
	
	blob[0]=bleob;
	blob[1]=bleob;
	blob[2]=bleob;
	
	blockob[0]=(int)floor(f[0]/blob[0]);
	blob[0]=f[0]/(float)blockob[0];
	
	
	blockob[1]=(int)floor(f[1]/blob[1]);
	blob[1]=f[1]/(float)blockob[1];
	
	
	blockob[2]=(int)floor(f[2]/blob[2]);
	blob[2]=f[2]/(float)blockob[2];
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Evaluation of bond variable ends!
	
	printf("Number of particles-    %d \n", N-nt[0]);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of cylindrical obstacle grids-    %d \n\n\n\n", nt[0]);
	printf("Radius of obstacle-    %lf\n", robs);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Density of obstacle-    %lf\n", obsd);
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
	else printf("Single particle diffusion applied.\n");
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
	
		printf("Aspect ratio for particle type %d is given by-   %lf\n", i, ar[i]);
		
		printf("Size factor for particle type %d is given by-   %lf\n", i, sar[i]);
	
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of particles for each type (aspect ratios) are given by \n");
	for (i=1; i<war; i++)
	{
		printf("Number of particles for particle type %d is given by-   %d\n", i, nt[i]);
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	
	
	
	printf("                                ******:::::::   Calculated parameters are given by   :::::::******\n");
	printf("\n");
	printf("Starting edges are given by-  %lf  %lf   %lf\n", f[0], f[1], f[2]);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of shell for each edge is given by-  %d  %d  %d\n", block[0], block[1], block[2]);
	printf("\n");
	printf("Shell length is given by-   %lf  %lf  %lf \n\n\n", bl[0], bl[1], bl[2]);
	printf("\n");
	
	
	printf("Number of shell at obstacle grid, for each edge is given by-  %d  %d  %d\n", blockob[0], blockob[1], blockob[2]);
	
	
	printf("Shell length for obstacle grid is given by-   %lf  %lf  %lf \n\n\n", blob[0], blob[1], blob[2]);
	
	pbc=ble+0.5;

	rho=obsd/(piee*robs*robs);
	sss=sqrt(1.0/rho)-2.0*robs;
	qqq=sqrt(2.0/rho)-2.0*robs;
	printf("Number fraction of the obstacles are given as- %lf, surface to surface closest distance- %lf, surface to surface diagonal distance- %lf \n\n\n", rho, sss, qqq);
	
	
	bdm=0.0;
	for (i=1; i<war; i++)
	{
		if (ar[0]>1)
		{
			if (rm[i]>bdm) bdm=rm[i];
		}
		else
		{
			if (r[i]>bdm) bdm=r[i];
		}
	
	}
	///printf("%lf %lf \n", 2.5*bdm, qqq);
	if (qqq<(2.3*bdm)) 
	{
		printf("Your given voids is so small to fit a particle so the code will not run. Please put higher value of obstacle density. Or increase the radius of obstalce.\n\n\n\n");
		printf("Thanks !!!!!!!!!\n\n\n");
		exit(0);
	}
	else 
	{
		printf("You have given valid obstacle density the code will run. Thanks!!!!!!!!\n\n\n\n");
	}
	
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
	
	int wpi;
	
	
	
	
	
	double** li=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		li[i]=malloc(9*sizeof(double));
	}
	
	
	double** fau_li=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		fau_li[i]=malloc(9*sizeof(double));
	}
	
	double** el=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		el[i]=malloc(9*sizeof(double));
	}
	
	double** ami=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		ami[i]=malloc(3*sizeof(double));
	}
	
	
	double** current_ami=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		current_ami[i]=malloc(3*sizeof(double));
	}
	
	
	double** fau_ami=malloc((N+1)*sizeof(double*));
	
	for (i=0;i<N+1;++i)
	{
		fau_ami[i]=malloc(3*sizeof(double));
	}
	
	int* pid=malloc((N+1)*sizeof(int));
	
	
	
	double** ldcxyz=malloc(N*sizeof(double*));
	double** ldcxy=malloc(N*sizeof(double*));
	double* ldcz=malloc(N*sizeof(double));
	double** lom=malloc(N*sizeof(double*));
	
	
	for (i=0; i<N; i++)
	{
		ldcxyz[i]=malloc(3*sizeof(double));
		ldcxy[i]=malloc(2*sizeof(double));
		lom[i]=malloc(3*sizeof(double));
	
	}
	
	
	
	
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
	
	
	
	int**** dri=malloc(block[0]*sizeof(int***));
	int*** wdr=malloc(block[0]*sizeof(int**));
		
	for (j=0; j<block[0]; j++)
	{
		dri[j]=malloc(block[1]*sizeof(int**));
		wdr[j]=malloc(block[1]*sizeof(int*));
		for (i=0; i<block[1]; i++)
		{
			dri[j][i]=malloc(block[2]*sizeof(int*));
			wdr[j][i]=malloc(block[2]*sizeof(int));
	
			for (k=0; k<block[2]; k++)
			{
				dri[j][i][k]=malloc(80*sizeof(int));
			}
		}	
	}
	
	int**** driob=malloc(blockob[0]*sizeof(int***));
	int*** wdrob=malloc(blockob[0]*sizeof(int**));
		
	for (j=0; j<blockob[0]; j++)
	{
		driob[j]=malloc(blockob[1]*sizeof(int**));
		wdrob[j]=malloc(blockob[1]*sizeof(int*));
		for (i=0; i<blockob[1]; i++)
		{
			driob[j][i]=malloc(blockob[2]*sizeof(int*));
			wdrob[j][i]=malloc(blockob[2]*sizeof(int));
	
			for (k=0; k<blockob[2]; k++)
			{
				driob[j][i][k]=malloc(20*sizeof(int));
			}
		}	
	}
	
	
	// Declaration of compresser variables starts!
	
	int vee;
	
	// Declaration of compresser variables ends!
	
	// Declaration of core variables ends!
	
	
	
	
	
	
	// Declaration of bonding variable starts!
	
	
	
	
	

	int  tnkx[27][3];
	
		
		
	
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
			nrdi[i][j][1]=0;
		
		}			
		
	}
	
	for (i=0; i<block[0]; i++)
	{
		
		for (j=0; j<block[1]; j++)
		{
			for (k=0; k<block[2]; k++)
			{
				wdr[i][j][k]=0;
				for (l=0; l<80; l++)
				{
					dri[i][j][k][l]=6000;
				}
			}	
		}
	}
	
	
	for (i=0; i<blockob[0]; i++)
	{
		
		for (j=0; j<blockob[1]; j++)
		{
			for (k=0; k<blockob[2]; k++)
			{
				wdrob[i][j][k]=0;
				for (l=0; l<20; l++)
				{
					driob[i][j][k][l]=6000;
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
			jtt=jtt+nt[j];
			if (i>=jtt)    jt=jt+1;		
		}
		pid[i]=jt;
		
		
		//printf("%d  %d\n", i, pid[i]);
	}
	pid[N]=pidres;
	
	
	
	
	double* rdtheta=malloc(war*sizeof(double));
	double* grxy=malloc(war*sizeof(double));
	double* grz=malloc(war*sizeof(double));
	double* grth=malloc(war*sizeof(double));
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Dynamical parameters are given as- \n");
	
	/// Pseudo gr for obstacle
	
		
	grz[0]=1.0;
	grxy[0]=1.0;
	grth[0]=1.0;
	
	for (jt=1; jt<war; jt++)
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
	
	
	

	
	// For the safety if few particles are out of consideration to fill the area for a particular density !!!!!!!!!!!
	
	
	
	it=(int)sqrt(nt[0]);
	s[0]=fe/(float)it;
	s[1]=fe/(float)it;
	
	//s[0]=1.5;
	//s[1]=1.5;
	//printf("%d \n", it);
	
	ep=s[1]/2.0;
	jt=0;
	for (i=0; i<it; i++)
	{
		//printf("%d\n",i);
		
		ec=s[0]/2.0;
		
		for(j=0; j<it; j++)
		{
			li[jt][6]=0.0;
			li[jt][7]=1.0;
			li[jt][8]=0.0;
			ami[jt][0]=ec;
			ami[jt][1]=0.5;
			ami[jt][2]=ep;
			ec=ec+s[0];
			
			xi=ami[jt][0];
			yi=ami[jt][1];
			zi=ami[jt][2];
	    		xk=(int)floor(xi/blob[0]);
	    		yk=(int)floor(yi/blob[1]);
	    		zk=(int)floor(zi/blob[2]);
	    	
			
			for (yk=0; yk<blockob[1]; yk++)
			{
				dum=wdrob[xk][yk][zk];
				driob[xk][yk][zk][dum]=jt;
				wdrob[xk][yk][zk]=wdrob[xk][yk][zk]+1;
			}
			
			jt=jt+1;
		}
		ep=ep+s[1];
	}
	
	// Parallel Cylindrical obstacle in 3d A carbon nanotube concept !!!!!!!!!!!!!!
	
	
	

	
	//printf("raam\n");
	
	
	
	
	
     













	//Preparing for simulation!
     
     
    
     
     
     
     
     
     
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
	nlies=2*nt[1];
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
	
	sitts=0;
	
	FILE* vu;
		
	int trdf=tphy1-1;
	

	char filename[255] = {0};

	
	
	// for tphy1=125 to tphy2=200 you will only need file 1.txt; not requuired previous file 0.txt containing the data for tphy1=0 to tphy2=100; 
	
	sitt=0;
	sitts=1;
	for (i=0; i<(int)floor((float)(0+1)/100.0); i++)
	{
		sitts=sitts+100;
	}

	/// above section ends here!
	
	sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));
	
	printf("%d %d\n", sitts, (int)floor((float)(tphy1+1)/100.0));
	vu=fopen(filename,"r");

	
	printf("dhdhdhddhdhddhdd  here is the watch do g %lf \n", ds[1][1]);

	
	
	
	
	double *** initial_ami=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		initial_ami[i]=malloc(N*sizeof(double*));
		for (j=0; j<N; j++)
		{
			initial_ami[i][j]=malloc(3*sizeof(double));
		
		}
	}
	
	
	double *** inic_ami=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		inic_ami[i]=malloc(N*sizeof(double*));
		for (j=0; j<N; j++)
		{
			inic_ami[i][j]=malloc(3*sizeof(double));
		
		}
	}
	
	
	double *** initial_li=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		initial_li[i]=malloc(N*sizeof(double*));
		for (j=0; j<N; j++)
		{
			initial_li[i][j]=malloc(9*sizeof(double));
		
		}
	}
	
	int* individual_clock=malloc(ftaverage*sizeof(int));
	
	
	for (i=0; i< ftaverage; i++)
	{
		individual_clock[i]=0;
	}
	
	
	
	
	target_time=tphy1;
	
	
	starter_flag=0;
	
	jt=(tphy2-tphy1);
	double** corre_orr=malloc(jt*sizeof(double*));
	for (i=0; i<jt; i++)
	{
		corre_orr[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<jt; i++)
	{
		for(j=0; j<2; j++)
		{
			corre_orr[i][j]=0.0;
		}
	}
	
	
	
	jt=(tphy2-tphy1);
	double** corre_trans=malloc(jt*sizeof(double*));
	for (i=0; i<jt; i++)
	{
		corre_trans[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<jt; i++)
	{
		for(j=0; j<2; j++)
		{
			corre_trans[i][j]=0.0;
		}
	}
	
	
	
	jt=(tphy2-tphy1);
	double** msd_trans=malloc(jt*sizeof(double*));
	for (i=0; i<jt; i++)
	{
		msd_trans[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<jt; i++)
	{
		for(j=0; j<2; j++)
		{
			msd_trans[i][j]=0.0;
		}
	}
	
	
	
	double** msd_or=malloc(jt*sizeof(double*));
	for (i=0; i<jt; i++)
	{
		msd_or[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<jt; i++)
	{
		for(j=0; j<2; j++)
		{
			msd_or[i][j]=0.0;
		}
	}
	
	
	
	
	
	
	
	
	radial_width=0.01; radial_range=5.0; 
	n_sec=(int)(radial_range/radial_width);

	double***  time_ave_histogram=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		time_ave_histogram[i]=malloc(n_sec*sizeof(double*));
		for (j=0; j<n_sec; j++)
		{
			time_ave_histogram[i][j]=malloc(2*sizeof(double));
		
		}
	
	}
	
	
	
	double***  time_ave_histobs=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		time_ave_histobs[i]=malloc(n_sec*sizeof(double*));
		for (j=0; j<n_sec; j++)
		{
			time_ave_histobs[i][j]=malloc(2*sizeof(double));
		
		}
	
	}
	
	
	
	
	
	
	double*** lor=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		lor[i]=malloc(N*sizeof(double*));
		for (j=0; j<N; j++)
		{
			lor[i][j]=malloc(3*sizeof(double));
		}
	}
	
	for (i=0; i<ftaverage; i++)
	{
		for(j=0; j<N; j++)
		{
			for (k=0; k<3; k++)
			{
				lor[i][j][k]=0.0;
			}
		}
	}
	
	
	
	
	
	
	for (i=0; i<ftaverage; i++)
	{
		for (j=0; j<n_sec; j++)
		{
			time_ave_histogram[i][j][0]=0.0;
			time_ave_histogram[i][j][1]=0.0;
			time_ave_histobs[i][j][0]=0.0;
			time_ave_histobs[i][j][1]=0.0;
		}
	}
	
	
	
	
	k_width=1.0/fe; k_range=40; 
	n_k_sec=(int)(k_range/k_width);

	double***  ta_structure_factor=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		ta_structure_factor[i]=malloc(n_k_sec*sizeof(double*));
		for (j=0; j<n_k_sec; j++)
		{
			ta_structure_factor[i][j]=malloc(2*sizeof(double));
		
		}
	
	}
	
	
	
	
	angle_width=0.1; angle_range=180; 
	n_a_sec=(int)(angle_range/angle_width);

	double***  ta_dir_theta=malloc(ftaverage*sizeof(double**));
	for (i=0; i<ftaverage; i++)
	{
		 ta_dir_theta[i]=malloc(n_a_sec*sizeof(double*));
		for (j=0; j<n_a_sec; j++)
		{
			 ta_dir_theta[i][j]=malloc(2*sizeof(double));
		
		}
	
	}
	
	
	
	
	
	
	for (i=0; i<ftaverage; i++)
	{
		for (j=0; j<n_k_sec; j++)
		{
			ta_structure_factor[i][j][0]=0.0;
			ta_structure_factor[i][j][1]=0.0;
		
		}
		
		for (j=0; j<n_a_sec; j++)
		{
			
			ta_dir_theta[i][j][0]=0.0;
			ta_dir_theta[i][j][1]=0.0;
		}
		
		
	}
	
	
	
	
	
	
	
	
	
	jt=(tphy2-sitts)+10;
	double** s_angle_time=malloc(jt*sizeof(double*));
	
	for (i=0; i<jt; i++)
	{
		s_angle_time[i]=malloc(2*sizeof(double));
	
	}
	
	for (i=0; i<jt; i++)
	{
		s_angle_time[i][0]=0.0;
		s_angle_time[i][1]=0.0;
	}
	
	double obs_dire[3]={0.0, 1.0, 0.0};
	
	
	for (i=0; i<N; i++)
	{
		current_ami[i][0]=0.0;
		current_ami[i][1]=0.0;
		current_ami[i][2]=0.0;
	
	}
	
	
	
	
	
	
	pbc=0.5*f[0];
	
	
	
	
	printf("%d  %d\n", tphy2, tphy1);
	
	
	
	
	
	rul=0;
	for(sitt=sitts; sitt<tphy2; sitt++)
	{ 
	
	// Initialization of bond variables starts!
	
		//printf("Here is your simulation %d  %d\n", sitt, sitts);
	
		
		
		
		fscanf(vu,"%lf", &f[0]);
		fscanf(vu,"%lf", &f[1]);
		fscanf(vu,"%lf", &f[2]);
		
		
		//printf("Go and yourself %lf\n", f[0]);
		
		for (i=0; i<N; i++)
		{
			
			
			//printf("Rise_above_hate!!!!!!!!%d %d\n", i, N);
			for(l=0;l<3;l++)
			{
				fscanf(vu,"%lf", &ami[i][l]);
			}
			
			
			for (j=0; j<9; j++) 
			{
				fscanf(vu,"%lf", &bdm);
			}
			
			
			for(l=0;l<9;l++)
			{
				fscanf(vu,"%lf", &li[i][l]);
			}
				//printf("Rise_above_hate!!!!!!!!%d  %d\n", i, N);
		
		}
			
		
		if(sitt==sitts)
		{
		
			for (i=0; i<N; i++)
			{
				fau_ami[i][0]=ami[i][0];
				fau_ami[i][1]=ami[i][1];
				fau_ami[i][2]=ami[i][2];
				
				
				for (j=0; j<9; j++)
				{
					fau_li[i][j]=li[i][j];
				}
			}
		
		}
		
		
		for (i=0; i<N; i++)
		{
			dx=ami[i][0]-fau_ami[i][0];
			dy=ami[i][1]-fau_ami[i][1];
			dz=ami[i][2]-fau_ami[i][2];
			
			if (dx < -pbc) dx=dx+f[0];
			else if (dx > pbc) dx=dx-f[0];
			
			if (dy < -pbc) dy=dy+f[1];
			else if (dy > pbc) dy=dy-f[1];
			
			if (dz < -pbc) dz=dz+f[2];
			else if (dz > pbc) dz=dz-f[2];
			
			
			current_ami[i][0]=current_ami[i][0]+dx;
			current_ami[i][1]=current_ami[i][1]+dy;
			current_ami[i][2]=current_ami[i][2]+dz;
		}
		
		
		
		
		
		
		
		
			
		//printf("Budda was born to decore the world!!!!!!\n");	
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
		
			
		//printf("Budda was born to decore the world!!!!!!\n");	
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
		
		printf("Physical time outer %d\n",(sitt+1));
			
		//printf("Budda was born to decore the world!!!!!!\n");	
		
		
					
			
			
		for (jt=1; jt<war; jt++)
		{
			for (i=0; i<3; i++)
			{
				for (j=0; j<3; j++)
				{
					rtm[i][j]=0.0;
					for (k=0; k<N; k++)
					{
						rpid1=pid[k];
						if (rpid1==jt)
						rtm[i][j]=rtm[i][j]+(li[k][6+i]*li[k][6+j]);
						//printf("raaam %lf \n",li[jt][k][6+j]);
					}
					rtm[i][j]=(3.0*rtm[i][j])/(2.0*(float)nt[jt]);
					if (i==j) rtm[i][j]=rtm[i][j]-(1.0/2.0);
					//printf("raaaddddddddddddddddddm %lf \n",rtm[i][j]);
				}
			
			}
			
			for (i=0; i<3; i++)
			{
				for (j=0; j<3; j++)
				{
					abrtm[i][j]=0.0;
					for (k=0; k<3; k++)
					{
						abrtm[i][j]=abrtm[i][j]+(rtm[i][k]*rtm[k][j]);
					}
				}
			
			}
			
			
			
			double p1, eig1, eig2, eig3, q, p2, B, rr, p; 
			
			p1 = rtm[0][1]*rtm[0][1] + rtm[0][2]*rtm[0][2] + rtm[1][2]*rtm[1][2];
			if (p1 == 0)
			{ 
				   printf("A is diagonal.");
				   eig1 = rtm[0][0];
				   eig2 = rtm[1][1];
				   eig3 = rtm[2][2];
			}
			else
			{
				q = (rtm[0][0]+ rtm[1][1]+ rtm[2][2])/3.0;               
				p2 = ( rtm[0][0] - q)*( rtm[0][0] - q) + ( rtm[1][1] - q)*( rtm[1][1] - q) + ( rtm[2][2] - q)*( rtm[2][2] - q) + 2.0 * p1;
				p = sqrt(p2 / 6.0);

				for (i=0; i<3; i++)
				{
					for (j=0; j<3; j++)
					{
						
						if (i==j)
						{
							abrtm[i][j]=(1.0/p)*(rtm[i][j]-q);
						}
						
						else
						{
							abrtm[i][j]=(1.0/p)*(rtm[i][j]);
						
						}
						
					}
				
				}
				
				rr=abrtm[0][0]*(abrtm[1][1]*abrtm[2][2]-abrtm[2][1]*abrtm[1][2])-abrtm[0][1]*(abrtm[1][0]*abrtm[2][2]-abrtm[2][0]*abrtm[1][2]);
				rr=rr+abrtm[0][2]*(abrtm[1][0]*abrtm[2][1]-abrtm[2][0]*abrtm[1][1]);
				rr=rr/2.0;
			   
				
				if (rr <= -1)
				{
					phi = piee / 3.0;
				}
				else if (rr >= 1)
				{
					phi = 0.0;
				}
				else
				{	
					phi = acos(rr) / 3.0;
				}

				
				eig1 = q + 2.0 * p * cos(phi);
				eig3 = q + 2.0 * p * cos(phi + (2.0*piee/3.0));
				eig2 = 3 * q - eig1 - eig3   ;  
			
			}
			
			
			
			
			
			
			a=1.0;
			b=0.0;
			c=(1.0/2.0)*(abrtm[0][0]+abrtm[1][1]+abrtm[2][2]);
			dd=rtm[0][0]*((rtm[1][1]*rtm[2][2])-(rtm[2][1]*rtm[1][2]));
			dd=dd+(rtm[0][1]*((rtm[2][0]*rtm[1][2])-(rtm[1][0]*rtm[2][2])));
			dd=dd+(rtm[0][2]*((rtm[1][0]*rtm[2][1])-(rtm[2][0]*rtm[1][1])));
			//if(rtm[0][0]+rtm[1][1]+rtm[2][2]>0.01) printf("raaam%lf\n",rtm[0][0]+rtm[1][1]+rtm[2][2]);
			
			// Determinants measurment with non-generalized scheme!
			
			//printf("%lf\n",dd);
			if (fabs(dd)<0.0001) dd=fabs(dd);
			sol[0]=pow(((4.0*dd)/a),(1.0/3.0));
			sol[1]=((-1)*sol[0])/2.0;
			sol[2]=sol[1];
			qua=pow(((4.0/3.0)*c),0.5);
			//printf(" Total eigen values for particle type %d is given by %lf %lf %lf %lf \n", jt, rtm[0][0], rtm[1][1], rtm[2][2]+rtm[0][0]+rtm[1][1], qua);
			
			printf(" Total eigen values for particle type %d is given by %lf \n", jt, eig1);
		
		
		
		
			for (i=0; i<3; i++)
			{
				for (j=0; j<3; j++)
				{
					
					if (i==j)
					{
						lt[i][j]=rtm[i][j]-eig1;
					}
					else
					{
						lt[i][j]=rtm[i][j];
					}
				
				}
			
			}
			
			
		
		
			lt[0][3]=0;
			lt[1][3]=0;
			lt[2][3]=0;
			cont=1;
			
			
			
			for (it=0; it<2; it++) 
			{
			    	dummy=0;
			    	if  (lt[it][it]==0.0)
				{
					for (jtt=0; jtt<3; jtt++) 
					{
				    		if (dummy == 0)	
						{
							if (lt[jtt][it] != 0.0) 
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
			    	for (jtt=cont; jtt<3; jtt++)
				{
					dumm=lt[jtt][it];
					for (mt=0; mt<4; mt++)
					{
				    		lt[jtt][mt]=lt[jtt][mt]-((lt[it][mt]/lt[it][it])*dumm);
					}
				}
			    	cont=cont+1;
			}
			
			
			for (i=0; i<3; i++)
			{
				if (sol[i]!=1.0) sol[i]=0.0;
			}
			
			if (nt[jt]>0)
			{
				for (i=0; i<3; i++)
				{
					trim=0;
					for (j=0; j<3; j++)
					{
						
						
							if (fabs(lt[i][j])>EPSILON)
								{
									trim=1; //printf("   raaam");
								}
								
								//printf("%lf     ", lt[i][j]);
							
					
					}
					
					//printf("\n ");
					if(trim==0) {sol[i]=1.0; trim=1; break;}
					else { trim=2;}
				}
				
				
				//NOW CHECK EXISTENCE OF SOLUTION OF EQUATION
				
				for (i=0; i<3; i++)
				{
					if (sol[i]!=1.0) sol[i]=0.0;
				}
			 	kt=3;
				for (it=trim; it>=0; it--)
				{    	
					sol[it]=lt[it][kt];
				    	for (jtt=0; jtt<3; jtt++)
					{
						if (it==jtt) continue;
						sol[it]=sol[it]-(lt[it][jtt]*sol[jtt]);
					}
				   	sol[it]=sol[it]/lt[it][it];
				}
				mr=sqrt(sol[0]*sol[0]+sol[1]*sol[1]+sol[2]*sol[2]);
				sol[0]=sol[0]/mr;
				sol[1]=sol[1]/mr;
				sol[2]=sol[2]/mr;
				
				printf("Eigen vector for the largest eigen values is given by, %lf   %lf    %lf\n", sol[0], sol[1], sol[2]);
		
		
			}
		
			s_angle_time[sitt][0]=eig1;
			
			bdm=obs_dire[0]*sol[0]+obs_dire[1]*sol[1]+obs_dire[2]*sol[2];
			bdm=180*acos(fabs(bdm))/piee;
			
			director_ore[0]=sol[0];
			director_ore[1]=sol[1];
			director_ore[2]=sol[2];
			s_angle_time[sitt][1]=bdm;
		
		}
			
			
			
			
			
			
			
		if(sitt>=(tphy1-1))
		{	
			
			
			
			
			
			
			
			
			
			
			
			printf("Physical time inner %d\n",(sitt+1));
			if (sitt==target_time && starter_flag<ftaverage)
			{
				for (i=0; i<N; i++)
				{
					initial_ami[starter_flag][i][0]=ami[i][0];
					initial_ami[starter_flag][i][1]=ami[i][1];
					initial_ami[starter_flag][i][2]=ami[i][2];
					
					inic_ami[starter_flag][i][0]=current_ami[i][0];
					inic_ami[starter_flag][i][1]=current_ami[i][1];
					inic_ami[starter_flag][i][2]=current_ami[i][2];
					
					for (j=0; j<9; j++)
					{
						initial_li[starter_flag][i][j]=li[i][j];
					
					}
				
				}
				
				
				//cph=0.0;
				//cpsh=0.0;
				///printf("Here is the loop %d", i);
				jt=nt[0];
				for (jdt=1; jdt<war; jdt++)
				{
					it=jt+nt[jdt];
					cpsh=-fe;
					for (i=0; i<3; i++)
					{
						
						printf("Here is the loop %d\n", i);
						cph=-fe;
						for (j=0; j<3; j++)
						{	
							cp=-fe;
							for (k=0; k<3; k++)
							{		
								for (l=jt; l<it; l++)
								{
									dx=ami[l][0];
									dy=ami[l][1];
									dz=ami[l][2];
									for(m=jt; m<it; m++)
									{
										mxi=dx-ami[m][0]-cp;
										myi=dy-ami[m][1]-cph;
										mzi=dz-ami[m][2]-cpsh;
										bdm=sqrt(mxi*mxi+myi*myi+mzi*mzi);
										if (bdm<radial_range)
										{
											n=(int)floor(bdm/radial_width);
											if(l!=m)	time_ave_histogram[starter_flag][n][1]=time_ave_histogram[starter_flag][n][1]+1;
										}
									}
									
								}
								
								for (l=0; l<nt[0]; l++)
								{
									for(lc=0; lc<100; lc++)
									{
										dx=ami[l][0];
										dz=ami[l][2];
										dy=f[1]*(float)lc/100.0;
										for(m=jt; m<it; m++)
										{
											mxi=dx-ami[m][0]-cp;
											myi=dy-ami[m][1]-cph;
											mzi=dz-ami[m][2]-cpsh;
											bdm=sqrt(mxi*mxi+myi*myi+mzi*mzi);
											if (bdm<radial_range)
											{
												n=(int)floor(bdm/radial_width);
												time_ave_histobs[starter_flag][n][1]=time_ave_histobs[starter_flag][n][1]+1;
											}
										}
									}
								}
								
								
									
								cp=cp+fe;
							}
							cph=cph+fe;
						}
						cpsh=cpsh+fe;
					}
					jt=it;
				}
				
				
				
				
				qua=0.5*k_width;
				for (mum=0; mum<n_k_sec; mum++)
				{
					dom0=0.0;
					dom1=0.0;
					
					for (tum=0; tum<100; tum++)
					{
						ra1=drand48();
						ra2=drand48();
						
						rdphi=tpiee*ra1;
						theta=acos(1-2.0*ra2);
						
						dx=qua*sin(theta)*cos(rdphi);
						dy=qua*sin(theta)*sin(rdphi);
						dz=qua*cos(theta);
						eb=0.0;
						teb=0.0;
						jt=nt[0];
						for (jdt=1; jdt<war; jdt++)
						{
							it=jt+nt[jdt];
						
							for (l=jt; l<it; l++)
							{
								///printf("%d\n", l);
								mxi=ami[l][0];
								myi=ami[l][1];
								mzi=ami[l][2];
								
								
								bdm=dx*mxi+dy*myi+dz*mzi;
									
								eb=eb+cos(bdm);
								teb=teb+sin(bdm);
								
							
							}
							
							eb=eb*eb;
							teb=eb+teb*teb;	
							jt=it;
						}
						
						dom1=dom1+teb;
						
						
					}
					
					dom1=dom1/(100.0*nt[1]);
					
					
					ta_structure_factor[starter_flag][mum][1]=ta_structure_factor[starter_flag][mum][1]+dom1;
					qua=qua+k_width; 
				
				}
				
				
				
				
				for (i=nt[0]; i<N; i++)
				{
					
						com=li[i][6]*director_ore[0]+ li[i][7]*director_ore[1] +li[i][8]*director_ore[2];
						if (drand48()<0.5)	com=-com;
						com=acos(com);
						com=com*angle_range/piee;
						n=(int)floor(com/angle_width);
						//printf("%d    %lf   %d\n", n, com, n_a_sec);
						ta_dir_theta[starter_flag][n][0]=ta_dir_theta[starter_flag][n][0]+1;
					
					for(j=i+1; j<N; j++ )
					{
						com=li[i][6]*li[j][6]+ li[i][7]*li[j][7] + li[i][8]*li[j][8];
						
						if (drand48()<0.5)	com=-com;
						com=acos(com);
						com=com*angle_range/piee;
						n=(int)floor(com/angle_width);
						//printf("%d    %lf   %d\n", n, com, n_a_sec);
						ta_dir_theta[starter_flag][n][1]=ta_dir_theta[starter_flag][n][1]+1;
					
					}
				}
				
				
				
				
				
				
				
				
				starter_flag=starter_flag+1;
				target_time=target_time+time_interval;
				
				printf("%d\n", target_time);
			
			}
			
			
			
			for (it=0; it<starter_flag; it++)
			{
				jt=nt[0];
				for (jdt=1; jdt<war; jdt++)
				{
						
					j=jt+nt[jdt];
					com=0.0;
					printf("buddha was born here!!!!!!!!!!!!!!!!!!!\n");
					cp=10.0;
					for (i=jt; i<j; i++)	
					{
						
						qua=0.0;
						xi=inic_ami[it][i][0]-current_ami[i][0];
						yi=inic_ami[it][i][1]-current_ami[i][1];
						zi=inic_ami[it][i][2]-current_ami[i][2];
							
						for (l=0; l<1000; l++)
						{
							ra1=drand48();
							ra2=drand48();
							
							rdphi=tpiee*ra1;
							theta=acos(1.0-2.0*ra2);
							
							dx=cp*sin(theta)*cos(rdphi);
							dy=cp*sin(theta)*sin(rdphi);
							dz=cp*cos(theta);
							
							
							
							qua = qua + cos(xi*dx+yi*dy+zi*dz);
							
							
							
						}
						qua=qua/1000.0;
						
						com=com+qua;
					}
					com=com/((float)nt[jdt]);
					jt=j;
			
				}	
				
				corre_trans[individual_clock[it]][0]=corre_trans[individual_clock[it]][0]+1.0;
				corre_trans[individual_clock[it]][1]=corre_trans[individual_clock[it]][1]+com;
				
				jt=nt[0];
						
				for (jdt=1; jdt<war; jdt++)
				{
						
					j=jt+nt[jdt];
					com=0.0;
					for (i=jt; i<j; i++)	
					{
						qua=(initial_li[it][i][6]*li[i][6])+(initial_li[it][i][7]*li[i][7])+(initial_li[it][i][8]*li[i][8]);
						qua=((3.0*qua*qua)-1.0)/2.0;
						com=com+qua;
					}
					
					
					com=com/((float)nt[jdt]);
					jt=j;
			
				}
				
				corre_orr[individual_clock[it]][0]=corre_orr[individual_clock[it]][0]+1.0;
				corre_orr[individual_clock[it]][1]=corre_orr[individual_clock[it]][1]+com;
				
				
				
				jt=nt[0];
						
				for (jdt=1; jdt<war; jdt++)
				{
						
					j=jt+nt[jdt];
					com=0.0;
					for (i=jt; i<j; i++)	
					{
						
						xi=inic_ami[it][i][0]-current_ami[i][0];
						yi=inic_ami[it][i][1]-current_ami[i][1];
						zi=inic_ami[it][i][2]-current_ami[i][2];
						
						
						
						qua=(xi*xi)+(yi*yi)+(zi*zi);
						com=com+qua;
					}
					
					
					com=com/((float)nt[jdt]);
					jt=j;
			
				}
				
				msd_trans[individual_clock[it]][0]=msd_trans[individual_clock[it]][0]+1.0;
				msd_trans[individual_clock[it]][1]=msd_trans[individual_clock[it]][1]+com;
				
				com=0.0;
				
				for (i=nt[0]; i<N; i++)
				{
					dx=fau_li[i][7]*li[i][8]-fau_li[i][8]*li[i][7];	
					dy=fau_li[i][8]*li[i][6]-fau_li[i][6]*li[i][8];
					dz=fau_li[i][6]*li[i][7]-fau_li[i][7]*li[i][6];
					
					if(individual_clock[it]==0)
					{
						dx=0.0;
						dy=0.0;
						dz=0.0;
					}
					
					lor[it][i][0]=lor[it][i][0]+dx;
					lor[it][i][1]=lor[it][i][1]+dy;
					lor[it][i][2]=lor[it][i][2]+dz;
					com=com+lor[it][i][0]*lor[it][i][0]+lor[it][i][1]*lor[it][i][1]+lor[it][i][2]*lor[it][i][2];
				}
				com=com/(float)(N-nt[0]);
				
				msd_or[individual_clock[it]][0]=msd_or[individual_clock[it]][0]+1.0;
				msd_or[individual_clock[it]][1]=msd_or[individual_clock[it]][1]+com;
				
				
				
				
				individual_clock[it]=individual_clock[it]+1;
			}
			
		
	
		}
	
		
		
		
	
		
		
		
		for (i=0; i<N; i++)
		{
			fau_ami[i][0]=ami[i][0];
			fau_ami[i][1]=ami[i][1];
			fau_ami[i][2]=ami[i][2];
			
			for (j=0; j<9; j++)
			{
				fau_li[i][j]=li[i][j];
			}
		}
	
		
		
		
		
		
		
		
		
		
		
		
		
			
		//printf("Budda was born to decore the world!!!!!!\n");	
	
		if (((float)(sitt)/100.0)==((float)floor((float)(sitt)/100.0)))
		{
		
			fclose(vu);
		
			FILE* vu;
			
			
			char filename[255] = {0};

			sprintf(filename, "%d.txt", (int)floor((float)(sitt)/100.0));

			
			vu=fopen(filename,"r");
		
	
		}
		
			
		//printf("Budda was born to decore the world!!!!!!\n");	
		

	}
	
	
	
	//VTK Zone starts!
	
	
	/// For the paraview image :::::

	
	
	FILE* mpu;
	mpu=fopen("Cord_particle.vtk","w");
	
	fprintf(mpu,"# vtk DataFile Version 3.0\n");
	fprintf(mpu,"Random data to test tensors\n");
	fprintf(mpu,"ASCII\n");
	fprintf(mpu,"DATASET POLYDATA\n");
	fprintf(mpu,"POINTS %d float\n", N-nt[0]);
	for (i=nt[0]; i<N; i++)
	{
		fprintf(mpu, "%lf   %lf   %lf\n", ami[i][0], ami[i][1], ami[i][2]);
	
	}
	
	
	fprintf(mpu,"\n");
	fprintf(mpu,"POINT_DATA %d\n", N-nt[0]);
	fprintf(mpu,"\n");
	fprintf(mpu,"\n");
	fprintf(mpu,"TENSORS non-spherical_ellipsoid float\n");
	fprintf(mpu,"\n");
	
	
	for (i=nt[0]; i<N; i++)
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
		
		fprintf(mpu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
		fprintf(mpu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
		fprintf(mpu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
		fprintf(mpu,"\n");
	}
	
	
	
	fclose(mpu);

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	printf("Here all the calculations start why !!!!!!\n");
	
	
	
	FILE*  sut;
	sut=fopen("S_evolution.txt", "w");
	
	for(i=sitts; i<tphy2; i++)
	{
		fprintf(sut, "%d    %lf    %lf\n", i, s_angle_time[i][0], s_angle_time[i][1]);
		
	}
	fclose(sut);
	
	
	
	
	
	
	
	
	FILE*  cut;
	cut=fopen("Correlation_trans.txt", "w");
	FILE*  cuo;
	cuo=fopen("Correlation_orr.txt", "w");
	FILE* muo;
	muo=fopen("MSD_trns.txt", "w");

	FILE* muor;
	muor=fopen("MSD_or.txt", "w");
	
	jt=tphy2-tphy1;
	for (i=0; i<jt; i++)
	{
		if (corre_trans[i][0]>0.0)
		fprintf(cut,"%d   %lf \n", i,  corre_trans[i][1]/corre_trans[i][0]);
		
		if (corre_orr[i][0]>0.0)
		fprintf(cuo,"%d   %lf \n", i,  corre_orr[i][1]/corre_orr[i][0]);
		
		if (msd_trans[i][0]>0.0)
		fprintf(muo,"%d   %lf \n", i,  msd_trans[i][1]/msd_trans[i][0]);
		
		if (msd_or[i][0]>0.0)
		fprintf(muor,"%d   %lf \n", i,  msd_or[i][1]/msd_or[i][0]);
	
	}
	
	
	fclose(cut);
	fclose(cuo);
	fclose(muo);
	fclose(muor);
	
	
	
	FILE* duo;
	
	
	duo=fopen("Diffusivity.txt", "w");
	
	fprintf(duo, "%lf      %lf \n", msd_trans[jt-1][1]/(msd_trans[jt-1][0]*(float)(jt-1)) ,  msd_or[jt-1][1]/(msd_or[jt-1][0]*(float)(jt-1)));
	
	fclose(duo);
	
	
	
	FILE* dug;
	dug=fopen("Director_theta.txt", "w");

	ec=0.0;
	for(i=0; i<n_a_sec; i++)
	{
		for (it=0; it<starter_flag; it++)
		{
			ec=ec+ta_dir_theta[it][i][0];
		}
	}
	
	cp=angle_width*0.5;
	for (i=0; i<n_a_sec; i++)
	{
		bdm=0.0;
		for (it=0; it<starter_flag; it++)
		{
			bdm=bdm+ta_dir_theta[it][i][0];
		}
		bdm=bdm/(ec*angle_width);
		
		fprintf(dug, "%lf   %lf \n", cp, bdm);
		
		
		cp=cp+angle_width;
	}
		
	fclose(dug);
	
	
	FILE* gug;
	gug=fopen("G_theta.txt", "w");
	
	ec=0.0;
	for(i=0; i<n_a_sec; i++)
	{
		for (it=0; it<starter_flag; it++)
		{
			ec=ec+ta_dir_theta[it][i][1];
		}
	}
	
	cp=angle_width*0.5;
	for (i=0; i<n_a_sec; i++)
	{
		bdm=0.0;
		for (it=0; it<starter_flag; it++)
		{
			bdm=bdm+ta_dir_theta[it][i][1];
		}
		bdm=bdm/(ec*angle_width);
		
		fprintf(gug, "%lf   %lf \n", cp, bdm);
		
		
		cp=cp+angle_width;
	}
		
	fclose(gug);
	
	
	
	
	
	FILE* cug;	
	cug=fopen("Rdf.txt", "w");
	cp=radial_width*0.5;
	for(i=0; i<n_sec; i++)
	{
		bdm=0.0;
		for (it=0; it<starter_flag; it++)
		{
			bdm=bdm+time_ave_histogram[it][i][1];
		
		}
		bdm=bdm/(float)(nt[1]*starter_flag);
		bdm=bdm/(4.0*piee*cp*cp*radial_width*fvf*(6.0/piee));
		
		fprintf(cug, "%lf   %lf\n", cp, bdm);
		cp=cp+radial_width;
	}
	fclose(cug);
	
	
	
	
	pfvf=f[0]*f[1]*f[2];
	FILE* cugobs;	
	cugobs=fopen("Rdf_obs.txt", "w");
	cp=radial_width*0.5;
	for(i=0; i<n_sec; i++)
	{
		bdm=0.0;
		for (it=0; it<starter_flag; it++)
		{
			bdm=bdm+time_ave_histobs[it][i][1];
		
		}
		bdm=bdm/(float)(nt[1]*starter_flag);
		bdm=bdm/(4.0*piee*cp*cp*radial_width*pfvf*(6.0/piee));
		
		fprintf(cugobs, "%lf   %lf\n", cp, bdm);
		cp=cp+radial_width;
	}
	fclose(cugobs);
	
	
	
	
	
	
	
	
	
	
	
	
	FILE* cugt;
	cugt=fopen("Rdf_time.txt", "w");
	cp=radial_width*0.5;
	printf("Number of sample:    %d\n", starter_flag);
	for (it=0; it<starter_flag; it++)
	{
		cp=radial_width*0.5;
		for(i=0; i<n_sec; i++)
		{
			bdm=0.0;
			
			bdm=bdm+time_ave_histogram[it][i][1];
		
		
			bdm=bdm/(float)(nt[1]);
			bdm=bdm/(4.0*piee*cp*cp*radial_width*fvf*(6.0/piee));
			
			fprintf(cugt, "%lf   %lf\n", cp, bdm);
			cp=cp+radial_width;
		}
		fprintf(cugt,"\n");
	}	
	fclose(cugt);
	
	
	
	
	
	
	
	
	
	
	
	
	
	FILE* cuss;	
	cuss=fopen("Structure_factor_direct.txt", "w");
	cp=0.5*k_width;
	for(i=0; i<n_k_sec; i++)
	{
		dom0=0.0;
		dom1=0.0;
		for (it=0; it<starter_flag; it++)
		{
			dom0=dom0+ta_structure_factor[it][i][1];
			
		
		}
		bdm=dom0/(float)(starter_flag);
		if(cp>0.6)
		fprintf(cuss, "%lf   %lf\n", cp, bdm);
		cp=cp+k_width;
	}	
	fclose(cuss);
	
	
	
	
	
	
	
	
	
	
	
	
	FILE* cus;
	cus=fopen("Static_structure_factor.txt", "w");
	cp=radial_width*0.5;
	printf("Number of sample:    %d\n", starter_flag);

	double sec_multi=10;	
	k_width=0.01;
	k_range=40.0;
	n_k_sec=(int)(k_range/k_width);
	double*** surface=malloc((n_sec*sec_multi)*sizeof(double**));
	for(i=0; i<(n_sec*sec_multi); i++)
	{
		surface[i]=malloc(n_k_sec*sizeof(double*));
		for (j=0; j<n_k_sec; j++)
		{
			surface[i][j]=malloc(2*sizeof(double));
		}
	}
	
	for(i=0; i<(n_sec*sec_multi); i++)
	{
		
		for (j=0; j<n_k_sec; j++)
		{
			surface[i][j][0]=0.0;
			surface[i][j][1]=0.0;
		}
	}		
	
	ec=4*fvf*6.0;
	for(i=0; i<(n_sec*sec_multi); i++)
	{
		if(i<n_sec)	
		{
			bdm=0.0;
			for (it=0; it<starter_flag; it++)
			{
				bdm=bdm+time_ave_histogram[it][i][1];
			
			}
			bdm=bdm/(float)(nt[1]*starter_flag);
			bdm=bdm/(4.0*piee*cp*cp*radial_width*fvf*(6.0/piee));
		}
		else 
		{
			bdm=1.0;
		}	
			//bdm=sin(tpiee*3*cp)+sin(tpiee*4*cp);
		
		ep=cp;	
		lem=1;
		
		for (j=0; j<n_k_sec; j++)
		{
		
			///surface[i][j][0]=ec*(bdm-1)*cos(tpiee*(float)j*k_width*(float)i/(float)n_sec)*cp*cp*radial_width;
			surface[i][j][1]=ec*(bdm-1)*sin((float)(j+0.5)*k_width*ep)*ep*radial_width/((float)(j+0.5)*k_width*lem);
			
			//if (j==0) printf("%lf \n", surface[i][j][1]);
		}
		cp=cp+radial_width;
		//printf("%lf \n", cp);
		
	}
	
	cp=0.5*k_width;
	for (j=0; j<n_k_sec; j++)
	{
		//dom1=0.0;
		if (cp>0.6)
		{
		
	
			dom2=0.0;
			for(i=1; i<(n_sec*sec_multi)-1; i++)
			{
				//dom1=dom1+surface[i][j][0];
				dom2=dom2+surface[i][j][1];
			}
			jt=(n_sec*sec_multi)-1;
			dom2=1.0+0.5*(surface[0][j][1]+surface[jt][j][1])+dom2;
			//bdm=sqrt(dom2*dom2);
			fprintf(cus, "%lf   %lf\n", cp, dom2);
		}
		cp=cp+k_width;
	}
	fclose(cus);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	jdt=tphy2-tphy1;
	j=2;
	double ** vati=malloc(jdt*sizeof(double*));
	for (i=0; i<jdt; i++)
	{
		vati[i]=malloc(j*sizeof(double));
		
	}
	for (i=0; i<jdt; i++)
	{
		
		for (jt=0; jt < j; jt++)
		{
			vati[i][jt]=0.0;
		}
	}
	
	
	printf("Stretched_exponential starts !!!!!!!!\n");
	jt=tphy2-tphy1;
	for (i=0; i<jt; i++)
	{
		vati[i][0]=i;
		if (corre_trans[i][0]>0.0)
		vati[i][1]=corre_trans[i][1]/corre_trans[i][0];
	}
	
	for (i=0; i<jt; i++)
	{
		vati[i][0]=i;
		if (vati[i][1]<0.0001) {vati[i][1]=0.0001; ptr=i+1; break;}
		vati[i][1]=log(vati[i][1]);
		printf("%lf   %lf \n", vati[i][0], vati[i][1]);
	}
	coef_1=0.0; 
	coef_2=0.0;
	stretched_fitting(&coef_1, &coef_2, &coef_3, ptr, vati);
	printf("%lf    %lf    %lf\n", coef_1, coef_2, coef_3);
	printf("%lf    %lf    %lf\n", exp(coef_1), -1/coef_2, coef_3);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
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
	
	
		
	for (i=0; i<block[0]; i++)
	{
		
		for (j=0; j<block[1]; j++)
		{
		
			for (k=0; k<block[2]; k++)
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





double stretched_fitting(double* coef_1, double* coef_2, double* coef_3, int ptr, double** vati)
{

	
	int  srk;
	double tinv[500][2];
	
	
	
	
	
	double slope;
	
	
	
	
	
	
	
	bdm=(vati[ptr-1][0]-vati[0][0])/500.0;
	
	
	
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
		
		slope=(vati[i][1]- vati[i-1][1])/ (vati[i][0]- vati[i-1][0]);
		printf("slope %lf  %lf\n", bdm,  slope);
		
		ec=vati[0][0];
		
		for(j=0; j< 500; j++)
		{
			ec=ec+bdm;
			if (vati[i-1][0]< ec && ec < vati[i][0])
			{
				tinv[j][0]=ec;
				tinv[j][1]=vati[i-1][1]+(ec-vati[i-1][0])*slope;
				
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
	ptr=500;
	
	
	
	
	
	
	
	
	
	
	
	
	
	ntt=1;
	
	double* sol=malloc((ntt+1)*sizeof(double));
	//printf("raaam1\n");
	
	double* dump = malloc((ntt+2)*sizeof(double));
	//printf("raaam1\n");
	
	
	double* bump=malloc((ntt+(2*ntt)+2)*sizeof(double));
	int wbump;
	wbump=0;
	
	
	double* trump=malloc((ntt+(2*ntt)+2)*sizeof(double));
	int wtrump;
	wbump=0;
	
	
	for(i=0; i<ntt+(2*ntt)+2; i++)
	{
		trump[i]=0;
	}
	
	for (j=0; j< ptr; j++)
	{
		wbump=0;
		qua=tinv[j][1];
		bump[wbump]=qua;
		wbump=wbump+1;
		for (i=0; i< ntt; i++)
		{
			qua=qua*tinv[j][0]; 
			bump[wbump]=qua;
			wbump=wbump+1;
		}	
	    	bump[wbump]=1.0;
	    	wbump=wbump+1;
	    	
	    	qua=tinv[j][0];
	    	
	    	bump[wbump]=qua;
		wbump=wbump+1;
	 
		for (i=0;   i<((ntt*2)-1); i++)
		{
			qua=qua*tinv[j][0];
			bump[wbump]=qua;
			wbump=wbump+1;
		}
		
		for(i=0; i< ntt+(2*ntt)+2; i++)
		{
			trump[i]=trump[i]+bump[i];
		}
		
		
		
	}
	
	double** lt=malloc((ntt+1)*sizeof(double*));
	for (i=0; i<ntt+1; i++ )
	{
		lt[i]=malloc((ntt+2)*sizeof(double));
	}
	
	
	
	for ( j=0; j< (ntt+1); j++)
	{
	
		for ( i=0; i< (ntt+1); i++)
		{
			k=j+i+ntt+1;
			lt[j][i]=trump[k];   
		}
		//printf("%d\n", i);
		lt[j][i]=trump[j];
	}
	    
	cont=1;     
	     
	for (it=0; it<ntt; it++) 
	{
	    	dummy=0;
	    	if  (lt[it][it]==0.0)
		{
			for (jtt=0; jtt<ntt+1; jtt++) 
			{
		    		if (dummy == 0)	
				{
					if (lt[jtt][it] != 0.0) 
					{
			    			mt=jt;
			    			dummy=dummy+1;
			    			
			    			for(i=0; i<ntt+2; i++)
			    			{
			    				dump[i]=lt[it][i];
						
			    				lt[it][i]=lt[mt][i];
						
			    				lt[mt][i]=dump[i];
						
			    			}
			    			
					}
				}
			}
			if (dummy==0) continue; 
		}
	    	for (jtt=cont; jtt<ntt+1; jtt++)
		{
			dumm=lt[jtt][it];
			for (mt=0; mt<ntt+2; mt++)
			{
		    		lt[jtt][mt]=lt[jtt][mt]-((lt[it][mt]/lt[it][it])*dumm);
			}
		}
	    	cont=cont+1;
	}
		
		
		
	for (i=0; i<ntt+1; i++)
	{
		sol[i]=0.0;
	}
	
	kt=ntt+1;
	for (it=ntt; it>=0; it--)
	{    	sol[it]=lt[it][kt];
	    	for (jt=0; jt<ntt+1; jt++)
		{
			if (it==jt) continue;
			sol[it]=sol[it]-(lt[it][jt]*sol[jt]);
		}
	   	sol[it]=sol[it]/lt[it][it];
	   	
	   	//printf("sdfsdf %d; %lf\n",it,  sol[it]);
	}
	*coef_1=sol[0];
	*coef_2=sol[1];
	*coef_3=1.0;
	//printf("%lf    %lf    %lf\n", coef_1, coef_2, coef_3);
	

}




