// README 
// Code compilation, dependencies, and files described in Ref. 1 below

/*********************/
/***** REFERENCE *****/
/*********************/

// 1. 
// Ding, E. A., Yokokura, T. J., Wang, R., Kumar, S.
// Dissecting neurofilament tail sequence-phosphorylation-structure 
// relationships with multicomponent reconstituted protein brushes
// Proc. Natl. Acad. Sci. 121. 2024.

// Code written by: Yokokura, T. J., Duan, C., Wang, R.

// 2.
// Yokokura, T. J., Ding, E. A., Duan, C., Kumar, R., Wang, R.
// Effects of ionic strength on the morphology, scattering, 
// and mechanical response of neurofilament-derived protein brushes. 
// Biomacromolecules 25, 328-337. 2024. 

//////////////////////
/***** INCLUDES *****/
//////////////////////

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include <string.h>

#define MaxIT 5e5  //Maximum iteration steps
#define Pi 3.141592653589
#define Pi_4 Pi*4.0

/*********************/
/***** FUNCTIONS *****/
/*********************/

/** Main **/
double freeE(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB, double *eta);
double getConc(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB);
void sovDif_CR(double *in,double **g, int K, double *W, int *Ns, int sign);
void solve_PB(double *phA, double **PHA, double **PHA_T);
double And_mix(double **WA, double *wB); 

/** Reporting **/
void write_ph(double *phA, double **PHA_T, double **PHA, double *phB);
void write_W(double **WA,double *wB,double *eta);
void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym);

/** Matrix Solve **/
void thomas(int n, double **a, double *b, double *x);
void inv(int ndim, double **A, double **invA);
void solve(int ndim, double **a, double *b, double *x);
void mult(int ndim, double **A, double *B, double *x);

/*******************/
/***** GLOBALS *****/
/*******************/

/**  Box   **/
int Nx = 150;	   // Grid number
double dx;         // Grid discretization
double lx;         // Size along x-axis 
double *rx,*rx_sq; // Position

/** Protein **/
int NF_N;          // Number of chain types (= 3 for L/M/H)
double *sigma_i;   // Grafting density of each L/M/H
double **Alpha_i, **Chi_i; // Charge and FH blocks for eahc chain type
double b01, b02, b03; // Kuhn segment of L, H
double v01, v02, v03; // Segment volume of L, H
int KMAX;          // Max number of blocks in a chain type
int *K_i;          // Number of blocks in a chain type
double Nm;         // Number of Kuhn segments for longest chain type
double **PHA,*phA; // Polymer density distributions

/** Solvent **/
double vs;         // Solvent volume
double mu_s;       // Solvent chemical potential
double zs;         // Solvent fugacity
double *phB;       // Solvent density distribution

/** Ion **/
int Z_minus,Z_plus;// Abs. value valency
double *rho_elec_plus,*rho_elec_minus,*rho_elec_polym; // Cation, anion, polymer charge densities
double fugac_plus,fugac_minus;// Ion fugacities
double c0_plus,c0_minus; // Bulk salt concentrations

/** Physical **/
double *eps_prof;  // Dielectric permittivity
double eps_0,eps_r_P,eps_r_S,eps_P,eps_S; // Permittivities
double Length,Temp,Charge,kbT;  // Units

/** Mathematical **/
int Ns0;           // Number of contour steps
double ds0;        // Contour step
double A_r;        // Random amplitude 
double integ_cons; // Integration constant for freeE

/** Protein helpers **/
double v01_ref, v02_ref, v03_ref, vs_ref; // Normalized volumes
int **Nm_i, **Ns_i;// Segment number and contour number of each chain type

/** Fields **/
double *Q1,Q2;          // Polymer and solvent partition functions
double **WA, *wB, *eta; // Polymer, solvent, incompressibility fields
double *pot_elec;       // Electric potential
double *field_elec,*field_elec_sq; // Electric field

/** Energies **/
double freeEnergy,freeOld;
double freeW,freeU,freeS,free_elec,freeDiff;
double free_elec_polym,free_elec_laplace,free_elec_ion;
double freeEnergy_bulk;

/** Iteration **/
int iter;                    // Main iteration number
double wand,wopt,wcmp;       // Coefficents of iteration
int and_it;                  // Max Anderson iteration number
double Sm1,Sm2;              // Error tolerance
double wopt_PB,Sm_PB,err_PB; // Poisson-Boltzmann numerical consts
int MaxIT_PB, iter_PB;       // Poisson-Boltzmann numerical consts
int and_NrMax=1;             // Anderson history record
double ***DAs, ***WAs;       // Anderson polymer differentials
double **DBs, **WBs;         // Anderson solvent differentials
double *Cs;                  // Anderson mixing coefficients
double **WAdiff, *wBdiff;    // Anderson fields

/** Filenames **/
char phname[30],phname_elec[30],Win[30],Wname[30],itname[30];

/** Legacy & Non-specific **/
double chi,b0,v0,p0,R_p,V_p,L_Deby_P,L_Deby_S,L_Bjer_S,rho_fix;
double a11,a12,a21,a22,b1,b2,det,det1,det2;
double A1,A2;

////////////////////
/***** MACROS *****/
////////////////////

/** Memory Allocation **/
#define INIT(size) (double *)calloc(size, sizeof(double))
#define INIT_2D(NAME, size) (double **)calloc(size, sizeof(double)); for (i=0; i<size; i++) NAME[i] = INIT(size)

/** Prints **/
#define CHECK_VAL(VAR) printf("%.6e ", VAR)
#define CHECK_VAL_N(VAR) printf(#VAR": %.6e\n", VAR)

/** Matrix Collapse Reference **/
#define _IJS(i,j,s) (j*(Nm_i[k][j]+1) + s)*Nx + i


int main(int argc, char **argv)
{
	double **PHA_T;
	double *wA1,*wA2,*wA3,*wA4,*wA5;
	double *phA1,*phA2,*phA3,*phA4,*phA5;
	double ves_dx, ves_r;
	double e1,e2,e3,e4,e5,e6,e7,e8;
	double xCneck, xCltot, xCmax;
	int i,j,k,in,iseed=-3;
	int vopt;
	double vcust;
	
	int i1,i2;
	FILE *fp, *fp1;
	char seq_file[50], seq[30], pattern[20];

	time_t ts;
	iseed=time(&ts); srand(iseed);

/****************************/
/***** READ INPUT FILES *****/
/****************************/

	/***** graft1.txt *****/
	fp1=fopen("graft1.txt","r");  // Input file
	fscanf(fp1,"%d",&in);         // Initialization options
	fscanf(fp1, "%s", seq_file);  // Sequence information

	fp = fopen(seq_file, "r");
	fscanf(fp, "%s\n", seq);
	NF_N = strlen(seq);
	sigma_i = INIT(NF_N);
	Q1 = INIT(NF_N);
	K_i = (int *)calloc(NF_N, sizeof(double));
	Nm_i = (int **)calloc(NF_N, sizeof(int));
	Ns_i = (int **)calloc(NF_N, sizeof(int));
	Alpha_i = (double **)calloc(NF_N, sizeof(double));
	Chi_i = (double **)calloc(NF_N, sizeof(double));

	for (i = 0; i < NF_N-1; i++) fscanf(fp1, "%lf,", &sigma_i[i]); 
	fscanf(fp1, "%lf\n", &sigma_i[NF_N-1]); 
	
	fscanf(fp1,"%lf",&lx);                   // box size
	fscanf(fp1,"%d",&Ns0);                   // Ns0
	fscanf(fp1,"%lf,%lf",&eps_r_P,&eps_r_S); // Polymer solvent dielectric constants
	fscanf(fp1,"%lf,%lf",&c0_plus,&c0_minus);// Bulk salt concentrations [M]
	fscanf(fp1,"%d,%d",&Z_plus,&Z_minus);    // abs(Valency)
	fscanf(fp1,"%lf,%lf",&Length,&Temp);     // Units 
	fscanf(fp1,"%s",Win);		             // input file name for W field;
	fscanf(fp1,"%s",Wname);		             // output file name for W field;
	fscanf(fp1,"%lf",&A_r);                  // random amplitude
	fscanf(fp1,"%lf,%lf",&wopt,&wcmp);       // freeE and inComp mixing parameters
	fscanf(fp1,"%lf,%d,%d",&wand, &and_it, &and_NrMax);// Anderson mixing parameters
	fscanf(fp1,"%lf,%lf",&Sm1,&Sm2);         // freeDiff and inCompmax thresholds
	fscanf(fp1,"%lf,%lf,%d",&wopt_PB,&Sm_PB,&MaxIT_PB);// PB solver parameters
	fscanf(fp1, "%lf, %lf, %lf\n", &xCmax, &xCltot, &xCneck); // Optional params for initialization
	fscanf(fp1, "[%lf, %lf], [%lf, %lf], [%lf, %lf]", &b01,&v01,&b02,&v02,&b03,&v03); // Kuhn length/segment pairs for each chain type
	fclose(fp1);

	printf("b01: %.2f ", b01);
	printf("v01: %.2f\n", v01);
	printf("b02: %.2f ", b02);
	printf("v02: %.2f\n", v02);
	printf("b03: %.2f ", b03);
	printf("v03: %.2f\n", v03);

	vs = 0.02994; //nm3; volume of water molecule
	v0 = vs;      //Reference volume of integration
	v01_ref = v01/v0;
	v02_ref = v02/v0;
	v03_ref = v03/v0;
	vs_ref  = vs/v0; 
	int orig_Nm; 	

	/***** seq_file *****/
	k = 0;
	Nm = 0; //Initialize max chain length
	if (strstr(seq, "L") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]*0.360/b01);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1
		printf("L: Nm = %d, sigma_L = %.6f\n", Nm_i[k][K_i[k]-1],sigma_i[k]);
		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}
	if (strstr(seq, "M") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]*0.360/b02);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1
		printf("M: Nm = %d, sigma_L = %.6f\n", Nm_i[k][K_i[k]-1],sigma_i[k]);
		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}
	if (strstr(seq, "H") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]*0.360/b03);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1
		
		printf("%.0f N_k = %d N_AA * 0.360 [nm/AA] / %.2f [nm/K]\n", Nm, orig_Nm, b03);
		printf("H: Nm = %d, sigma_H = %.6f\n", Nm_i[k][K_i[k]-1],sigma_i[k]);
		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}
	fclose(fp);

/*****************/
/***** SETUP *****/
/*****************/

	/** Reporting **/
	snprintf(pattern, sizeof(char)*20, "_a%03.0fc%03.0f_x%03.0fs%03.0f", Alpha_i[0][0]*100, c0_plus*1000, Chi_i[0][0]*100, sigma_i[0]*100);
	snprintf(phname, sizeof(char)*30, "ph%s.dat", pattern);
	snprintf(phname_elec, sizeof(char)*30, "el%s.dat", pattern);
	snprintf(Wname, sizeof(char)*30, "W%s.dat", pattern);
	snprintf(itname, sizeof(char)*30, "it%s.dat", pattern);

	printf("%s\n", pattern);

	/** Physical constants **/
	Charge=1.6e-19;
	kbT=1.38e-23*Temp;
	c0_plus*=6.022e23*1e3*pow(Length,3.0);  // dimensionless cation concentration //
	c0_minus*=6.022e23*1e3*pow(Length,3.0);  // dimensionless anion concentration //
	eps_0=8.854e-12*kbT*Length/(Charge*Charge);  // dimensionless permittivity in vaccum //
	eps_P=eps_0*eps_r_P;  // dimensionless dielectric constant of polymer //
	eps_S=eps_0*eps_r_S;  // dimensionless dielectric constant of solvent //
	L_Bjer_S=1.0/(Pi_4*eps_0*eps_r_S);  // dimensionless Bjerrum length in solvent //
	L_Deby_P=sqrt(eps_0*eps_r_P/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  // dimensionless Debye length in polymer //
	L_Deby_S=sqrt(eps_0*eps_r_S/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  // dimensionless Debye length in solvent //

	mu_s=-1.0;
	zs=exp(mu_s);  
	fugac_plus=c0_plus; 
	fugac_minus=c0_minus;

	/** Polymer discretization **/
	ds0=1.0*Nm/Ns0; 
	CHECK_VAL_N(ds0);
	
	for (k=0;k<NF_N;k++) {
		printf("%d : [", k);
		for (j=0;j<K_i[k]-1;j++) {
			Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
			printf("%d ",Ns_i[k][j]);
		}
		Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
		printf("%d]\n",Ns_i[k][j]);
	}
	
	//// size of polymer with max radius //
	V_p=Nm*v0;
	R_p=pow(0.75/Pi*V_p,1.0/3);
	///////////////////////
	

	/** Grid discretization **/
	dx=17.0/150; //Arbitrary dx that gives reasonable results
	Nx = roundl(lx/dx);
	printf("Nx: %d\n", Nx);
	integ_cons=dx/v0; 
	
	/** Array Initialization **/
	PHA = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) PHA[i] = calloc(Nx*K_i[i], sizeof(double));
	PHA_T = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) PHA_T[i] = calloc(Nx, sizeof(double));
	WA = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WA[i] = calloc(Nx*K_i[i], sizeof(double));

	wB = INIT(Nx); phA = INIT(Nx); phB = INIT(Nx); eta = INIT(Nx);
	rx = INIT(Nx); rx_sq = INIT(Nx);
	rho_elec_minus = INIT(Nx); rho_elec_plus = INIT(Nx); rho_elec_polym = INIT(Nx); pot_elec = INIT(Nx);
	field_elec = INIT(Nx); field_elec_sq = INIT(Nx); eps_prof = INIT(Nx);

	for(i=0;i<Nx;i++) { rx[i]=i*dx; rx_sq[i]=rx[i]*rx[i]; } 

	freeEnergy_bulk=0.0;  // bulk free energy of ions //
	for(i=0;i<Nx;i++)  freeEnergy_bulk-=c0_plus+c0_minus; 
	freeEnergy_bulk*=integ_cons*v0;

	/** Electric field analytical solution **/
	// Mostly unused
	
	V_p=Nm*v0;
	R_p=pow(0.75/Pi*V_p,1.0/3);

	R_p=pow(0.75/Pi*1*Nm*v0,1.0/3);  // actual radius for each np //

	rho_fix=Alpha_i[0][0]/v0;

	a11=2.0*eps_P*(sinh(R_p/L_Deby_P)-R_p/L_Deby_P*cosh(R_p/L_Deby_P));
	a12=eps_S*(R_p/L_Deby_S+1.0)*exp(-R_p/L_Deby_S);
	a21=2.0*sinh(R_p/L_Deby_P);
	a22=exp(-R_p/L_Deby_S);

	b1=0.0;
	b2=0.5*rho_fix*R_p/c0_plus;

	det=a11*a22-a12*a21;
	det1=b1*a22-b2*a12;
	det2=a11*b2-a21*b1;

	A1=det1/det;
	A2=det2/det;

	for(i=1;i<Nx;i++)
	{
		if(rx[i]<R_p)
		{
			pot_elec[i]=A1/rx[i]*(exp(-rx[i]/L_Deby_P)-exp(rx[i]/L_Deby_P))+0.5*rho_fix/c0_plus; // spherical analytical solution
		}
		else
		{
			pot_elec[i]=A2/rx[i]*exp(-rx[i]/L_Deby_S);
		}
	}
	pot_elec[0]=pot_elec[1];
	for(i=0;i<Nx;i++) pot_elec[i]=-pot_elec[i]; // negative charge on polymer //

	printf("Analytical over\n");

/**********************************/
/***** INITIALIZATION OPTIONS *****/
/**********************************/

	/** Read from previous run **/
	if(in == -1) { 
		fp=fopen(Win,"r"); 
		for(i=0;i<Nx;i++) {
			for (k=0; k<NF_N; k++) {
				fscanf(fp, "[ ");
				for (j=0; j<K_i[k]-1; j++) {fscanf(fp, "%lf ", &WA[k][j*Nx + i]); WA[k][j*Nx + i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));}
				fscanf(fp, "%lf ] ", &WA[k][j*Nx + i]); 
				WA[k][j*Nx + i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			}
			fscanf(fp, "%lf %lf %lf\n", &wB[i], &eta[i], &pot_elec[i]);
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			eta[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			pot_elec[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		}
		fclose(fp);

		// electric field 
		i=0;i2=i+1;
		field_elec[i]=field_elec[i2];

		for(i=1;i<Nx-1;i++)
		{
			i1=i-1; i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}

		i=Nx-1;
		field_elec[i]=field_elec[i2];

		for(i=0;i<Nx;i++) field_elec_sq[i]=field_elec[i]*field_elec[i];

	}
	/** Homogeneous polymer density **/
	else if (in == 0) { // Using xCmax
		if (xCmax == 0)  xCmax = 0.80; //default
		for (i=1; i<Nx; i++){ //Respect no pen at x = 0
			for (k=0; k<NF_N; k++) {
				PHA_T[k][i] = 0;
				PHA[k][0*Nx + i] = xCmax * Nm_i[k][0]/Nm;
				for (j=1; j<K_i[k]; j++){
					PHA[k][j*Nx + i] = xCmax*((Nm_i[k][j]-Nm_i[k][j-1])/Nm);
					PHA_T[k][i] += PHA[k][j*Nx + i];
				}
				PHA_T[k][i] += PHA[k][0*Nx + i];
			}
		}
	
		for (i=0; i<Nx; i++){
			phA[i] = 0;
			for (k=0; k<NF_N; k++) phA[i] += PHA_T[k][i]; 

			phB[i]=1.0-phA[i]; 
			eps_prof[i]=phA[i]*eps_P+(1.0-phA[i])*eps_S;

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++){	
				WA[k][j*Nx + i] = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i];
			}
			
			wB[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wB[i] += Chi_i[k][j]*PHA[k][j*Nx + i];			

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++)	WA[k][j*Nx + i] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		
			pot_elec[i]=0.0;
			field_elec[i]=0.0;
			field_elec_sq[i]=0.0;
		}
		//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec, rho_elec_minus, rho_elec_polym); exit(0);
	}
	/** Assume multi-layer morphology **/
	else { //xC init, in==no. peaks; not super useful for xC = 1
		if (xCmax==0)  xCmax = 0.80; //defaults
		if (xCltot==0) xCltot = lx/2.0;
		if (xCneck==0) xCneck = xCltot/(2.0*in);
		double w;
		int c;

		w = (xCltot - (in-1)*xCneck)/in; 
		for (i=0; i<Nx; i++){
			for (k=0; k<NF_N; k++) {
				PHA_T[k][i] = 0;
				PHA[k][0*Nx + i] = 0; for (c=1; c<=in; c++) PHA[k][0*Nx + i] += xCmax*(Nm_i[k][0]/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
				for (j=1; j<K_i[k]; j++){
					PHA[k][j*Nx + i] = 0;
					for (c=1; c<=in; c++) PHA[k][j*Nx + i] += xCmax*((Nm_i[k][j]-Nm_i[k][j-1])/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
					PHA_T[k][i] += PHA[k][j*Nx + i];
				}
				PHA_T[k][i] += PHA[k][0*Nx + i];
			}
			phA[i] = 0;
			for (k=0; k<NF_N; k++) phA[i] += PHA_T[k][i]; 

    		phB[i]=1.0-phA[i]; 
			eps_prof[i]=phA[i]*eps_P+(1.0-phA[i])*eps_S;

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++){	
				WA[k][j*Nx + i] = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i];
			}
			
			wB[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wB[i] += Chi_i[k][j]*PHA[k][j*Nx + i];			

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++)	WA[k][j*Nx + i] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		
			pot_elec[i]=0.0;
			field_elec[i]=0.0;
			field_elec_sq[i]=0.0;
		}
		//printf("for out\n");
		//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec, rho_elec_minus, rho_elec_polym); exit(0);
	}
	printf("Init over\n");

	freeE(phA, PHA, PHA_T, phB, WA, wB, eta);

}

/*******************************/
/***** Reporting Functions *****/
/*******************************/

void write_ph(double *phA, double **PHA_T, double **PHA, double *phB)
{// Report polymer and solvent density distributions
	int i, j, k;

	FILE *fp=fopen(phname,"w");

	for(i=0;i<Nx;i++) 
	{
		fprintf(fp,"%10.3e ",rx[i]);
		fprintf(fp, "%10.4e ", phA[i]);
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "%10.4e ", PHA_T[k][i]);	
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", PHA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", PHA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.4e\n", phB[i]);
	}

	fclose(fp);
}

void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym)
{// Report electric potential and charge density distributions
	int i;

	FILE *fp=fopen(phname_elec,"w");

	for(i=0;i<Nx;i++)
	{	
		fprintf(fp,"%10.3e %10.4e %10.4e %10.4e %10.4e\n",rx[i],-pot_elec[i],rho_elec_plus[i],rho_elec_minus[i],rho_elec_polym[i]);
	}

	fclose(fp);
}


void write_W(double **WA,double *wB,double *eta) 
{// Report fields
	int i,j,k;

	FILE *fp=fopen(Wname,"w");

	for(i=0;i<Nx;i++) 
	{
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", WA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", WA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.5e %10.5e %10.5e\n", wB[i], eta[i], pot_elec[i]);
	}

	fclose(fp);
}

/**************************/
/***** MAIN FUNCTIONS *****/
/**************************/

double freeE(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB, double *eta)
{// Main loop, calculate equilibrium free energy
	int i,j,k,n,maxIter;
	int is;
	int MAXMAX;
	double chi_kj;
	double and_err;

	double inCompMax;
	double psum,fpsum,wpref, free_elec_polym_step,eta1,eta2,eta3,eta4,eta5,freeU_step;

	/** Local arrays **/
	double **WAnew;
	WAdiff = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAdiff[i] = calloc(Nx*K_i[i], sizeof(double));
	WAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAnew[i] = calloc(Nx*K_i[i], sizeof(double));

	double *wBnew;

	wBdiff=(double *)malloc(sizeof(double)*Nx);
	wBnew=(double *)malloc(sizeof(double)*Nx);

	DAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){ 
		DAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (k=0; k<NF_N; k++)	DAs[n][k] = (double *)calloc(K_i[k]*Nx, sizeof(double));	
	}
	WAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){
		WAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (k=0; k<NF_N; k++) WAs[n][k] = (double *)calloc(K_i[k]*Nx, sizeof(double));
	}
	DBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) DBs[n] = (double *)calloc(Nx, sizeof(double));
	WBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) WBs[n] = (double *)calloc(Nx, sizeof(double));

	maxIter=MaxIT; iter=0; freeEnergy=0.0;

	printf("freeE in\n");

/** START OF ITERATION **/
	do
	{
		iter=iter+1;

		// Get polymer density and fields
		getConc(phA, PHA, PHA_T, phB, WA, wB); 
		for(i=0;i<Nx;i++) eps_prof[i]=phA[i]*eps_P+(1-phA[i])*eps_S;

		// Solve Poisson-Boltzmann
		solve_PB(phA, PHA, PHA_T);  

		freeW=0.0; freeU=0.0; freeS=0.0; free_elec_polym=0.0;
		free_elec_laplace=0.0; free_elec_ion=0.0; inCompMax=0.0;

		for(i=0;i<Nx;i++)
		{// Integrate free energy
			psum=1.0-phA[i]-phB[i];

			fpsum=fabs(psum);
			if(fpsum>inCompMax){inCompMax=fpsum; MAXMAX = i;}

				eta1 = 0; eta2 = 0; 
				for (k=0;k<NF_N;k++){	
					eta1 += K_i[k];
					for (j=0;j<K_i[k];j++){
						eta2 += (1-phA[i])*Chi_i[k][j] - pot_elec[i]*Alpha_i[k][j] + Chi_i[k][j]*PHA[k][j*Nx + i] - WA[k][j*Nx + i];
					}	
				}
				eta[i] = 1/(eta1 + 1) * (eta2 - wB[i]);
				
				freeU_step = 0;
				for (k=0;k<NF_N;k++){
					for (j=0;j<K_i[k];j++){ 
						if (k==0) freeU_step += Chi_i[k][j]*PHA[k][j*Nx+i]/v01_ref;
						if (k==1) freeU_step += Chi_i[k][j]*PHA[k][j*Nx+i]/v02_ref;
						if (k==2) freeU_step += Chi_i[k][j]*PHA[k][j*Nx+i]/v03_ref;
					}
				}
				freeU += freeU_step * phB[i]/vs_ref;

			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) freeW -= WA[k][j*Nx+i] * PHA[k][j*Nx+i];
			freeW -= wB[i]*phB[i];

			free_elec_polym_step = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++){
					if (k==0) free_elec_polym_step += Alpha_i[k][j]*PHA[k][j*Nx+i]/v01;
					if (k==1) free_elec_polym_step += Alpha_i[k][j]*PHA[k][j*Nx+i]/v02;
					if (k==2) free_elec_polym_step += Alpha_i[k][j]*PHA[k][j*Nx+i]/v03;
			}
			free_elec_polym -= free_elec_polym_step*pot_elec[i];

			free_elec_laplace-=0.5*eps_prof[i]*field_elec_sq[i];
			free_elec_ion-=(rho_elec_plus[i]+rho_elec_minus[i]);
		}

		// Integration constants
		freeU*=integ_cons; freeW*=integ_cons; 
		free_elec_polym*=integ_cons*v0; free_elec_laplace*=integ_cons*v0; 
		free_elec_ion*=integ_cons*v0; free_elec_ion-=freeEnergy_bulk; 
		free_elec=free_elec_polym+free_elec_laplace+free_elec_ion; 

		freeS = -zs*Q2;
		for (k=0;k<NF_N;k++) freeS -= sigma_i[k]*log(Q1[k]);  

		freeOld=freeEnergy;
		freeEnergy=freeU+freeW+freeS+free_elec+freeEnergy_bulk; 
		freeDiff=fabs((freeEnergy-freeOld)/freeEnergy);

		// Iterate, perform Anderson mixing or simple mixing
		and_err = And_mix(WA, wB); 
	
		// Report 
		if(iter%10==0)printf(" %5d : %.8e, %.8e, %.8e, %.8e, %d\n", 
			iter, freeEnergy, freeDiff, inCompMax, and_err, MAXMAX);

		if(iter%100==0||iter>=maxIter)
		{

			FILE *fp=fopen("printout.dat","w");
			fprintf(fp,"%d %10.3e %10.5e %10.5e %10.5e %10.5e\n",
				iter,lx,freeEnergy,freeDiff,inCompMax,and_err);
			fclose(fp);

			fp=fopen(itname,"a");
			fprintf(fp,"%d %10.5e %10.5e %10.5e %10.5e\n",
				iter,freeEnergy,freeDiff,inCompMax, and_err);
			fclose(fp);
		}
		if(iter%100==0)
		{
			write_ph(phA, PHA_T, PHA,phB);
			write_W(WA,wB,eta);
			write_elec(pot_elec,rho_elec_plus,rho_elec_minus,rho_elec_polym);
			
		}

	}while(iter<maxIter&&(and_err>1e-3||inCompMax>Sm1||freeDiff>Sm2||iter<123));

/** END OF ITERATION **/

	// Report
	write_ph(phA, PHA_T, PHA, phB);
	write_W(WA,wB,eta);
	write_elec(pot_elec,rho_elec_plus,rho_elec_minus,rho_elec_polym);

	
	FILE *fp=fopen("printout.dat","w");
	fprintf(fp,"%d %10.3e %10.5e %10.5e %10.5e %10.5e\n",
		iter,lx,freeEnergy,freeDiff,inCompMax,and_err);
	fclose(fp);

	#if !defined(_USE_BULK)
	fp = fopen("(0)_ProgramOver.dat", "w");
	fclose(fp);
	#endif

	fp=fopen(itname,"a");
	fprintf(fp,"%d %10.5e %10.5e %10.5e %10.5e\n",
		iter,freeEnergy,freeDiff,inCompMax,and_err);
	fclose(fp);
}


double getConc(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB)
{// Get polymer and solvent densities and fields

	int i,j,k,iz,s,s1;
	double fflA;

	double *qInt;
	double **QA, **QcA;
	double **qA, **qcA, *wA;
	int s0;

/*******************/
/***** POLYMER *****/
/*******************/

	// Propagators
	QA = (double **)calloc(NF_N, sizeof(double));
	for (k=0; k<NF_N; k++) QA[k] = INIT(K_i[k] * (Ns_i[k][K_i[k]-1]+1) * Nx); 
	QcA = (double **)calloc(NF_N, sizeof(double));
	for (k=0; k<NF_N; k++) QcA[k] = INIT(K_i[k] * (Ns_i[k][K_i[k]-1]+1) * Nx);

	qInt = INIT(Nx);
	qA = (double **)calloc(Nx, sizeof(double));
	qcA = (double **)calloc(Nx, sizeof(double));
	int *Ns;
	
	// Solve MDE for each chain type
	for (k=0; k<NF_N; k++){ 
		wA = INIT(Nx*K_i[k]);
		Ns = calloc(Ns_i[k][K_i[k]-1], sizeof(int));
		for (j=0; j<K_i[k]; j++) {
			for (i=0; i<Nx; i++) wA[j*Nx + i] = WA[k][j*Nx+i];
			Ns[j] = Ns_i[k][j];
		}

		//Begin forwards prop
		for (i=0;i<Nx;i++){
			qInt[i] = 0.0;  //Graft at epsilon = 1*dx, discretized Dirac delta distribution
			qA[i] = INIT(Ns[K_i[k]-1] + 1); 
		} 
		qInt[1] = 1.0/dx;

		if (k==0) b0 = b01; //b0 global used in sovDif_CR
		if (k==1) b0 = b02; 
		if (k==2) b0 = b03; 
		sovDif_CR(qInt, qA, K_i[k], wA, Ns, 1); 

		for (j=0; j<K_i[k]; j++){
			if (j==0) s0 = 0;
			else s0 = Ns_i[k][j-1];

			for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++){
				QA[k][_IJS(i,j,s)] = qA[i][s];
			}
		}
		
		//Begin backwards prop
		for (i=0;i<Nx;i++){
			qInt[i] = 1.0; //Graft at epsilon = 1*dx, discretized Dirac delta distribution
			qcA[i] = INIT(Ns[K_i[k]-1] + 1); 
		} 
		qInt[0] = 0.0;

		sovDif_CR(qInt, qcA, K_i[k], wA, Ns, -1); 

		for (j=0; j<K_i[k]; j++){
			if (j==0) s0 = 0;
			else s0 = Ns_i[k][j-1];

			for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++){
				QcA[k][_IJS(i,j,s)] = qcA[i][s];
			}
		}

		for (i=0;i<Nx;i++) free(qcA[i]);
		for (i=0;i<Nx;i++) free(qA[i]);
		free(wA); free(Ns);
	}

	for (k=0;k<NF_N;k++){
		Q1[k]=0.0;
		for(i=0;i<Nx;i++)
		{
			Q1[k]+=QA[k][_IJS(i, (K_i[k]-1) , (Ns_i[k][K_i[k]-1]) )]; //integral of q(r, N)
		}
		if (k==0) Q1[0]*=dx/v01;	
		if (k==1) Q1[1]*=dx/v02;
		if (k==2) Q1[2]*=dx/v03;

		fflA=sigma_i[k]*ds0/Q1[k];

		for(i=0;i<Nx;i++){
			for (j=0;j<K_i[k];j++){
				PHA[k][j*Nx + i] = 0;

				if (j==0) {
					PHA[k][0*Nx+i] += 0.50*QA[k][_IJS(i,0, 0 )]*QcA[k][_IJS(i,0, 0)];
					for (s=1;s<Ns_i[k][0];s++) PHA[k][0*Nx+i] += QA[k][_IJS(i,0,s)]*QcA[k][_IJS(i,0,s)];
					PHA[k][0*Nx+i] += 0.50*QA[k][_IJS( i, 0, Ns_i[k][0] )]*QcA[k][_IJS( i, 0, Ns_i[k][0] )];
				}
				else{
					PHA[k][j*Nx+i] += 0.50*QA[k][_IJS(i,j, Ns_i[k][j-1] )]*QcA[k][_IJS(i,j, Ns_i[k][j-1] )];
					for (s=Ns_i[k][j-1]+1;s<Ns_i[k][j];s++) PHA[k][j*Nx+i] += QA[k][_IJS(i,j,s)]*QcA[k][_IJS(i,j,s)];
					PHA[k][j*Nx+i] += 0.50*QA[k][_IJS( i, j, Ns_i[k][j] )]*QcA[k][_IJS( i, j, Ns_i[k][j] )];
				}
			}
		}
		
		for (i=0; i<Nx;i++){
			PHA_T[k][i] = 0;
			for (j=0;j<K_i[k];j++){
				PHA[k][j*Nx+i] *= fflA;
				PHA_T[k][i] += PHA[k][j*Nx+i];
			}
		}
	}

	for (i=0;i<Nx;i++) {
		phA[i] = 0;
		for(k=0;k<NF_N;k++) phA[i] += PHA_T[k][i];	
	}

/*******************/
/***** SOLVENT *****/
/*******************/

  	Q2=0.0;
  	for(i=0;i<Nx;i++)
  	{
  		phB[i]=zs*exp(-wB[i]);
  		Q2+=exp(-wB[i]);
  	}
  	Q2*=dx/vs;


	free(qInt); 
	for(k=0;k<NF_N;k++){free(QA[k]); free(QcA[k]);}
	free(QA); free(QcA); free(qA); free(qcA);
}

void sovDif_CR(double *in,double **g, int K, double *W, int *Ns, int sign)
{// Solve modified diffusion equation

	int i,j,k,is;
	double *Aq,*Bq;
	double *func_in,*func_out,**Mij;

	Aq=(double *)malloc(sizeof(double)*Nx);
	Bq=(double *)malloc(sizeof(double)*Nx);
	func_in=(double *)malloc(sizeof(double)*Nx);
	func_out=(double *)malloc(sizeof(double)*Nx);

	// Tridiagonal coeff matrix
	Mij=(double **)malloc(sizeof(double)*Nx);
  	for(i=0;i<Nx;i++) Mij[i]=(double *)malloc(sizeof(double)*3);


	for(i=0;i<Nx;i++) Aq[i] = -ds0*b0*b0; //Cartesian 1D

	if(sign==1)  // forward //
	{
		is=0;
		for(i=0;i<Nx;i++) g[i][is]=in[i];  // initial value //
		
		for(is=1;is<=Ns[K-1];is++)
		{
			for (j=0; j<K; j++) if (is <= Ns[j]){ //First less than is belonging block
				for (i=0; i<Nx; i++) Bq[i] = 2*ds0*b0*b0 + 6*ds0*dx*dx*W[j*Nx+i]; 
				break;
			}

			for(i=0;i<Nx;i++)
			{
				if(i==0)
				{
					Mij[i][0]= 0.0;
					Mij[i][1]= 1.0; //Dirichlet
					Mij[i][2]= 0.0;
				}
				else if(i==Nx-1)
				{
					Mij[i][0]= 0.0;
					Mij[i][1]= 1.0; //Dirichlet
					Mij[i][2]= 0.0;
				}
				else
				{
					Mij[i][0] = Aq[i];
					Mij[i][1] = 12*dx*dx + Bq[i];
					Mij[i][2] = Aq[i];
				}
			}

			for(i=0;i<Nx;i++)
			{
				if(i==0) func_in[i]=0.0; //No penetration
				else if(i==Nx-1) func_in[i]=0.0; //Infinite 
				else 
				{
					func_in[i] = -Aq[i]*g[i-1][is-1]  //Cartesian 1D
							+ (12*dx*dx - Bq[i]) *g[i][is-1]
							-Aq[i]*g[i+1][is-1];
				}
			}

			thomas(Nx,Mij,func_in,func_out);
			
			for(i=0;i<Nx;i++) g[i][is]=func_out[i];

		}
	}
	if(sign == -1)  // backward //
	{

		is=Ns[K-1]; 
		for(i=0;i<Nx;i++) g[i][is]=in[i];  // initial value //
			
		for(is=Ns[K-1]-1;is>=0;is--)
		{
			for (j=0; j<K; j++) if (is <= Ns[j]){ //First less than is belonging block
				for (i=0; i<Nx; i++) Bq[i] = 2*ds0*b0*b0 + 6*ds0*dx*dx*W[j*Nx+i]; 
				break;
			}
			for(i=0;i<Nx;i++)
			{
				if(i==0)
				{
					Mij[i][0]= 0.0;
					Mij[i][1]= 1.0; //Dirichlet
					Mij[i][2]= 0.0;
				}
				else if(i==Nx-1)
				{
					Mij[i][0]= 0.0;
					Mij[i][1]= 1.0; //Dirichlet
					Mij[i][2]= 0.0;
				}
				else
				{
					Mij[i][0] = Aq[i];
					Mij[i][1] = 12*dx*dx + Bq[i];
					Mij[i][2] = Aq[i];
				}
			}

			for(i=0;i<Nx;i++)
			{
				if(i==0) func_in[i]=0.0; //No penetration
				else if(i==Nx-1) func_in[i]=0.0; //Infinite 
				else
				{
					func_in[i] = -Aq[i]*g[i-1][is+1]  //Cartesian 1D
							+ (12*dx*dx - Bq[i]) *g[i][is+1]
							-Aq[i]*g[i+1][is+1];			
				}
			}

			thomas(Nx,Mij,func_in,func_out);

			for(i=0;i<Nx;i++) g[i][is]=func_out[i];
		}
	}

	free(Aq);
	free(Bq);

	free(func_in);
	free(func_out);

	for(i=0;i<Nx;i++) free(Mij[i]);
	free(Mij);

}

void solve_PB(double *phA, double **PHA, double **PHA_T)
{// Solve Poisson-Boltmann Equation

	int i, j, k; 
	int i1,i2;

	double err_max,err_check;
	double *pot_elec_old;
	double *V_in,*V_out,**Aij;

	pot_elec_old = INIT(Nx);
	V_in = INIT(Nx); V_out = INIT(Nx);

	Aij = (double **)malloc(sizeof(double)*Nx);
	for(i=0;i<Nx;i++) Aij[i] = INIT(3);

	iter_PB=0;
	do  
	{
		iter_PB=iter_PB+1;

		for(i=0;i<Nx;i++) pot_elec_old[i]=pot_elec[i];

		for(i=0;i<Nx;i++)
		{
			rho_elec_plus[i] = Z_plus*fugac_plus*exp(-Z_plus*pot_elec_old[i]); 
			rho_elec_minus[i] = Z_minus*fugac_minus*exp(Z_minus*pot_elec_old[i]);
			rho_elec_polym[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++){
				if (k==0) rho_elec_polym[i] += 1/v01*Alpha_i[k][j]*PHA[k][j*Nx+i];
				if (k==1) rho_elec_polym[i] += 1/v02*Alpha_i[k][j]*PHA[k][j*Nx+i];
				if (k==2) rho_elec_polym[i] += 1/v03*Alpha_i[k][j]*PHA[k][j*Nx+i];
			}
		}

		Aij[0][1] = 1.0; Aij[0][2] = -1.0; //Neumann
		for(i=1;i<Nx-1;i++)
		{
			Aij[i][0] = 0.25*(eps_prof[i+1]-eps_prof[i-1]) - eps_prof[i]; 
   			Aij[i][1] = 2.0*eps_prof[i] + dx*dx* (Z_minus*rho_elec_minus[i] + Z_plus*rho_elec_plus[i]);
			Aij[i][2] = -0.25*(eps_prof[i+1]-eps_prof[i-1]) - eps_prof[i];
		}
   		Aij[Nx-1][1] = 1.0; //Dirichlet

   		V_in[0] = 0.0; //Symmetry
		for(i=1;i<Nx-1;i++)
		{
			V_in[i] = dx*dx * (rho_elec_plus[i]*(1+Z_plus*pot_elec_old[i]) - rho_elec_minus[i]*(1-Z_minus*pot_elec_old[i])-rho_elec_polym[i]);
		}
		V_in[Nx-1] = 0.0; //Infinite

		thomas(Nx,Aij,V_in,V_out);

		for(i=0;i<Nx;i++) pot_elec[i]=wopt_PB*V_out[i]+(1.0-wopt_PB)*pot_elec_old[i];  // Linear mixing
		
		err_max=0.0;
		for(i=0;i<Nx;i++)
		{
			err_check= (pot_elec[i]-pot_elec_old[i])*(pot_elec[i]-pot_elec_old[i]);
			if(err_check>err_max)err_max=err_check;
		}
		err_max=sqrt(err_max);

	}while(err_max>Sm_PB&&iter_PB<MaxIT_PB);
	
	//printf("iter=%d,err_max_PB=%10.5e,iter_PB=%d\n",iter,err_max,iter_PB);  // check //

	err_PB=err_max;

	///// get electric field /////
	for(i=0;i<Nx;i++)
	{
		if((i>0)&&(i<Nx-1))
		{
			i1=i-1;
			i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}
	}

	i=0;i2=i+1;
	field_elec[i]=field_elec[i2];

	i=Nx-1;i2=i-1;
	field_elec[i]=field_elec[i2];


	for(i=0;i<Nx;i++)
	{
		field_elec_sq[i]=field_elec[i]*field_elec[i];
	}


	free(pot_elec_old); free(V_out); free(V_in);
	for(i=0;i<Nx;i++) free(Aij[i]);
	free(Aij);

}


void thomas(int n, double **a, double *b, double *x) 
{
    /*Thomas Algorithm for tridiagonal systems -- Hoffman 1.8.3*/
    //Input: 
        //n: number of rows in a
        //a: reduced nx3 matrix (changed in-place)
        //b: homogeneity vector from ax = b
    //Output:
        //a: nx3 simplified through Gaussain elim (changed in-place)
        //x: nx1 solution vector to ax = b (changed in-place)
    int i, f = n-1;
    double em; 

    //Forw1rd Elimination
    for (i = 1; i < n; i++)
    {
        em = a[i][0]/a[i-1][1];
        a[i][0] = em;
        a[i][1] = a[i][1] - em*a[i-1][2];
        b[i] = b[i] - a[i][0]*b[i-1];
    }

    //Back Substitution
    x[f] = b[f]/a[f][1];
    for (i = f-1; i >= 0; i--)
    {
        x[i] = (b[i] - a[i][2]*x[i+1])/a[i][1];
    }
    return;
}


double And_mix(double **WA, double *wB){ 
//Anderson Mixing
	//Input globals: PHA, Chi_i, Alpha_i, pot_elec, eta 
	//Output updated WA and wB
	//Updated globals DAs WAs
	double **WAdiff, **WAnew, *wBdiff, *wBnew, **DAnew, *DBnew;
	double **u, **u_temp, *v;
	double psum, lambda;
	double amax, bmax;
	double and_err;

	int i,j,k, n, m;
	int and_Nr, end;
	WAdiff = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAdiff[i] = calloc(Nx*K_i[i], sizeof(double));
	WAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAnew[i] = calloc(Nx*K_i[i], sizeof(double));
	wBdiff=(double *)malloc(sizeof(double)*Nx);
	wBnew=(double *)malloc(sizeof(double)*Nx);
	DAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) DAnew[i] = calloc(Nx*K_i[i], sizeof(double));
	DBnew = (double *)calloc(Nx, sizeof(double)); 

	// Histories
	for (n=1; n<=and_NrMax; n++){ 
		for (k=0; k<NF_N; k++){ 
			for (j=0;j<K_i[k];j++){
				for (i=0;i<Nx;i++){ 
					DAs[n-1][k][j*Nx+i] = DAs[n][k][j*Nx+i];
					DBs[n-1][i] = DBs[n][i];
					WAs[n-1][k][j*Nx+i] = WAs[n][k][j*Nx+i];
					WBs[n-1][i] = WBs[n][i];
				}
			}
		}	
	}

	//Calculate and Update
	end = and_NrMax;

	amax = 0.0; 
	bmax = 0.0;
	for (i=0;i<Nx;i++){
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) {
			WAnew[k][j*Nx+i]    = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i] - eta[i];
			WAdiff[k][j*Nx+i]   = WAnew[k][j*Nx+i] - WA[k][j*Nx+i];
			DAs[end][k][j*Nx+i] = WAdiff[k][j*Nx+i]; //Update last val
			WAs[end][k][j*Nx+i] = WA[k][j*Nx+i]; //Update last val
			amax += WAdiff[k][j*Nx+i]*WAdiff[k][j*Nx+i] / (WA[k][j*Nx+i]*WA[k][j*Nx+i]);
			}
		wBnew[i] = 0;
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wBnew[i] += Chi_i[k][j]*PHA[k][j*Nx + i];
		wBnew[i] -= eta[i];
		wBdiff[i] = wBnew[i]-wB[i];		
		DBs[end][i] = wBdiff[i]; //Update last val 
		WBs[end][i] = wB[i]; //Update last val
		bmax += wBdiff[i]*wBdiff[i] / (wB[i]*wB[i]);
	}	
	and_err = amax+bmax;
	and_err = pow(and_err, 0.50);
	
	if (iter%and_it!=0 || and_NrMax == 0){
		for (i=0; i<Nx; i++){
			psum = 1.0-phA[i]-phB[i];
			for (k=0;k<NF_N;k++) for (j=0; j<K_i[k]; j++) WA[k][j*Nx+i] += wopt*(WAdiff[k][j*Nx+i]- wcmp*psum);
			wB[i] += wopt*(wBdiff[i]-wcmp*psum);
		}
		for (k=0;k<NF_N;k++){free(WAdiff[k]); free(WAnew[k]); free(DAnew[k]);}
		free(WAdiff); free(WAnew); free(DAnew);
		free(wBdiff); free(wBnew); free(DBnew);
		return and_err;
	}	

	and_Nr = fmin(iter-1, and_NrMax); 
	Cs = (double *)calloc(and_Nr, sizeof(double)); 

	//Store histories
	u = INIT_2D(u, and_Nr);
	u_temp = INIT_2D(u_temp, and_Nr);
	v = INIT(and_Nr);
	for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) for (i=0;i<Nx;i++){
		for (m=1;m<=and_Nr;m++) {
			for (n=1;n<=and_Nr;n++) v[m-1] += (DAs[end][k][j*Nx+i] - DAs[end-m][k][j*Nx+i])*DAs[end][k][j*Nx+i];
		}
	}	
	for (i=0;i<Nx;i++) for (m=1;m<=and_Nr;m++){
 		for (n=1;n<=and_Nr;n++) u[m-1][n-1] += (DBs[end][i]-DBs[end-m][i])*(DBs[end][i]-DBs[end-n][i]);
		v[m-1] += (DBs[end][i]-DBs[end-m][i])*DBs[end][i];
	}

	// Compute Coeffs	
	inv(and_Nr, u, u_temp); //u ^-1	
	for (m=0;m<and_Nr;m++) for (n=0;n<and_Nr;n++) u[n][m] = u_temp[m][n]; //uT^-1
	mult(and_Nr, u, v, Cs); //CN = uT^-1 * v
printf("CAS: ");
for (n=0;n<and_Nr;n++) printf("%.6e ", Cs[n]);
printf("\n");

	//Half step
	for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) for (i=0;i<Nx;i++){
		WAnew[k][j*Nx+i] = WAs[end][k][j*Nx+i]; //k+1/2
		DAnew[k][j*Nx+i] = DAs[end][k][j*Nx+i];
		for (n=1;n<=and_Nr;n++) {
			WAnew[k][j*Nx+i] += wand*Cs[end-n]*(WAs[end-n][k][j*Nx+i]-WAs[end][k][j*Nx+i]);
			DAnew[k][j*Nx+i] += wand*Cs[end-n]*(DAs[end-n][k][j*Nx+i]-DAs[end][k][j*Nx+i]);
		}
	}
	for (i=0;i<Nx;i++){
		wBnew[i] = WBs[end][i]; //k+1/2
		DBnew[i] = DBs[end][i];
		for (n=1; n<=and_Nr; n++){
			wBnew[i] += wand*Cs[end-n]*(WBs[end-n][i] - WBs[end][i]);
			DBnew[i] += wand*Cs[end-n]*(DBs[end-n][i] - DBs[end][i]);
		}
	}	

	//Update fields
	lambda =1.0; //1.0-pow(0.9, iter/200);
	for (i=0;i<Nx;i++){
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++){
			WA[k][j*Nx+i] = WAnew[k][j*Nx+i] + lambda * DAnew[k][j*Nx+i];
		}
		wB[i] = wBnew[i] + lambda * DBnew[i];
	}

	//Free
	for (k=0;k<NF_N;k++){free(WAdiff[k]); free(WAnew[k]); free(DAnew[k]);}
	free(WAdiff); free(WAnew); free(DAnew);
	free(wBdiff); free(wBnew); free(DBnew);
	for (n=0;n<and_Nr;n++) {free(u[n]); free(u_temp[n]);}
	free(v);

	return and_err;
}


void inv(int ndim, double **A, double **invA){
//Inverse using Doolittle LU Factorization
//Adapted from Hoffman
//Output inverse of A into invA
	int i,j,k;
	double em;
	double **l, **u;
	double *b, *x, *col; 
	b = INIT(ndim); x = INIT(ndim); col = INIT(ndim);
	u = INIT_2D(u, ndim);
	l = INIT_2D(l, ndim);
	//LU Decomp
	for (k=0;k<ndim-1;k++){
		for (i=k+1;i<ndim;i++){
			em = A[i][k] / A[k][k]; 	
			A[i][k] = em;
			for (j=k+1;j<ndim;j++){
				A[i][j] = A[i][j] - em*A[k][j];
			}
		}	
	}	
	//A-1 = L-1 * U-1
	for (i=0;i<ndim;i++){
		l[i][i] = 1.0;
		for (j=i-1;j>=0;j--) l[i][j] = A[i][j];
		for (j=i;j<ndim;j++) u[i][j] = A[i][j];
	}
	for (i=0;i<ndim;i++){
		for (j=0;j<ndim;j++){ 
			if (j==i) b[j] = 1.0;
			else      b[j] = 0.0;
		}
		solve(ndim, l, b, x);	
		solve(ndim, u, x, col);
		for (j=0;j<ndim;j++) invA[i][j] = col[j];
	}
	return;
}
void solve(int ndim, double **a, double *b, double *x){
	double *bp; bp = INIT(ndim);
	int i,j,k;
	int n = ndim-1;
	bp[0]= b[0];
	for (i=1;i<ndim;i++){
		bp[i] = b[i];
		for (j=0;j<i;j++){
			bp[i] = bp[i] - a[i][j]*bp[j];
		}
	}	
	x[n] = bp[n]/a[n][n];
	for (i=n-1;i>=0;i--){
		x[i] = bp[i];
		for (j=n;j>i;j--){
			x[i] = x[i]-a[i][j]*x[j];
		}
		x[i] = x[i]/a[i][i];
	}
}
void mult(int ndim, double **A, double *B, double *x){
//x = Ab: Square matrix by vector
	int i, j;
	for (i=0;i<ndim;i++){
		x[i] = 0.0;
		for (j=0;j<ndim;j++) x[i] += A[i][j]*B[j];
	}
	return;
}
