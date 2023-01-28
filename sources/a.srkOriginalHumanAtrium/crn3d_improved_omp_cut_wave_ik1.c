/*
ssh yddraig@safleprofi.grid.cf.ac.uk

*/

/*
0 = normal
1 = AF by tony:
3 = Beta blocker by Tony
4 = beta blocker AND AF by Tony
(Ip now included).
*/

#define AF 1
#define qten 0 // 0 without q10, 1 with q10.
#define STIM 5*20.0 // is 20.0 for cell model. 3*20 is for control, 5*20 is for hannah AF.
#define PacingFrequency 325.0

#define first_or_continue 0 // 0 for a first run, 1 for continuing runs

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>
#include<malloc.h>

#define oops(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) {oops("error: malloc()"); }

#define beats 40
#define NODES  235
#define NODESY 269 // y nodes
#define NODESZ 298 // z nodes.
#define DX 0.33 // mm: was 0.1 mm
#define DY 0.33
#define DZ 0.33

#define sigma_l 0.03125
#define sigma_t 0.03125

/* the constants are global */  

double const kqten = 3.0;
  
// double const    vcell = 20100.0; /* um3 */   // this variable is never used.
double const    vi =  0.68*20100.0;    
double const    vup = 0.0552*20100.0;  
double const    vrel = 0.0048*20100.0;
double const    T = 310; /* 37 Celcius */  
double const    Tfac = 3.0;   
double const    Csp = 1.0e+6; /* pF/cm2 */  
double const    Cm = 100.0; /* pF */  
double const    F = 96.4867; /* coul/mmol */  
double const    R = 8.3143; /* J K-1 mol-1 */  
double const    kb = 5.4; /* mM */  
double const    nab = 140; /* mM */  
double const    cab = 1.8; /* mM */  
double const    nac = 140.;
double const    cac = 1.8;  
double const    kc = 5.4;    
double const    gna = 7.8; /* nS/pF */   
double const    gto = 0.1652; /* nS/pF */   
double const    gkr = 0.029411765; /* nS/pF */    
double const    gks = 0.12941176; /* nS/pF */  
double const    gcaL = 0.12375; /* nS/pF */   
double const    ErL = 65.0; /* mV */  
double const    gk1 = 0.09; /* nS/pF */   
double const    gbna = 0.0006744375; /* nS/pF */   
double const    gbk = 0.0;  
double const    gbca = 0.001131; /* nS/pF */   
double const    inakbar = 0.59933874; /* pA/pF */    
double const    kmnai = 10.0; /* mM */  
double const    kmko = 1.5; /* mM */
double const    icapbar = 0.275; /* pA/pF */
double const    kmcap = 0.0005; /* mM */
double const    knacalr = 1600.0; /* pA/pF */
double const    kmnalr = 87.5; /* mM */
double const    kmcalr = 1.38; /* mM */
double const    ksatlr = 0.1;
double const    gammalr = 0.35;
double const    trpnbar = 0.070; /* mM */
double const    cmdnbar = 0.050; /* mM */
double const    csqnbar = 10; /* mM */
double const    kmcsqn = 0.8; /* mM */
double const    kmcmdn = 0.00238; /* mM */
double const    kmtrpn = 0.0005; /* mM */
double const    grelbar = 30.0; /* ms-1 */
double const    kmup = 0.00092; /* mM */
double const    iupbar = 0.005; /* mM/ms */
double const    caupmax = 15.0; /* mM */
double const    kupleak = 0.005/15.0; /* ms-1 */  
double const    tautr = 180.0; /* ms */  

double stimtime;
double resting_at_ap;

double courtemach(double dt, double t,int counter1);
void reseted(int counter1);

double *pulse;  
double *new_v;
double **state;
double *sol;

int *san = NULL;
int **nbd = NULL;

double solution[NODESZ][NODESY][NODES];
char tissue[NODESZ][NODESY][NODES];

int filecounter;
int cutted = 0;
double start_t = 0.0;

int main()
{ 
    FILE *out,*geometry,*neighbours,*cont_params,*final_dump,*timedata,*margots_time;
    
    double t = 0.0;
    float start, end;
    double S1S2 = 0.0;

    int i,h,g,counter,counter1,z_max = -100;
    double writetime,margotstime = 0.0;
 
    double d2udx2, d2udy2, d2udz2;
    double Xinv, Yinv, Zinv;
    
//  double dt = 0.005; /* time step in milliseconds */  
    double dt = 0.05;  /* in accordance with what we did for the first round of runs */  
    char *str;
   
    int temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11;
    double temp12;

    start = 0.0; 
    end = 0.0;
    i = 0;
    h = 0;
    g = 0;
  
    writetime = 0;
    counter = 0;
    
    if(AF==0){
    S1S2 = 350.0; // horace
    }
    
    if(AF==1){
    S1S2 = 290.0;
    }

    Zinv = 1.0/(DZ*DZ);
    Yinv = 1.0/(DY*DY);
    Xinv = 1.0/(DX*DX);

geometry = fopen("atrium.mat","r");
fread(&tissue,sizeof(unsigned char),NODES*NODESY*NODESZ,geometry);
fclose(geometry);

counter = 0;

start = omp_get_wtime();

for(g=0;g<NODESZ;g++)
for(h=0;h<NODESY;h++)
    for(i=0;i<NODES;i++)
    if(tissue[g][h][i]>0){
	    if(g>z_max) z_max = g;
     counter++; /* now counter does not change */
     }

MALLOC(san,sizeof(int)*counter);
MALLOC(sol,sizeof(double)*counter); // for the solution

for(i=0;i<counter;i++)
    san[i] = 0;

counter1 = 0;

for(g=0;g<NODESZ;g++)
for(h=0;h<NODESY;h++)
for(i=0;i<NODES;i++){
solution[g][h][i] = 0;
if(tissue[g][h][i]>0){
if(g>=(z_max-10)) // tissue > 79 is SAN. The current thing is max z.
san[counter1] = 1;
counter1++;
}
}

if(counter1!=counter){
	printf("Read in neighbours. Take a look at crn3d_improved_omp_homogenous.c program\n");
	exit(EXIT_FAILURE);
}

counter1 = 0;

/* this point on, dynamics has to follow whats done in neighbourhood.c */

MALLOC(nbd,sizeof(int*)*counter);
for (i=0;i<counter;i++)
{
    MALLOC(nbd[i],sizeof(int)*11);
}
 
for(i=0;i<counter;i++)
    for(h=0;h<11;h++) 
	 nbd[i][h] = -1; /* if there is no neighbour, then the value is -1.*/

neighbours = fopen("neighbours.dat","r");

counter1 = 0;

while(fscanf(neighbours,"%d %d %d %d %d %d %d %d %d %d %d",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6,&temp7,&temp8,&temp9,&temp10,&temp11)!=EOF){
nbd[counter1][0] = temp1; // itself
nbd[counter1][1] = temp2; // x+1
nbd[counter1][2] = temp3; // x-1
nbd[counter1][3] = temp4; // y+1
nbd[counter1][4] = temp5; // y-1
nbd[counter1][5] = temp6; // z+1
nbd[counter1][6] = temp7; // z-1
nbd[counter1][7] = temp8;
nbd[counter1][8] = temp9;
nbd[counter1][9] = temp10;
nbd[counter1][10] = temp11;
counter1++;
}
fclose(neighbours);

if(counter1!=counter){
printf("Read in neighbours. Take a look at crn3d_improved_omp_homogenous.c program\n");
exit(EXIT_FAILURE);
}

counter1 = 0;

MALLOC(pulse,sizeof(double*)*counter);
MALLOC(state,sizeof(double*)*counter);
MALLOC(new_v,sizeof(double*)*counter);

for(i=0;i<counter;i++){
    MALLOC(state[i],sizeof(double)*24);
}

for(counter1=0;counter1<counter;counter1++){
    for(i=0;i<24;i++) 
    	state[counter1][i] = 0.0;     	
        new_v[counter1]    = 0.0;
}

if(first_or_continue==0){ 
counter1=0;
stimtime = 0.0;
  for(g=0;g<NODESZ;g++)
   for(h=0;h<NODESY;h++)
    for(i=0;i<NODES;i++){
    if(tissue[g][h][i]>0)
      {
	state[counter1][ 0] = -8.118e+01; // U    
	state[counter1][ 1] = 2.908e-03;  // m  
	state[counter1][ 2] = 9.649e-01;  // h  
	state[counter1][ 3] = 9.775e-01;  // j  
	state[counter1][ 4] = 1.367e-04;  // d  
	state[counter1][ 5] = 9.996e-01;  // f  
	state[counter1][ 6] = 3.296e-05;  // xr  
	state[counter1][ 7] = 1.869e-02;  // xs  
	state[counter1][ 8] = 1.117e+01;  // nai
	state[counter1][ 9] = 1.013e-04;  // cai
	state[counter1][10] = 1.390e+02;  // ki  
	state[counter1][11] = 1.488e+00;  // caup  
	state[counter1][12] = 1.488e+00;  // carel  
	state[counter1][13] = 3.043e-02;  // oa  
	state[counter1][14] = 9.992e-01;  // oi
	state[counter1][15] = 4.966e-03;  // ua   
	state[counter1][16] = 9.986e-01;  // ui  
	state[counter1][17] = 7.755e-01;  // fca  
	state[counter1][18] = 2.042e-03;  // cmdn  
	state[counter1][19] = 1.180e-02;  // trpn  
	state[counter1][20] = 6.504e+00;  // csqn
	state[counter1][21] = 2.350e-112; // w  
	state[counter1][22] = 1.000e+00;  // v  
	state[counter1][23] = 9.992e-01;  // u  
	new_v[counter1] = 100.0; // making sure that the variable has values.
	pulse[counter1] = 0.0;
	counter1++;
      }
    }
}

  counter1 = 0;
  start = omp_get_wtime();

if(first_or_continue==1){ // read in all the data here

	cont_params = fopen("cont_params.dat","r");
	fscanf(cont_params,"%lf %lf %lf %d %d",&start_t,&stimtime,&writetime,&filecounter,&cutted);
	fclose(cont_params);

	final_dump = fopen("crn_state.dat","r");
	for(counter1=0;counter1<counter;counter1++){
	
	for(i=0;i<24;i++){
	fscanf(final_dump,"%lf ",&temp12);
	state[counter1][i] = temp12;
	}
	
	}
	fclose(final_dump);
counter1 = 0;
}

margotstime = 0.0;

for (t=start_t;t<=start_t+10000.0;t=t+dt)
{
#pragma omp parallel shared(tissue,nbd,state,new_v,pulse,san,solution,counter,Xinv,Yinv,Zinv,t,dt,stimtime) private(counter1,d2udx2,d2udy2,d2udz2) default(none)
{	

#pragma omp for schedule(static)
for(counter1=0;counter1<counter;counter1++){    
   pulse[counter1] = 0.0;
   if(t<=2&&san[counter1]==1) pulse[counter1] = STIM;
}

#pragma omp for schedule(static)
for(counter1=0;counter1<counter;counter1++){
    d2udx2 = 0.0;
    d2udy2 = 0.0;
    d2udz2 = 0.0;

    if(nbd[counter1][1]>=0&&nbd[counter1][2]>=0)	
	d2udx2 = (state[nbd[counter1][1]][0] - 2.0*state[nbd[counter1][0]][0] + state[nbd[counter1][2]][0])*Xinv;
    if(nbd[counter1][3]>=0&&nbd[counter1][4]>=0)
	d2udy2 = (state[nbd[counter1][3]][0] - 2.0*state[nbd[counter1][0]][0] + state[nbd[counter1][4]][0])*Yinv;
    if(nbd[counter1][5]>=0&&nbd[counter1][6]>=0)	
	d2udz2 = (state[nbd[counter1][5]][0] - 2.0*state[nbd[counter1][0]][0] + state[nbd[counter1][6]][0])*Zinv;

/* apply the no flux boundary condition */

    if(nbd[counter1][5]==-1&&nbd[counter1][6]>=0){
	d2udz2 = (state[nbd[counter1][6]][0] + state[nbd[counter1][6]][0] - 2.0*state[nbd[counter1][0]][0])*Zinv;
    }
    if(nbd[counter1][6]==-1&&nbd[counter1][5]>=0){
	d2udz2 = (state[nbd[counter1][5]][0] + state[nbd[counter1][5]][0] - 2.0*state[nbd[counter1][0]][0])*Zinv;
    }
    if((nbd[counter1][5]==-1)&&(nbd[counter1][5]==-1))
    {  
	d2udz2  = 0.0;
    }
    if(nbd[counter1][3]==-1&&nbd[counter1][4]>=0){
	d2udy2 = (state[nbd[counter1][4]][0] + state[nbd[counter1][4]][0] - 2.0*state[nbd[counter1][0]][0])*Yinv;
    }
    if(nbd[counter1][4]==-1&&nbd[counter1][3]>=0){
	d2udy2 = (state[nbd[counter1][3]][0] + state[nbd[counter1][3]][0] - 2.0*state[nbd[counter1][0]][0])*Yinv;
    }
    if((nbd[counter1][3]==-1)&&(nbd[counter1][4]==-1)){
	d2udy2 = 0.0;
    }
    if(nbd[counter1][1]==-1&&nbd[counter1][2]>=0){
	d2udx2 = (state[nbd[counter1][2]][0] + state[nbd[counter1][2]][0] - 2.0*state[nbd[counter1][0]][0])*Xinv;
    }
    if(nbd[counter1][2]==-1&&nbd[counter1][1]>=0){
	d2udx2 = (state[nbd[counter1][1]][0] + state[nbd[counter1][1]][0] - 2.0*state[nbd[counter1][0]][0])*Xinv;
    }
    if((nbd[counter1][1]==-1)&&(nbd[counter1][2]==-1)){
	d2udx2 = 0.0;
    }
     new_v[counter1] = state[counter1][0]+dt*(	
	sigma_l*d2udx2 +
	sigma_l*d2udy2 +
	sigma_l*d2udz2 - courtemach(dt,t,counter1) + pulse[counter1]
	);    
} // end of counter1 diffusion loop

/* update the voltage variable */
#pragma omp for schedule(static)
 for(counter1=0;counter1<counter;counter1++) 
 	state[counter1][0] = new_v[counter1];
} // END OF PARALLEL

// cut wave.
counter1 = 0;

if((t>S1S2)&&cutted==0){
	for(g=0;g<NODESZ;g++)
	   for(h=0;h<NODESY;h++)
	     for(i=0;i<NODES;i++){
	       if(tissue[g][h][i]>0){
	 	if(h>=(NODESY/2+40)){
       	        	state[counter1][0] = 50.0;
		}
	    counter1++;
	   }
	}	
cutted = 1;
}

counter1 = 0;

if(writetime>=2.5){

// write just the state[0] and then reconstruct the data as post-proc. especially useful for Horace, or even sun02??
/*
str = malloc(32*sizeof(char));
sprintf(str,"court%d.dat",filecounter);
out = fopen(str,"w");
free(str);

for(counter1=0;counter1<counter;counter1++)
	fprintf(out,"%10.10f\n",state[counter1][0]);
fclose(out);
*/

/*
for(counter1=0;counter1<counter;counter1++)
sol[counter1] = state[counter1][0];
counter1 = 0;
str = malloc(32*sizeof(char));
sprintf(str,"court%d.bin",filecounter);
out = fopen(str,"wb");
free(str);
fwrite(sol,sizeof(double),counter,out);
fclose(out);
*/

counter1 = 0;
for(g=0;g<NODESZ;g++)
  for(h=0;h<NODESY;h++)
  for(i=0;i<NODES;i++){
  if(tissue[g][h][i]>0){
  	solution[g][h][i] = state[counter1][0];
   counter1++;
}
else
solution[g][h][i] = -150.0;
}

str = malloc(32*sizeof(char));
sprintf(str,"court%d.bin",filecounter);
out = fopen(str,"wb");
free(str);
fwrite(solution,sizeof(double),NODES*NODESY*NODESZ,out);
fclose(out);

/*
counter1 = 0;
for(g=0;g<NODESZ;g++)
  for(h=0;h<NODESY;h++)
  for(i=0;i<NODES;i++){
  if(tissue[g][h][i]>0){
   	solution[g][h][i] = state[counter1][9];
   counter1++;
  }
  else
  solution[g][h][i] = -1.0;
}

str = malloc(32*sizeof(char));
sprintf(str,"courtccai%d.bin",filecounter);
out = fopen(str,"wb");
free(str);
fwrite(solution,sizeof(double),NODES*NODESY*NODESZ,out);
fclose(out);
*/

writetime = 0.0;
filecounter++;

/*
write the state into a data file and then read it in at the next startup.
When you apply S2 at 240 ms, it seems to play up. In any case each
data file takes 5 hours on 8 procs to write. So make a backup of
the state at each file producing time.
*/

/*
ascii data file data for time, counter, 
*/

cont_params = fopen("cont_params.dat","w");
fprintf(cont_params,"%f %f %f %d %d",t,stimtime,writetime,filecounter,cutted);
fclose(cont_params);

/*
the state array: size of which is counter*24
*/

final_dump = fopen("crn_state.dat","w");
for(counter1=0;counter1<counter;counter1++)
{
for(i=0;i<24;i++)
fprintf(final_dump,"%lf\n",state[counter1][i]);
}
fclose(final_dump);

end = omp_get_wtime();

timedata = fopen("timetaken.dat","a+");
fprintf(timedata,"Elapsed time = %f \n",end - start);
fclose(timedata);

} // end of writing.

writetime = writetime + dt;
margotstime = margotstime + dt;
stimtime = stimtime + dt;

if(stimtime>=PacingFrequency+2.0) stimtime = 0.0;

}   // end of time loop.

end = omp_get_wtime();

timedata = fopen("timetaken.dat","a+");
fprintf(timedata,"total: Elapsed time = %f \n",end - start);
fclose(timedata);

return(0);
}
 
double courtemach(double dt, double t,int counter1)  
{ 
	double     V = state[counter1][0];
	double     m = state[counter1][1];
	double     hh = state[counter1][2];
	double     j = state[counter1][3];
	double     d = state[counter1][4];
	double     f = state[counter1][5];
	double     xr = state[counter1][6];
	double     xs = state[counter1][7];
	double     nai = state[counter1][8];
	double     cai = state[counter1][9];
	double     ki = state[counter1][10];  
	double     caup = state[counter1][11];
	double     carel = state[counter1][12];
	double     oa = state[counter1][13];
	double     oi = state[counter1][14];
 	double     ua = state[counter1][15];
 	double     ui = state[counter1][16];
 	double     fca = state[counter1][17];
	double     cmdn = state[counter1][18];
	double     trpn = state[counter1][19];
	double     csqn = state[counter1][20];
	double     u = state[counter1][21];
	double     v = state[counter1][22];  
  	double     w = state[counter1][23];  
   
	/* computed quantities */  
   
	double Ek, Ena, Eca,ina,icaL,ibna,ibk,ibca,inak,ik1;  
	double ito,gkur,ikur,ikr,iks;  
	double inaca,naidot,kidot;  
	double caidot,irel;  
	double itr,iup,iupleak;  
	double fnak,caflux,icap,sigma;  
	double itot;

	/* utility variables */  
  
	double a,b,tau,inf,vshift;  
  
	/* compute Ek */     
	Ek = 26.71*log(kc/ki);  
	/* compute Ena */  
	Ena = 26.71*log(nac/nai);  
	/* compute Eca */
	Eca = 13.35*log(cac/cai);
	/* compute sodium current */
	ina = Cm*gna*m*m*m*hh*j*(V-Ena);  
	/* compute transient outward current */

	ito = Cm*gto*oa*oa*oa*oi*(V-Ek);   

	/* compute ultra-rapid potassium current */  
	gkur = 0.005+0.05/(1+exp((V-15)/-13));  
    
	ikur = Cm*gkur*ua*ua*ua*ui*(V-Ek);
 
	/* compute the rapid delayed outward rectifier K current */ 

      	ikr = Cm*gkr*xr*(V-Ek)/(1+exp((V+15)/22.4));
         	    
	/* compute the slow delayed outward rectifier K current */  
	iks = Cm*gks*xs*xs*(V-Ek); 
//	ik = ikr + iks;  
 
  	/* compute calcium current */  

 	icaL = Cm*gcaL*d*f*fca*(V-ErL);   
  
  	/* update the fca gate immediately */  
  
 	inf = 1/(1+pow(cai/0.00035,1.0));   
 	tau = 2.0;   
 	fca = inf + (fca-inf)*exp(-dt/tau);   

 	/* compute time independent potassium current: Ik1 Kir2.1 gene mutation data */  

	if(AF==0) // control
 	ik1 = Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));   
  if(AF==1) // AF
 	ik1 = 3.5*Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));

 
 	/* compute ibna background current */  
  
 	ibna = Cm*gbna*(V-Ena);  
  
 	/* compute potassium background current */  
  
 	ibk = Cm*gbk*(V-Ek);  
  
 	/* compute ibca background current */  
 
 	ibca = Cm*gbca*(V-Eca);  
    
 	/* compute inak sodium-potassium pump current, LR-style */  
  
 	sigma = (exp(nac/67.3)-1)/7.0;    
 	fnak = 1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0365*sigma*exp(-V*F/(R*T)));  
 	
	inak = Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));  
	  
 	/* compute icap calcium pump current LR-style */  
  
 	icap = Cm*icapbar*cai/(cai+kmcap);   
  
 	/* compute inaca exchanger current LR-style */
  
 	inaca = Cm*knacalr/(pow(kmnalr,3.0)+pow(nac,3.0))/(kmcalr+cac)/(1+ksatlr*exp((gammalr-1)*V*F/(R*T)))*(nai*nai*nai*cac*exp(V*gammalr*F/(R*T))-nac*nac*nac*cai*exp(V*(gammalr-1)*F/(R*T)));

 	/* compute naidot sodium concentration derivative */  

 	naidot = (-3*inak-3*inaca-ibna-ina)/(F*vi);  
    
 	/* compute kidot potassium concentration derivative */  

 	kidot = (2*inak-ik1-ito-ikur-ikr-iks-ibk)/(F*vi);  
  
 	/* calcium buffer dynamics */  
  
   	cmdn = cmdnbar*cai/(cai+kmcmdn);  
   	trpn = trpnbar*cai/(cai+kmtrpn); 
   	csqn = csqnbar*carel/(carel+kmcsqn);  
  
   	/* SR calcium handling */  
  
   	irel = grelbar*u*u*v*w*(carel-cai);
   	iup = iupbar/(1+kmup/cai);
   	iupleak = kupleak*caup;
   	itr = (caup-carel)/tautr;  
 
   	/* compute caidot calcium concentration derivative */  
   	/* using steady-state buffer approximation */  
  
   	caidot = ((2*inaca-icap-icaL-ibca)/(2*F*vi)+(iupleak-iup)*vup/vi+irel*vrel/vi)/(1+trpnbar*kmtrpn/(cai*cai+2*cai*kmtrpn+kmtrpn*kmtrpn)+cmdnbar*kmcmdn/(cai*cai+2*cai*kmcmdn+kmcmdn*kmcmdn));   
  
   	/* update caup calcium in uptake compartment */  
 
   	caup = caup + dt*(iup-itr*vrel/vup-iupleak);  
 
   	/* update carel calcium in release compartment */  
  
    	carel = carel + dt*((itr-irel)/(1+csqnbar*kmcsqn/(carel*carel+2*carel*kmcsqn+kmcsqn*kmcsqn)));  

   	/* update all concentrations */  
  
   	nai = nai + dt*naidot;  
   	ki = ki + dt*kidot;  
   	cai = cai + dt*caidot;  
   
   	/* update ina m gate */  
  
   	a = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));  
    	if (fabs(V+47.13) < 1e-10) a = 3.2; /* denominator = 0 */  
     	b = 0.08*exp(-V/11.0);  	
  	tau = 1/(a+b); 
  	inf = a*tau;
    	m = inf + (m-inf)*exp(-dt/tau);  
  
   	/* update ina h gate */  
 
   	if (V >= -40.0)  
	{  
	    a  = 0.0;
	    b = 1/(0.13*(1+exp((V+10.66)/-11.1)));  
	}  
   	else  
	{  
	    a = 0.135*exp((V+80)/-6.8);  
	    b = 3.56*exp(0.079*V)+3.1e5*exp(0.35*V);
	}  
  
   	tau = 1/(a+b); 
   	inf = a*tau;  
   	hh = inf + (hh-inf)*exp(-dt/tau);  

   	/* update ina j gate */  
   
   	if (V >= -40.0)
	{  
	    a  = 0.0;
	    b = 0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)));  
	}   
   	else
	{
	    a = (-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))*(V+37.78)/(1+exp(0.311*(V+79.23)));
	    b = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));
	}  
	
   	tau = 1/(a+b); 
   	inf = a*tau;
   	j = inf + (j-inf)*exp(-dt/tau);  
  
/*
I have removed the vshift from these variables to make it clearer to 
compare with the Courtmanche paper.
*/

/* some of the parameters are modifed to AF without comments - check with original 
Courtmanche paper. */

   	/* update oa ito gate */   
   	/* corrected for 37 deg */  
   	/* define an voltage shift due to junctional potential  
     	  and effect of Cd++ */  
//   	vshift = -10;

   	a = 0.65/(exp((V+10.0)/-8.5)+exp((V-30.0)/-59.0));  
   	b = 0.65/(2.5+exp((V+82.0)/17.0));   
   	tau = 1/(a+b);
   	if(qten==1) tau = tau/kqten;

   	inf = 1/(1+exp(-(V+20.47)/17.54));  // ito activation
   
   	oa = inf + (oa-inf)*exp(-Tfac*dt/tau);  
   	/* update oi and ois ito gate */
   	/* corrected for 37 deg */    
   	/* define an voltage shift due to junctional potential and effect of Cd++ */  
//   	vshift = -10;  
   	a = 1/(18.53+exp((V+113.7)/10.95));  
   	b = 1/(35.56+exp(-(V+1.26)/7.44));  
   	tau = 1/(a+b); 
   		if(qten==1) tau = tau/kqten;  
   	inf = 1/(1+exp((V+43.1)/5.3));  
   
   	oi = inf + (oi-inf)*exp(-Tfac*dt/tau);  
 
 
   	/* update ua ikur gate */  
   	/* corrected for 37 deg */ 
   	/* define an voltage shift due to junctional potential and effect of Cd++ */  
  
   	vshift = -10;
   	a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp(-(V-vshift-40.0)/59.0));
   	b = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
   	tau = 1/(a+b);
   		if(qten==1) tau = tau/kqten;
   	inf = 1/(1+exp(-(V-vshift+20.3)/9.6));
   	ua = inf + (ua-inf)*exp(-Tfac*dt/tau);  
 
   	/* update ui ikur gate */ 
   	/* corrected for 37 deg */ 
   	/* define an voltage shift due to junctional potential  
     	  and effect of Cd++ */  
   
   	vshift = -10;
   	a = 1.0/(exp(-(V-vshift-195.)/28.)+21.);
   	b = 1.0/(exp(-(V-vshift-168.)/16.));
   	tau = 1/(a+b);
   		if(qten==1) tau = tau/kqten;
   	inf = 1/(1+exp((V-vshift-109.45)/27.48));  
   	ui = inf + (ui-inf)*exp(-Tfac*dt/tau);  
    
   	/* update the xr ikr gate */  
  
 	vshift = 0;
   	a = 0.0003*(V-vshift+14.1)/(1-exp(-(V-vshift+14.1)/5.));
   	b = 0.000073898*(V-vshift-3.3328)/(exp((V-vshift-3.3328)/5.1237)-1.0);
    	if (fabs(V-vshift+14.1) < 1e-10) a = 0.0015; /* denominator = 0 */
  	if (fabs(V-vshift-3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */  
   	tau = 1/(a+b);
   	inf = 1/(1+exp(-(V-vshift+14.1)/6.5));
   	xr = inf + (xr-inf)*exp(-dt/tau);  
    
   	/* update the xs  iks gate */  
  
   	vshift = 0;  
   	a = 0.00004*(V-vshift-19.9)/(1-exp(-(V-vshift-19.9)/17.));
   	b = 0.000035*(V-vshift-19.9)/(exp((V-vshift-19.9)/9)-1.);
   	if (fabs(V-vshift-19.9) < 1e-10) /* denominator = 0 */
	    {
	      a = 0.00068;
	      b = 0.000315;  
	    }  
   	/* tau reduced by 50% as described in manuscript of courtmanche. */  
   	tau = 0.5/(a+b);
   	inf = sqrt(1/(1+exp((V-vshift-19.9)/-12.7)));
   	xs = inf + (xs-inf)*exp(-dt/tau);  

   	/* update icaL d gate */  
  
   	vshift = 0;
 	a = 1/(1+exp((V-vshift+10)/-6.24));
     	tau = a*(1-exp((V-vshift+10)/-6.24))/(0.035*(V-vshift+10));
   	if (fabs(V-vshift+10) < 1e-10) tau = a*4.579; /* denominator = 0 */
   	inf = 1/(1+exp((V-vshift+10)/-8));
   	d = inf + (d-inf)*exp(-dt/tau);
 
   	/* update icaL f gate */  
 
//   	vshift = 0;
   	inf = exp(-(V+28)/6.9)/(1+exp(-(V+28)/6.9));   	
   	tau = 9.0/(0.0197*exp(-0.0337*0.0337*(V+10.0)*(V+10.0))+0.02);  

   	f = inf + (f-inf)*exp(-dt/tau);   
   
   	/* update the SR gating variables */   
   	/* caflux is expected in umoles/ms, hence the factor of 1000 */
   	/* 1e-15 is used to scale the volumes! */   
 
   	vshift=0;
   	caflux = 1e3*(1e-15*vrel*irel-1e-15*(0.5*icaL-0.1*2*inaca)/(2*F));
   	inf = 1/(1+exp(-(caflux-3.4175e-13-vshift)/13.67e-16));
   	tau = 8.0;
   	u = inf + (u-inf)*exp(-dt/tau);
   	inf = 1-1/(1+exp(-(caflux-6.835e-14-vshift)/13.67e-16));
   	tau = 1.91+2.09/(1+exp(-(caflux-3.4175e-13-vshift)/13.67e-16));  
   	v = inf + (v-inf)*exp(-dt/tau);
   	inf = 1-1/(1+exp(-(V-40)/17.0));
   	tau = 6.0*(1-exp(-(V-7.9)/5.0))/(1+0.3*exp(-(V-7.9)/5.0))/(V-7.9); 
   	if (fabs(V-7.9) < 1e-10) tau = 6*0.2/1.3;  
   
   	w = inf + (w-inf)*exp(-dt/tau);  

   	/* update membrane voltage */  

//	V = V - dt*(ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca + pulse(t))/Cm;  
	itot = (ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca)/Cm; 
//	state[ 0] = V ;
	state[counter1][ 1] = m ;
	state[counter1][ 2] = hh ;
	state[counter1][ 3] = j ;
	state[counter1][ 4] = d ; 
	state[counter1][ 5] = f ;
	state[counter1][ 6] = xr ;  
	state[counter1][ 7] = xs ;  
	state[counter1][ 8] = nai ;  
	state[counter1][ 9] = cai ;  
  	state[counter1][10] = ki ;  
  	state[counter1][11] = caup ;  
  	state[counter1][12] = carel ;  
  	state[counter1][13] = oa ;
  	state[counter1][14] = oi ;
	state[counter1][15] = ua ;
	state[counter1][16] = ui ;
	state[counter1][17] = fca ;  
	state[counter1][18] = cmdn ;  
  	state[counter1][19] = trpn ;  
  	state[counter1][20] = csqn ;  
	state[counter1][21] = u ;
	state[counter1][22] = v ;  
  	state[counter1][23] = w ;  

	return(itot);

}




void reseted(int counter1){

	state[counter1][ 0] = -8.118e+01; // U    
	state[counter1][ 1] = 2.908e-03;  // m  
	state[counter1][ 2] = 9.649e-01;  // h  
	state[counter1][ 3] = 9.775e-01;  // j  
	state[counter1][ 4] = 1.367e-04;  // d  
	state[counter1][ 5] = 9.996e-01;  // f  
	state[counter1][ 6] = 3.296e-05;  // xr  
	state[counter1][ 7] = 1.869e-02;  // xs  
	state[counter1][ 8] = 1.117e+01;  // nai
	state[counter1][ 9] = 1.013e-04;  // cai
	state[counter1][10] = 1.390e+02;  // ki  
	state[counter1][11] = 1.488e+00;  // caup  
	state[counter1][12] = 1.488e+00;  // carel  
	state[counter1][13] = 3.043e-02;  // oa  
	state[counter1][14] = 9.992e-01;  // oi
	state[counter1][15] = 4.966e-03;  // ua   
	state[counter1][16] = 9.986e-01;  // ui  
	state[counter1][17] = 7.755e-01;  // fca  
	state[counter1][18] = 2.042e-03;  // cmdn  
	state[counter1][19] = 1.180e-02;  // trpn  
	state[counter1][20] = 6.504e+00;  // csqn
	state[counter1][21] = 2.350e-112;  // w  
	state[counter1][22] = 1.000e+00;  // v  
	state[counter1][23] = 9.992e-01;  // u  
	new_v[counter1] = -10.0; // making sure.

}
