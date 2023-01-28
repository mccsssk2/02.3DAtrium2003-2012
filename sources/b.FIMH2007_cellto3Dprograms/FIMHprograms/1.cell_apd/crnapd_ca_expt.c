/*
Sanjay Kharche

FIMH 2007 program.

1. Make a table with the following data:
Tony AF: i(K1) up 90%; i(CaL) down 64%; i(to) down 65%; i(K,sus) up 12, i(p) down 12%.
Bosch AF:
  
*/

// #define AF 0

/*
AF = 0 is for control
AF = 1 is for AF Tony
AF = 2 is for AF Bosch
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>

#define beats 10

#define delta 200.0

/* the constants are global */  
  
double const    vcell = 20100.0; /* um3 */  
double const    vi =  0.68*20100.0;    
double const    vup = 0.0552*20100.0;  
double const    vrel = 0.0048*20100.0;
double const    T = 310; /* 37 Celcius */  
double const    Tfac = 3;   
double const    Csp = 1e+6; /* pF/cm2 */  
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
  
  
double apd[100],dis[100],told[100],tnew[100];  
int apdcounter,imax,someothercounter,numS1,keych;
double eold,stimtime;
double resting_at_ap;

double S2times[] = {1000,-1,-1};

double Ek, Ena, Eca,ina,icaL,ik,ibna,ibk,ibca,inak,ik1;  
double ito,gkur,ikur,ikr,iks;  
double inaca,naidot,kidot; 

double dvdtmax = -100.0;

void courtemach(double dt, double t, double state[24], double dvdt);  
double pulse( double x);  

int AF;

int main()  
{  
  double t;
  double dt= 0.005; /* time step in milliseconds */  
  double BCL = 1000.0;
  double state[24];  
  int outputcounter = 0;

  FILE *out;  
  FILE *out_two;
  char *str;
  
  double dvdt,vmin, vmax;

  for(AF=0;AF<3;AF++){

  someothercounter = 0;
  outputcounter = 0;
do{  

for(apdcounter=0;apdcounter<101;apdcounter++){
apd[100]=0;
dis[100]=0;
told[100]=0;
tnew[100]=0;
}

  apdcounter = 0;
  stimtime = 0.0;
  
  dvdtmax = -1000.0;
  vmin = 1000.0;
  vmax = -1000.0;
  
	/* steady-state initial conditions as used in the manuscript */  
  
    
	state[ 0] = -8.118e+01; // U    
	state[ 1] = 2.908e-03;  // m  
	state[ 2] = 9.649e-01;  // h  
	state[ 3] = 9.775e-01;  // j  
  
	state[ 4] = 1.367e-04;  // d  
  
	state[ 5] = 9.996e-01;  // f  
  
	state[ 6] = 3.296e-05;  // xr  
  
	state[ 7] = 1.869e-02;  // xs  
  
	state[ 8] = 1.117e+01;  // nai  
  
	state[ 9] = 1.013e-04;  // cai  
  
	state[10] = 1.390e+02;  // ki
	state[11] = 1.488e+00;  // caup
	state[12] = 1.488e+00;  // carel
	state[13] = 3.043e-02;  // oa
	state[14] = 9.992e-01;  // oi 
	state[15] = 4.966e-03;  // ua 
	state[16] = 9.986e-01;  // ui  
	state[17] = 7.755e-01;  // fca  
	state[18] = 2.042e-03;  // cmdn  
  
	state[19] = 1.180e-02;  // trpn  
  
	state[20] = 6.504e+00;  // csqn  
  
	state[21] = 2.350e-112;  // w  
  
	state[22] = 1.000e+00;  // v  
  
	state[23] = 9.992e-01;  // u  
  
	t  = 0.0;  
  

numS1 = 0;	

str = malloc(32*sizeof(char));
sprintf(str,"ap_profile%d.dat",AF);
out_two = fopen(str,"w");
free(str);

	for (t=0;t<=beats*BCL+S2times[someothercounter]+1000;)  
	{  
	eold = state[0];
	    courtemach(dt, t, state,dvdt);
		t = t + dt;

	if(stimtime>=0&&stimtime<=dt){ 
   	resting_at_ap = state[0];
        told[apdcounter] = t;
      
      }

		stimtime = stimtime + dt;
  
  
  if(stimtime>=BCL&&numS1<beats-1){ stimtime = 0.0; numS1++; /*printf("%d \n",numS1);*/}
  
  if((numS1==(beats-1))&&(stimtime>S2times[someothercounter])){stimtime = 0.0; numS1++;}

  outputcounter++;
	   /* measure apd di */
if ((t>=beats*BCL-200)&&(outputcounter%200==0))
	fprintf(out_two,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t %10.10f\t%10.10f\t%10.10f\t %10.10f\t %10.10f\t %10.10f\t %10.10f\t%d\n",t,state[0],state[9],icaL,ik1,ito,ikur,inaca,inak,iks,ikr,numS1);
   
if(eold>resting_at_ap*0.9&&state[0]<=resting_at_ap*0.9){ tnew[apdcounter] = t; apdcounter++;}

if(numS1>9&&numS1<=10){
 if(state[0]<vmin) vmin = state[0];
        if(state[0]>vmax) vmax = state[0];
//      if(dvdtmax<(state[0] - eold)/dt) dvdtmax = (- state[0] + eold)/dt;
//        printf("%f %f\n",dvdt,dvdtmax);
// 	printf("%f \n",(state[0] - eold)/dt);
}
  
	}   // end of time loop.

  /*END OF MAIN LOOP*/
	imax = apdcounter;


	str = malloc(32*sizeof(char));
	sprintf(str,"apd90%d.dat",AF);
	out = fopen(str,"w");  
	free(str);

	if(out==NULL)
	{  
		printf("  can't open %s \n" , "test.dat");  
		return 0;  
  	}

	for(apdcounter=1;apdcounter<imax;apdcounter++){
	  apd[apdcounter] = tnew[apdcounter] - told[apdcounter];
	  dis[apdcounter] = told[apdcounter] - tnew[apdcounter-1];
	  if(apdcounter>=imax-1){
	    fprintf(out,"%f\t%f\t%d\n",apd[apdcounter],dis[apdcounter],apdcounter);
	    fprintf(out,"%f\t%f\t%d\n",vmin,vmax,dvdtmax);
	    printf("%f \n",dvdtmax);
	    }
	}

    fclose(out_two);
	fclose(out);  
  someothercounter++;
  
/* end of restitution loop. */  
  }while(S2times[someothercounter]>0);

  } // end of AF loop
  
return 0;  
}  

  
void courtemach(double dt, double t, double state[24], double dvdt)  
{  
 
  
	double     V = state[0];  
 	double     m = state[1];  
	double     h = state[2];
	double     j = state[3];  
	double     d = state[4];  
	double     f = state[5];  
	double     xr = state[6];  
	double     xs = state[7];  
	double     nai = state[8];  
	double     cai = state[9];  
	double     ki = state[10];  
	double     caup = state[11];  
	double     carel = state[12];  
	double     oa = state[13];  
	double     oi = state[14];  
 	double     ua = state[15];  
 	double     ui = state[16];  
 	double     fca = state[17];  
	double     cmdn = state[18];    
	double     trpn = state[19];  
	double     csqn = state[20];  
	double     u = state[21];  
	double     v = state[22];  
	double     w = state[23];  

	/* computed quantities */  
  /*
	double Ek, Ena, Eca,ina,icaL,ik,ibna,ibk,ibca,inak,ik1;  
  
  
  
	double ito,gkur,ikur,ikr,iks;  
  
  
  
	double inaca,naidot,kidot;  
  
  */
  
	double caidot,irel;  
  
  
  
	double itr,iup,iupleak;  
  
  
  
	double fnak,caflux,icap,sigma;  
  
  
  
	/* utility variables */  
  
  
  
	double a,b,tau,inf,vshift;  
  
  
  
  
  
 
	/* compute Ek */  
  
  
  
	Ek = 26.71*log(kc/ki);  
  
	/* compute Ena */  
  
  
  
	Ena = 26.71*log(nac/nai);  
  
  
  
    
  
  
  
	/* compute Eca */  
  
  
  
	Eca = 13.35*log(cac/cai);  
  
	/* compute sodium current */  
  
	ina = Cm*gna*m*m*m*h*j*(V-Ena);  
 
	/* compute transient outward current */  
  
 /* corrections added: all from the powerpoint data v*/
	if(AF==0)
	ito = Cm*gto*oa*oa*oa*oi*(V-Ek);   
	if(AF==1) // tony
	ito = 0.35*Cm*gto*oa*oa*oa*oi*(V-Ek);   
	if(AF==2) // bosch
	ito = 0.16*Cm*gto*oa*oa*oa*oi*(V-Ek);   


	/* compute ultra-rapid potassium current */  

  
	gkur = 0.005+0.05/(1+exp((V-15)/-13));  
    

	/* sustained outward current */

	if(AF==0||AF==2) // control or Bosch AF
	ikur = Cm*gkur*ua*ua*ua*ui*(V-Ek);  
   if(AF==1) // Tony AF
	ikur = 1.12*Cm*gkur*ua*ua*ua*ui*(V-Ek); 

	ikr = Cm*gkr*xr*(V-Ek)/(1+exp((V+15)/22.4));   
  
	/* compute the slow delayed outward rectifier K current */  
  
	iks = Cm*gks*xs*xs*(V-Ek);  
 
  
	ik = ikr + iks;  
 
  
  	/* compute calcium current */  
  
  
  if(AF==0) // control
  	icaL = Cm*gcaL*d*f*fca*(V-ErL);  
  if(AF==1) // Tony AF
  	icaL = 0.36*Cm*gcaL*d*f*fca*(V-ErL);  // see the powerpoint
 
  if(AF==2) // Bosch AF
  	icaL = 0.36*Cm*gcaL*d*f*fca*(V-ErL);  // see the powerpoint
 

  	/* update the fca gate immediately */  
  
 	inf = 1/(1+pow(cai/0.00035,1.0));   
  	
	if(AF==0||AF==1) // control or Tony AF
	tau = 2.0;     	
	if(AF==2) // Bosch AF
	tau = 1.62*2.0;   
  	
	fca = inf + (fca-inf)*exp(-dt/tau);   
   
 	/* compute time independent potassium current */  
 
	if(AF==0) // control
 	ik1 = Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));   
  if(AF==1) // Tony AF
 	ik1 = 1.9*Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));
  if(AF==2) // Bosch AF
 	ik1 = 3.35*Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));
 
 
  /* compute ibna background current */  
  	ibna = Cm*gbna*(V-Ena);  
  
 	/* compute potassium background current */  
  	ibk = Cm*gbk*(V-Ek);  
   
 	/* compute ibca background current */  
  	ibca = Cm*gbca*(V-Eca);  
   
 	/* compute inak sodium-potassium pump current, LR-style */  
  	sigma = (exp(nac/67.3)-1)/7.0;  
   
 	fnak = 1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0365*sigma*exp(-V*F/(R*T)));  

	/* THIS IS YOUR IP */
	if(AF==0||AF==2) // control or Bosch AF
   	inak = Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));  
	if(AF==1) // In AF
		inak = 0.88*Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));  

 	/* compute icap calcium pump current LR-style */  
  
 	icap = Cm*icapbar*cai/(cai+kmcap);   
    
 	/* compute inaca exchanger current LR-style */  
    
 	inaca = Cm*knacalr/(pow(kmnalr,3.0)+pow(nac,3.0))/(kmcalr+cac)/  
   	  (1+ksatlr*exp((gammalr-1)*V*F/(R*T)))*  
   	   (nai*nai*nai*cac*exp(V*gammalr*F/(R*T))-nac*nac*nac*cai*  
   	    exp(V*(gammalr-1)*F/(R*T)));    
  
   	/* compute naidot sodium concentration derivative */  
    
 	naidot = (-3*inak-3*inaca-ibna-ina)/(F*vi);  
    
 	/* compute kidot potassium concentration derivative */  
  
 	kidot = (2*inak-ik1-ito-ikur-ikr-iks-ibk)/(F*vi);  
  
  
 	/* calcium buffer dynamics */  
  
   	cmdn = cmdnbar*cai/(cai+kmcmdn);  
  
   	trpn = trpnbar*cai/(cai+kmtrpn);  
  
   	csqn = csqnbar*carel/(carel+kmcsqn);  
  
   	/* SR calcium handling */  

if(t<10000+delta)  
   	irel = grelbar*u*u*v*w*(carel-cai);
   else
   	irel = grelbar*(carel - cai)*exp(-(t-10000-delta)/4.0)*(1.0 - exp(-(t-10000-delta)/4.0));
   	
   	iup = iupbar/(1+kmup/cai);  
  
   	iupleak = kupleak*caup;  
  
   	itr = (caup-carel)/tautr;

   	/* compute caidot calcium concentration derivative */  
  	/* using steady-state buffer approximation */  
  
   	caidot = ((2*inaca-icap-icaL-ibca)/(2*F*vi)+  
 	    (iupleak-iup)*vup/vi+irel*vrel/vi)/  
      (1+trpnbar*kmtrpn/(cai*cai+2*cai*kmtrpn+kmtrpn*kmtrpn)+  
	       cmdnbar*kmcmdn/(cai*cai+2*cai*kmcmdn+kmcmdn*kmcmdn));   
  
  	/* update caup calcium in uptake compartment */  
  
   	caup = caup + dt*(iup-itr*vrel/vup-iupleak);  
  
   	/* update carel calcium in release compartment */  
  
   	carel = carel + dt*((itr-irel)/(1+csqnbar*kmcsqn/  
				  (carel*carel+2*carel*kmcsqn+kmcsqn*kmcsqn)));  
  
   	/* update all concentrations */  
   	nai = nai + dt*naidot;  
   	ki = ki + dt*kidot;  
   	cai = cai + dt*caidot;  

   	/* update ina m gate */  
   	a = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));  
  
  
  
   	if (fabs(V+47.13) < 1e-10) a = 3.2; /* denominator = 0 */  
  
  
  
   	b = 0.08*exp(-V/11);  
  
  
  
   	tau = 1/(a+b); inf = a*tau;  
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
  
  
  
   	tau = 1/(a+b); inf = a*tau;  
  
  
  
   	h = inf + (h-inf)*exp(-dt/tau);  
  
  
  
    
  
  
  
   	/* update ina j gate */  
  
  
  
   	if (V >= -40.0)  
  
  
  
    {  
  
  
  
      a  = 0.0;  
  
  
  
      b = 0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)));  
  
  
  
    }  
  
  
  
   	else  
  
  
  
    {  
  
  
  
      a = (-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))*(V+37.78)/  
  
  
  
	(1+exp(0.311*(V+79.23)));  
  
  
  
      b = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));  
  
  
  
    }  
  
  
  
   	tau = 1/(a+b); inf = a*tau;  
  
  
  
   	j = inf + (j-inf)*exp(-dt/tau);  
  
  
  
    
  
  
  
   	/* update oa ito gate */  
  
  
  
   	/* corrected for 37 deg */  
  
  
  
   	/* define an voltage shift due to junctional potential  
  
  
  
    	  and effect of Cd++ */  
  
  
  
   	vshift = -10;  
  
  
  
   	a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));  
  
  
  
   	b = 0.65/(2.5+exp((V-vshift+72.0)/17.0));   
  
  
  
   	tau = 1/(a+b);   
  
  
  
   	inf = 1/(1+exp((V-vshift+10.47)/-17.54));  
  
  
  
   	oa = inf + (oa-inf)*exp(-Tfac*dt/tau);  
  
  
  
  
  
  
  
   	/* update oi and ois ito gate */  
  
  
  
   	/* corrected for 37 deg */  
  
  
  
   	/* define an voltage shift due to junctional potential  
  
  
  
     and effect of Cd++ */  
  
  
  
   	vshift = -10;  
  
  
  
   	a = 1/(18.53+exp((V-vshift+103.7)/10.95));  
  
  
  
   	b = 1/(35.56+exp((V-vshift-8.74)/-7.44));  
  
  
  
   	tau = 1/(a+b);   
  
  
  
   	inf = 1/(1+exp((V-vshift+33.1)/5.3));  
  
  
  
   	oi = inf + (oi-inf)*exp(-Tfac*dt/tau);  
  
  
  
  
  
  
  
   	/* update ua ikur gate */  
  
  
  
   	/* corrected for 37 deg */  
  
  
  
   	/* define an voltage shift due to junctional potential  
  
  
  
    	  and effect of Cd++ */  
  
  
  
   	vshift = -10;  
  
  
  
   	a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));  
  
  
  
   	b = 0.65/(2.5+exp((V-vshift+72.0)/17.0));   
  
  
  
   	tau = 1/(a+b);   
  
  
  
   	inf = 1/(1+exp((V-vshift+20.3)/-9.6));  
  
  
  
   	ua = inf + (ua-inf)*exp(-Tfac*dt/tau);  
  
  
  
  
  
  
  
   	/* update ui ikur gate */  
  
  
  
   	/* corrected for 37 deg */  
  
  
  
   	/* define an voltage shift due to junctional potential  
  
  
  
    	  and effect of Cd++ */  
  
  
  
   	vshift = -10;  
  
  
  
   	a = 1/(exp((V-vshift-195)/-28)+21);  
  
  
  
   	b = 1/(exp((V-vshift-168)/-16));  
  
  
  
   	tau = 1/(a+b);   
  
  
  
   	inf = 1/(1+exp((V-vshift-109.45)/27.48));  
  
  
  
   	ui = inf + (ui-inf)*exp(-Tfac*dt/tau);  
  
  
  
  
  
  
  
   	/* update the xr ikr gate */  
  
  
  
   	vshift = 0;  
  
  
  
   	a = 0.0003*(V-vshift+14.1)/(1-exp((V-vshift+14.1)/-5));  
  
  
  
   	b = 0.000073898*(V-vshift-3.3328)/(exp((V-vshift-3.3328)/5.1237)-1);  
  
  
  
   	if (fabs(V-vshift+14.1) < 1e-10) a = 0.0015; /* denominator = 0 */   
  
  
  
   	if (fabs(V-vshift-3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */   
  
  
  
   	tau = 1/(a+b);   
  
  
  
   	inf = 1/(1+exp((V-vshift+14.1)/-6.5));  
  
  
  
   	xr = inf + (xr-inf)*exp(-dt/tau);  
  
  
  
    
  
  
  
   	/* update the xs ikr gate */  
  
  
  
   	vshift = 0;  
  
  
  
   	a = 0.00004*(V-vshift-19.9)/(1-exp((V-vshift-19.9)/-17));  
  
  
  
   	b = 0.000035*(V-vshift-19.9)/(exp((V-vshift-19.9)/9)-1);  
  
  
  
   	if (fabs(V-vshift-19.9) < 1e-10) /* denominator = 0 */  
  
  
  
    {  
  
  
  
      a = 0.00068;   
  
  
  
      b = 0.000315;   
  
  
  
    }  
  
  
  
   	/* tau reduced by 50% as described in manuscript */  
  
  
  
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
  
  
  
   	vshift = 0;  
  
  
  
   	inf = exp(-(V-vshift+28)/6.9)/(1+exp(-(V-vshift+28)/6.9));  
  
  
  
   	tau = 1.5*2*3/(0.0197*exp(-0.0337*0.0337*(V-vshift+10)*  
  
  
  
			    (V-vshift+10))+0.02);  
  
  
  
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
  
  	dvdt = -(ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca + pulse(stimtime))/Cm;

  if(numS1>9&&numS1<=10)
      if(dvdtmax<(dvdt)) dvdtmax = dvdt;
  
   	V = V + dt*dvdt;  
	state[ 0] = V ;  
	state[ 1] = m ;  
	state[ 2] = h ;  
	state[ 3] = j ;  
	state[ 4] = d ;  
	state[ 5] = f ;  
	state[ 6] = xr ;  
	state[ 7] = xs ;  
	state[ 8] = nai ;  
	state[ 9] = cai ;  
  	state[10] = ki ;  
  	state[11] = caup ;  
  	state[12] = carel ;  
  	state[13] = oa ;  
  	state[14] = oi ;    
	state[15] = ua ;  
	state[16] = ui ;  
	state[17] = fca ;  
	state[18] = cmdn ;  
  	state[19] = trpn ;  
  	state[20] = csqn ;  
	state[21] = u ;  

	state[22] = v ;  
  	state[23] = w ;  

	return ;  
}  

double pulse( double x)
{
	double y;
	if(x>=0&&x<=2) y = -2e3;
		else 
	y = 0.0;

	return y ;
}  
