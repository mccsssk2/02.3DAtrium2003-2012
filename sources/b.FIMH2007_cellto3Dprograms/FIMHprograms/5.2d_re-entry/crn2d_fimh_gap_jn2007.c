/*

Sanjay Kharche
FIMH 2007 program

Control, Tony AF and Bosch AF

  For uniform Gap Junction remodelled.

*/

// #define AF 3 

#define qten 0 // 0 without q10, 1 with q10.

#define STIM 5*20.0 // is 20.0 for cell model.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include<malloc.h>
 
#define beats 40
#define NODES 375
#define NODESY 375 // y nodes
#define HX 0.1 // mm: was 0.1 mm
#define D 0.3*0.03125 // mm^2 / ms : FIMH value.
// #define D 0.04846 // mm^2 / ms  : value to give CV = 0.33 m/s
// #define D 0.0375 // value that was used originally.
#define length 2 // stimulation length

/* the constants are global */  

double const kqten = 3.0;
  
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

double stimtime;
double resting_at_ap;

double courtemach(double dt, double t, double state[NODES+1][NODESY+1][24], int n, int ny);
double pulse[NODES+1][NODESY+1];  
double new_v[NODES+1][NODESY+1];
double state[NODES+1][NODESY+1][24];
double to_write[NODES+1][NODESY+1];
double old_oi[NODES+1][NODESY+1];
int filecounter;

int AF;

double runtime   = 0.0;

int main()
{ 
    FILE *out,*out1,*gnufile, *timetaken, *tipdata,*position, *diagonal,*ecg;
    
    char *dirname;
 
 double sum = 0.0, distance = 0.0, d2udx2 = 0.0, d2udy2 = 0.0, a =0.0011; // radius of cell mm
 
    double t = 0.0;
    float start, end;
    double BCL = 500.0;
    int n,ny,counter,tipvalue,tip_x,tip_y;
    double writetime,finaltime,diagonaltime=0.0;
    double S1S1;
   
double dt=0.05; /* time step in milliseconds: was 0.005 */  
char *str;
int index;

for(AF=0;AF<3;AF++){

/* start the s1s1 loop here. */

dirname = malloc(64*sizeof(char));
sprintf(dirname,"mkdir AF%d\n",AF);
system(dirname);
free(dirname);

if(AF==0){ // control
S1S1 = 350.0;
runtime = 3000.0;
}

if(AF==1){ // AF Tony
S1S1 = 150.0;
runtime = 3000;
}

if(AF==2){ // AF Bosch
S1S1 = 130;
runtime = 3000;
}


filecounter = 0;

    for(n=0;n<=NODES;n++)
    for(ny=0;ny<=NODESY;ny++){
	state[n][ny][ 0] = -8.118e+01; // U    
	state[n][ny][ 1] = 2.908e-03;  // m  
	state[n][ny][ 2] = 9.649e-01;  // h  
	state[n][ny][ 3] = 9.775e-01;  // j  
	state[n][ny][ 4] = 1.367e-04;  // d  
	state[n][ny][ 5] = 9.996e-01;  // f  
	state[n][ny][ 6] = 3.296e-05;  // xr  
	state[n][ny][ 7] = 1.869e-02;  // xs  
	state[n][ny][ 8] = 1.117e+01;  // nai
	state[n][ny][ 9] = 1.013e-04;  // cai
	state[n][ny][10] = 1.390e+02;  // ki  
	state[n][ny][11] = 1.488e+00;  // caup  
	state[n][ny][12] = 1.488e+00;  // carel  
	state[n][ny][13] = 3.043e-02;  // oa  
	state[n][ny][14] = 9.992e-01;  // oi
	old_oi[n][ny] =  9.992e-01;  // oi  
	state[n][ny][15] = 4.966e-03;  // ua   
	state[n][ny][16] = 9.986e-01;  // ui  
	state[n][ny][17] = 7.755e-01;  // fca  
	state[n][ny][18] = 2.042e-03;  // cmdn  
	state[n][ny][19] = 1.180e-02;  // trpn  
	state[n][ny][20] = 6.504e+00;  // csqn
	state[n][ny][21] = 2.350e-112;  // w  
	state[n][ny][22] = 1.000e+00;  // v  
	state[n][ny][23] = 9.992e-01;  // u  
	new_v[n][ny] = -10.0; // making sure.
    }
	t  = 0.0;  

	writetime = 0;
	diagonaltime = 0.0;
	counter = 1;
	finaltime = 0.0;

	start = omp_get_wtime();
	for (t=0;t<=runtime;)
	{
/*
make 1D simulation first, then copy to all other nodes. note that new_v is NOT initialised now.
*/

if(finaltime<S1S1){
// 1. pulse is going from left to right.
	ny = 0;
	    for(n=0;n<=NODES;n++){
		pulse[n][ny] = 0.0;
		if(t<=2&&n<2) pulse[n][ny] = STIM; // this is what it takes in the single cell model.
	    }

	ny = 0;

 for(n=1;n<NODES;n++)
 {
  new_v[n][ny] = state[n][ny][0]+dt*(D*(state[n+1][ny][0]+state[n-1][ny][0]-2*state[n][ny][0])/(HX*HX)+
      	       - courtemach(dt,t,state,n,ny)+ pulse[n][ny]);
 }

/* boundary conditions */
	    ny = 0;
	    new_v[0][ny] = new_v[1][ny];
	    new_v[NODES][ny] = new_v[NODES-1][ny];

/* update VW counter */

	    ny = 0;
	    if(state[NODES/2][ny][0]<=-60.0&&new_v[NODES/2][ny]>-60.0) counter++;
	    
	    if(counter>1) finaltime = finaltime + dt;

// copy new to old.
	    ny = 0;
	    for(n=0;n<=NODES;n++)
		state[n][ny][0] = new_v[n][ny];

/* copy ALL state to rest of the tissue */

	for(ny=1;ny<=NODESY;ny++) 
		for(n=0;n<=NODES;n++)
			for(index=0;index<=23;index++) 
				state[n][ny][index] = state[n][0][index];

} // end of section 1. this section should run real fast.
else
{
// 2.
    for(ny=0;ny<=NODESY;ny++)
    for(n=0;n<=NODES;n++){
	pulse[n][ny] = 0.0;
	if(AF==0)
	if(finaltime>=S1S1&&finaltime<=(S1S1+dt)&&ny>(NODESY/2+20)) state[n][ny][0] = 100;
	if(AF==1)
	if(finaltime>=S1S1&&finaltime<=(S1S1+dt)&&ny>(NODESY/2+20)) state[n][ny][0] = 100;
	if(AF==2)
	if(finaltime>=S1S1&&finaltime<=(S1S1+dt)&&ny>(NODESY/2+20)) state[n][ny][0] = 100;
}

// start of parallel section.
// the dt is being passed to courtmanche, so it HAS to be shared by all.
#pragma omp parallel shared (state,new_v,pulse,t,dt) private (n,ny) default(none)
{
#pragma omp for schedule(static)
    for(ny=1;ny<NODESY;ny++)
    for(n=1;n<NODES;n++){
	new_v[n][ny] = state[n][ny][0]+dt*(D*(state[n+1][ny][0]+state[n-1][ny][0]-2*state[n][ny][0])/(HX*HX)+
D*(state[n][ny+1][0]+state[n][ny-1][0]-2*state[n][ny][0])/(HX*HX) - courtemach(dt,t,state,n,ny)+ pulse[n][ny]);
}

} // end of parallel.

// boundary conditions. 

for(ny=0;ny<=NODESY;ny++){
    new_v[0][ny] = new_v[1][ny];
    new_v[NODES][ny] = new_v[NODES-1][ny];
}

for(n=0;n<=NODES;n++){
    new_v[n][0] = new_v[n][1];
    new_v[n][NODESY] = new_v[n][NODESY-1];
}


// this has to be after the counter increment.
    for(ny=0;ny<=NODESY;ny++) 
	for(n=0;n<=NODES;n++){
	state[n][ny][0] = new_v[n][ny];
		}

finaltime = finaltime + dt;

} // end of if - else

/*
write diagonal data at a frequency of more than 8 Hz.
*/

if(diagonaltime>=2.5){ // 10 ms

str = malloc(32*sizeof(char));
sprintf(str,"AF%d/diagonal.dat",AF,filecounter);
diagonal = fopen(str,"a+");
free(str);

// diagonal = fopen("diagonal.dat","a+");
fprintf(diagonal,"%f ",t);
for(n=0;n<=NODES;n=n+10) fprintf(diagonal,"%f ",state[n][ny][0]);
fprintf(diagonal,"\n");
fclose(diagonal);
diagonaltime = 0.0;

/*
compute ecg
*/
sum = 0.0;
  for(ny=0;ny<=NODESY;ny++) 
    for(n=0;n<=NODES;n++){
    d2udx2 = (state[n+1][ny][0] + state[n-1][ny][0] - 2.0*state[n][ny][0]);
    d2udy2 = (state[n][ny+1][0] + state[n][ny-1][0] - 2.0*state[n][ny][0]);
    distance = sqrt((n-600)*(n-600)+ny*ny);
    sum -= (d2udx2+d2udy2);
    }
sum = a*a*sum/(64.0*HX*HX*distance);

// ecg = fopen("ecg.dat","a+");
str = malloc(32*sizeof(char));
sprintf(str,"AF%d/ecg.dat",AF,filecounter);
ecg = fopen(str,"a+");
free(str);
fprintf(ecg,"%f %15.10f \n",t,100000*sum);
fclose(ecg);
}
   
if(writetime>=2.5){	
   str = malloc(32*sizeof(char));
   sprintf(str,"AF%d/court2d%d.dat",AF,filecounter);
   out = fopen(str,"w");
   free(str);
   for(ny=0;ny<=NODESY;ny++){
   for(n=0;n<=NODES;n++){
   fprintf(out,"%g ",state[n][ny][0]);
  }
  fprintf(out,"\n");         
}
fclose(out);


   str = malloc(32*sizeof(char));
   sprintf(str,"AF%d/court_cai2d%d.dat",AF,filecounter);
   out = fopen(str,"w");
   free(str);
   for(ny=0;ny<=NODESY;ny++){
   for(n=0;n<=NODES;n++){
   fprintf(out,"%g ",state[n][ny][9]);
  }
  fprintf(out,"\n");         
}
fclose(out);


/*
check if the to_write array is actually populated!!
*/

/*
  for(ny=0;ny<=NODESY;ny++) 
	for(n=0;n<=NODES;n++){
	to_write[n][ny] = state[n][ny][0];
		}
	  str = malloc(32*sizeof(char));
	  sprintf(str,"court2d%d.bin",filecounter);
		out = fopen(str,"wb");
		free(str);
		fwrite(to_write,sizeof(double),(NODES+1)*(NODESY+1),out);
	    fclose(out);
	    */
	writetime = 0.0;
	filecounter++;
} // end of writing.

 t = t+dt;
 writetime = writetime + dt;
 diagonaltime = diagonaltime + dt;
}   // end of time loop.


} // END OF AF LOOPs

end = omp_get_wtime();
timetaken = fopen("time.dat","w");
fprintf(timetaken,"%f %f\n",t,end-start);
return(0);
}
 
double courtemach(double dt, double t, double state[NODES+1][NODESY+1][24], int n, int ny)
{ 
double temporary = 0.0;

	double     V = state[n][ny][0];
	double     m = state[n][ny][1];
	double     h = state[n][ny][2];
	double     j = state[n][ny][3];
	double     d = state[n][ny][4];
	double     f = state[n][ny][5];
	double     xr = state[n][ny][6];
	double     xs = state[n][ny][7];
	double     nai = state[n][ny][8];
	double     cai = state[n][ny][9];
	double     ki = state[n][ny][10];  
	double     caup = state[n][ny][11];
	double     carel = state[n][ny][12];
	double     oa = state[n][ny][13];
	double     oi = state[n][ny][14];
 	double     ua = state[n][ny][15];
 	double     ui = state[n][ny][16];
 	double     fca = state[n][ny][17];
	double     cmdn = state[n][ny][18];
	double     trpn = state[n][ny][19];
	double     csqn = state[n][ny][20];
	double     u = state[n][ny][21];
	double     v = state[n][ny][22];  
  	double     w = state[n][ny][23];  
   
	/* computed quantities */  
   
    double Ek, Ena, Eca,ina,icaL,ik,ibna,ibk,ibca,inak,ik1;  
	double ito,gkur,ikur,ikr,iks;  
	double inaca,naidot,kidot;  
	double caidot,irel;  
	double itr,iup,iupleak;  
	double fnak,caflux,icap,sigma;  
	double itot;
	double a,b,tau,inf,vshift;  
  
  
   
	Ek = 26.71*log(kc/ki);  
	Ena = 26.71*log(nac/nai);  
    Eca = 13.35*log(cac/cai);  
	ina = Cm*gna*m*m*m*h*j*(V-Ena);  
    gkur = 0.005+0.05/(1+exp((V-15)/-13)); 
	ikr = Cm*gkr*xr*(V-Ek)/(1+exp((V+15)/22.4));
	iks = Cm*gks*xs*xs*(V-Ek); 
	ik = ikr + iks;  
	inf = 1/(1+pow(cai/0.00035,1.0));   
		
  	if(AF==0||AF==1) // control or Tony AF
	tau = 2.0;
	if(AF==2) // Bosch AF
	tau = 1.62*2.0;   
 
 	fca = inf + (fca-inf)*exp(-dt/tau);   
	ibna = Cm*gbna*(V-Ena);  
 	ibk = Cm*gbk*(V-Ek);  
  	ibca = Cm*gbca*(V-Eca);  
        sigma = (exp(nac/67.3)-1)/7.0;
 	fnak = 1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0365*sigma*exp(-V*F/(R*T)));  
	icap = Cm*icapbar*cai/(cai+kmcap);   	
 	inaca = Cm*knacalr/(pow(kmnalr,3.0)+pow(nac,3.0))/(kmcalr+cac)/(1+ksatlr*exp((gammalr-1)*V*F/(R*T)))*(nai*nai*nai*cac*exp(V*gammalr*F/(R*T))-nac*nac*nac*cai*exp(V*(gammalr-1)*F/(R*T)));
		
	if(AF==0)
	ito = Cm*gto*oa*oa*oa*oi*(V-Ek);   
	if(AF==1) // tony
	ito = 0.35*Cm*gto*oa*oa*oa*oi*(V-Ek);   
	if(AF==2) // bosch
	ito = 0.16*Cm*gto*oa*oa*oa*oi*(V-Ek);   

   if(AF==0||AF==2) // control or Bosch AF
   ikur = Cm*gkur*ua*ua*ua*ui*(V-Ek);  
   if(AF==1) // Tony AF
   ikur = 1.12*Cm*gkur*ua*ua*ua*ui*(V-Ek); 
   	
        if(AF==0) // control
  	icaL = Cm*gcaL*d*f*fca*(V-ErL);  
        if(AF==1) // Tony AF
  	icaL = 0.36*Cm*gcaL*d*f*fca*(V-ErL);  // see the powerpoint
        if(AF==2) // Bosch AF
  	icaL = 0.36*Cm*gcaL*d*f*fca*(V-ErL);  // see the powerpoint
  	   	
	if(AF==0) // control
 	ik1 = Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));   
        if(AF==1) // Tony AF
 	ik1 = 1.9*Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));
        if(AF==2) // Bosch AF
 	ik1 = 3.35*Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80)));
	
	if(AF==0||AF==2) // control or Bosch AF
   	inak = Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));  
	if(AF==1) // In Tony AF
	inak = 0.88*Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));  

  	 
 
  	
  
 
 	naidot = (-3*inak-3*inaca-ibna-ina)/(F*vi);  	 
 	kidot = (2*inak-ik1-ito-ikur-ikr-iks-ibk)/(F*vi);  
  
 	
   	cmdn = cmdnbar*cai/(cai+kmcmdn);  
   	trpn = trpnbar*cai/(cai+kmtrpn);  
   	csqn = csqnbar*carel/(carel+kmcsqn);  
  
   	
   	irel = grelbar*u*u*v*w*(carel-cai); 
   	iup = iupbar/(1+kmup/cai); 
   	iupleak = kupleak*caup; 
   	itr = (caup-carel)/tautr;  
 
   	
  
   	caidot = ((2*inaca-icap-icaL-ibca)/(2*F*vi)+(iupleak-iup)*vup/vi+irel*vrel/vi)/(1+trpnbar*kmtrpn/(cai*cai+2*cai*kmtrpn+kmtrpn*kmtrpn)+cmdnbar*kmcmdn/(cai*cai+2*cai*kmcmdn+kmcmdn*kmcmdn));   
  
   	
 
   	caup = caup + dt*(iup-itr*vrel/vup-iupleak);  
 
   	 
  
    	carel = carel + dt*((itr-irel)/(1+csqnbar*kmcsqn/(carel*carel+2*carel*kmcsqn+kmcsqn*kmcsqn)));  

   	
  
   	nai = nai + dt*naidot;  
   	ki = ki + dt*kidot;  
   	cai = cai + dt*caidot;  
   
    
  
   	a = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));  
    	if (fabs(V+47.13) < 1e-10) a = 3.2; /* denominator = 0 */ 
     	b = 0.08*exp(-V/11);  
  	tau = 1/(a+b); inf = a*tau;    
    	m = inf + (m-inf)*exp(-dt/tau);  
  
   	 
 
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
	
   	tau = 1/(a+b); inf = a*tau;  
   	j = inf + (j-inf)*exp(-dt/tau);  
  
  
   	vshift = -10;  

   	a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));  
   	b = 0.65/(2.5+exp((V-vshift+72.0)/17.0));   
   	tau = 1/(a+b);
   	inf = 1/(1+exp((V-vshift+10.47)/-17.54)); 
   	oa = inf + (oa-inf)*exp(-Tfac*dt/tau);  
 
 
   	vshift = -10;  
   	a = 1/(18.53+exp((V-vshift+103.7)/10.95));  
   	b = 1/(35.56+exp((V-vshift-8.74)/-7.44));  
   	tau = 1/(a+b);   
   	inf = 1/(1+exp((V-vshift+33.1)/5.3));  
   	oi = inf + (oi-inf)*exp(-Tfac*dt/tau);  
  
   	  
  
   	vshift = -10;
   	a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));
   	b = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
   	tau = 1/(a+b);
   	inf = 1/(1+exp((V-vshift+20.3)/-9.6));
   	ua = inf + (ua-inf)*exp(-Tfac*dt/tau);  
 
   	
   
   	vshift = -10;
   	a = 1/(exp((V-vshift-195)/-28)+21);
   	b = 1/(exp((V-vshift-168)/-16));
   	tau = 1/(a+b);
   	inf = 1/(1+exp((V-vshift-109.45)/27.48));  
   	ui = inf + (ui-inf)*exp(-Tfac*dt/tau);  
    
   	  
  
 	vshift = 0;
   	a = 0.0003*(V-vshift+14.1)/(1-exp((V-vshift+14.1)/-5));
   	b = 0.000073898*(V-vshift-3.3328)/(exp((V-vshift-3.3328)/5.1237)-1);
    	if (fabs(V-vshift+14.1) < 1e-10) a = 0.0015; /* denominator = 0 */
  	if (fabs(V-vshift-3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */  
   	tau = 1/(a+b);
   	inf = 1/(1+exp((V-vshift+14.1)/-6.5));
   	xr = inf + (xr-inf)*exp(-dt/tau);  
    
   	
  
   	vshift = 0;  
   	a = 0.00004*(V-vshift-19.9)/(1-exp((V-vshift-19.9)/-17));
   	b = 0.000035*(V-vshift-19.9)/(exp((V-vshift-19.9)/9)-1);
   	if (fabs(V-vshift-19.9) < 1e-10) /* denominator = 0 */
    {
      a = 0.00068;
      b = 0.000315;  
    }  
 
   
   	tau = 0.5/(a+b);
   	inf = sqrt(1/(1+exp((V-vshift-19.9)/-12.7)));
   	xs = inf + (xs-inf)*exp(-dt/tau);  
    
   	  
  
   	vshift = 0;
 	a = 1/(1+exp((V-vshift+10)/-6.24));
     	tau = a*(1-exp((V-vshift+10)/-6.24))/(0.035*(V-vshift+10));
   	if (fabs(V-vshift+10) < 1e-10) tau = a*4.579; /* denominator = 0 */
   	inf = 1/(1+exp((V-vshift+10)/-8));
   	d = inf + (d-inf)*exp(-dt/tau);
 
   	 
 
   	vshift = 0;
   	inf = exp(-(V-vshift+28)/6.9)/(1+exp(-(V-vshift+28)/6.9));
   	tau = 1.5*2*3/(0.0197*exp(-0.0337*0.0337*(V-vshift+10)*(V-vshift+10))+0.02);  
 
   	f = inf + (f-inf)*exp(-dt/tau);   
   
   	  
 
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

//	V = V - dt*(ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca + pulse(stimtime))/Cm;
	itot = (ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca)/Cm;
//	state[ 0] = V ;
	state[n][ny][ 1] = m ;
	state[n][ny][ 2] = h ;
	state[n][ny][ 3] = j ;
	state[n][ny][ 4] = d ;
	state[n][ny][ 5] = f ;
	state[n][ny][ 6] = xr ;
	state[n][ny][ 7] = xs ;
	state[n][ny][ 8] = nai ;
	state[n][ny][ 9] = cai ;
  	state[n][ny][10] = ki ;
  	state[n][ny][11] = caup ;
  	state[n][ny][12] = carel ;
  	state[n][ny][13] = oa ;
  	state[n][ny][14] = oi ;
	state[n][ny][15] = ua ;
	state[n][ny][16] = ui ;
	state[n][ny][17] = fca ;
	state[n][ny][18] = cmdn ;
  	state[n][ny][19] = trpn ;
  	state[n][ny][20] = csqn ;
	state[n][ny][21] = u ;
	state[n][ny][22] = v ;  
  	state[n][ny][23] = w ;

	return(itot);

}

