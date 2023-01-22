/*
To get filament of the 3D simulation based on RHC algorithm

This program needs:

1. Input files in binary format (big endian)
2. atrium.mat
3. enough memory to store n and n+1 binary data files.

i = x direction
h = y direction
g = z direction

as much as possible, follow the format of the 3D code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>
#include<malloc.h>

#define oops(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) {oops("error: malloc()"); }

#define NODES  235
#define NODESY 269 // y nodes
#define NODESZ 298 // z nodes.
#define DX 0.33 // mm: was 0.1 mm
#define DY 0.33
#define DZ 0.33

#define Viso -20.0 // voltage crossover line

#define tolerance 0.01 // low tolerance for testing, for the proper run, use tolerance = 0.00001 smallest value the determinant can take

char tissue[NODESZ][NODESY][NODES];
double old_voltage[NODESZ][NODESY][NODES];
double new_voltage[NODESZ][NODESY][NODES];

int main(){

FILE *geometry, *prev, *next;
char *str;

int i,h,g; // reserviert fur x,y,z

int success = 0;

double p,q,r,pnew,qnew,rnew,dp,dq,dr; // p is between 0,1 for x, q is between 0 and 1 for y, r is between 0 and 1 for z
double pguess,qguess,rguess; // set to 0.5 always: hopefully the function is simple enough that this works

double a,b,c,d, ainv,binv,cinv,dinv,det; // matrix variables

double x,y,z; // the coordinates if there is a filament location

/* x1 at time n, x2 at time n, x1 at time n1, etc. */
double x1n,x2n,x3n,x4n,x1n1,x2n1,x3n1,x4n1;


int filecounter,iterativesolution,isolcounter,finalfile;

/* initialisations and data read in */

finalfile = 1290; // at rate of 2.5 ms, this is more than 3 seconds.

iterativesolution = 1000000;

/* read in geometry: if there is no tissue at a location case  */
geometry = fopen("atrium.mat","r");
fread(&tissue,sizeof(unsigned char),NODES*NODESY*NODESZ,geometry);
fclose(geometry);



/* start outermost loop */
for(filecounter = 706; filecounter <=706 /* finalfile - 1 */; filecounter++){

// check that this works
str = malloc(32*sizeof(char));
sprintf(str,"court%d.bin",filecounter);
prev = fopen(str,"rb");
free(str);
fread(&old_voltage,sizeof(double),NODES*NODESY*NODESZ,prev);

// check that this works: output a couple of files to screen
str = malloc(32*sizeof(char));
sprintf(str,"court%d.bin",filecounter+1);
prev = fopen(str,"rb");
free(str);
fread(&new_voltage,sizeof(double),NODES*NODESY*NODESZ,prev);

// see if it really works
/*
for(g=1;g<NODESZ-1;g++)
  for(h=1;h<NODESY-1;h++)
  for(i=1;i<NODES-1;i++)
  if(tissue[g][h][i]>0){
printf("%d %f %f \n",tissue[g][h][i],old_voltage[g][h][i],new_voltage[g][h][i]);
}

printf("read in the data files \n");
exit(0);
*/

for(g=1;g<NODESZ-1;g++)
  for(h=1;h<NODESY-1;h++)
  for(i=1;i<NODES-1;i++)
  if(tissue[g][h][i]>0){

// do the z = g*delta plane first

x1n = - old_voltage[g][h][i] + old_voltage[g][h][i+1];
x2n = - old_voltage[g][h][i] + old_voltage[g][h+1][i];
x3n = old_voltage[g][h][i] - old_voltage[g][h][i+1] + old_voltage[g][h+1][i+1] - old_voltage[g][h+1][i];;
x4n = Viso - old_voltage[g][h][i];

x1n1 = - new_voltage[g][h][i] + new_voltage[g][h][i+1];
x2n1 = - new_voltage[g][h][i] + new_voltage[g][h+1][i];
x3n1 = new_voltage[g][h][i] - new_voltage[g][h][i+1] + new_voltage[g][h+1][i+1] - new_voltage[g][h+1][i];;
x4n1 = Viso - new_voltage[g][h][i];

// take a guess
pguess = 0.5;
qguess = 0.5;
rguess = 0.5;

p = -1; // impossible value
q = -1; // impossible value
r = -1; // impossible value

// the iterative loop starts here

for(isolcounter = 0; isolcounter < iterativesolution; isolcounter++){

a = x1n  + qguess*x3n;
b = x2n  + pguess*x3n;
c = x1n1 + qguess*x3n1;
d = x2n1 + pguess*x3n1;

det = a*d - b*c;

if(fabs(det)<=tolerance){
p = -1;
q = -1;
r = -1;
success = 0;
break;
}
else
{
// compute the inverse of the 2x2 matrix
/*
the original matrix is

a  b
c  d

and the inverse is

d/det  -c/det
-b/det  a/det
*/

ainv =  d/det;
binv = -c/det;
cinv = -b/det;
dinv =  a/det;

dp = ainv*pguess + binv*qguess;
dq = cinv*pguess + dinv*qguess;
} // end of fabs tolerance if

if(fabs(dp)<tolerance&&fabs(dq)<tolerance)
{
	success = 1;
	break;
}
else
{
	pguess = pguess + dp;
	qguess = qguess + dq;
}

if(fabs(pguess)>1.0||fabs(qguess)>1.0) break;

printf("iterative loop: %d %f %f %f %f %f %d %d %d \n",isolcounter,dp,dq,pguess,qguess,det,i,h,g);

} // the iterative loops ends here

if(success==1){
if(pguess>=0.0&&pguess<1.0)
p = pguess;
if(qguess>=0.0&&qguess<1.0)
q = qguess;

if(p>=0.0&&p<1.0&&q>=0.0&&q<1.0){
x = i*DX + p;
y = h*DY + q;
z = g*DZ;
printf("%f %f %f\n",x,y,z); // found a filament location!
exit(0);
}
}  // end of if success==1


}  // end of g,h,i loops


} // end of filecounter analysis loop


return 0;
}
