/*
to extract data from simulation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>
#include<malloc.h>

#define if_voltage_calcium 0 // 0 for voltage, 1 for calcium

#define oops(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) {oops("error: malloc()"); }

#define NODES  235
#define NODESY 269 // y nodes
#define NODESZ 298 // z nodes.
#define DX 0.33 // mm: was 0.1 mm
#define DY 0.33
#define DZ 0.33

char tissue[NODESZ][NODESY][NODES];
double voltage[NODESZ][NODESY][NODES];

int main(){

FILE *geometry, *prev,*output;
char *str;
int i,h,g; // reserviert fur x,y,z
int filecounter = 200;
int somecounter = 0;

/* read in geometry: if there is no tissue at a location case  */
geometry = fopen("atrium.mat","r");
fread(&tissue,sizeof(unsigned char),NODES*NODESY*NODESZ,geometry);
fclose(geometry);

for(filecounter=0;filecounter<1200;filecounter++){

str = malloc(32*sizeof(char));

if(if_voltage_calcium==0)
sprintf(str,"gunzip court%d.bin.gz",filecounter);
if(if_voltage_calcium==1)
sprintf(str,"gunzip courtccai%d.bin.gz",filecounter);


system(str);
free(str);
str = malloc(32*sizeof(char));

if(if_voltage_calcium==0)
sprintf(str,"court%d.bin",filecounter);
if(if_voltage_calcium==1)
sprintf(str,"courtccai%d.bin",filecounter);

prev = fopen(str,"rb");
free(str);
fread(&voltage,sizeof(double),NODES*NODESY*NODESZ,prev);

/* put your voltage output here */

/*
printf("%d %f %f %f %f %f %f %f %f %f %f %f \n",filecounter,2.505*(filecounter+1),voltage[104][127][23],voltage[134][127][22],voltage[166][124][15]
,voltage[217][124][15],voltage[95][149][28],voltage[101][98][33],voltage[28][95][132],voltage[68][57][19]
,voltage[144][64][182],voltage[190][55][129]);
*/

output = fopen("ap_profiles.dat","a+");
fprintf(output,"%d %f %f %f %f %f %f %f %f %f %f %f \n",filecounter,2.505*(filecounter+1),voltage[104][127][23],voltage[134][127][22],voltage[166][124][15]
,voltage[217][124][15],voltage[95][149][28],voltage[101][98][33],voltage[28][95][132],voltage[68][57][19]
,voltage[144][64][182],voltage[190][55][129]);
fclose(output);

str = malloc(32*sizeof(char));
if(if_voltage_calcium==0)
sprintf(str,"gzip court%d.bin",filecounter);
if(if_voltage_calcium==1)
sprintf(str,"gzip courtccai%d.bin",filecounter);
system(str);
free(str);

}

return 0;
}
