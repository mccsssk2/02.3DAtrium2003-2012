/*
To covert the new binary files in doulble format to
the ones in int format with 0 for nothing, 2 for tissue,
and 4 for propagation.
This code is assuming that we are working with big endian data
i.e. binary files from sun 02 and not horace.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>
#include<malloc.h>

#define NODES  235
#define NODESY 269 // y nodes
#define NODESZ 298 // z nodes.

#define FILE_START 115
#define FILE_END 115

double old_solution[NODESZ][NODESY][NODES];
int solution[NODESZ][NODESY][NODES];
char tissue[NODESZ][NODESY][NODES];

int main(){

int filecounter,i,h,g;

FILE *input, *output,*geometry;

char *str;

geometry = fopen("atrium.mat","r");
fread(&tissue,sizeof(unsigned char),NODES*NODESY*NODESZ,geometry);
fclose(geometry);


for(filecounter=FILE_START;filecounter<=FILE_END;filecounter++){

str = malloc(32*sizeof(char));
sprintf(str,"court%d.bin",filecounter);
input = fopen(str,"rb");
free(str);
fread(&old_solution,sizeof(double),NODES*NODESY*NODESZ,input);
fclose(input);


/*
for(g=0;g<NODESZ;g++)
  for(h=0;h<NODESY;h++)
  for(i=0;i<NODES;i++){
if(old_solution[g][h][i]!=0.0)  printf("%f \n",old_solution[g][h][i]);
}
*/


for(g=0;g<NODESZ;g++)
  for(h=0;h<NODESY;h++)
  for(i=0;i<NODES;i++){
if(tissue[g][h][i]>0){
if(old_solution[g][h][i]<-25)
solution[g][h][i] = 2;
else
solution[g][h][i] = 4;
}
else
solution[g][h][i] = 0;
}

/*
for(g=0;g<NODESZ;g++)
  for(h=0;h<NODESY;h++)
  for(i=0;i<NODES;i++){
if(solution[g][h][i]!=0) printf("%d \n",solution[g][h][i]);
}
*/

str = malloc(32*sizeof(char));
sprintf(str,"NEWcourt%d.bin",filecounter);
output = fopen(str,"wb");
free(str);
fwrite(solution,sizeof(int),NODES*NODESY*NODESZ,output);
fclose(output);
}


return 0;
}
