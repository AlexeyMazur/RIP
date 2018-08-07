/*
------------------------------------------------------------------------------

    Processing of experimental Repeat-Induced Point mutation (RIP) data

    Copyright (C) 2018  Alexei Mazour

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

------------------------------------------------------------------------------

    Calculation of Partitioned RIP Propensity (PRP) profiles from multiple
    sequence alignments produced by program Clustal

    Compilation:       gcc -o mk1prp mk1kprp.c
         Usage:       mk1prp XLA2.aln > XLA2.out
   
    Parameters:
                Ref  - reference sequence in the input data.
                Lww1 - window size for data averaging.
                Lww2 - window size for profile smoothing.
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define WR fprintf(stdout,

int Lww1=10;
int Lww2=20;
char *Ref="A_REF";

double *ayoo1(double *indd,int lein,int ww)
{
  static double *rz=NULL;
  double u0;
  int i,j,up,ow,ww2=ww/2,ww2p=ww2+ww%2;

  rz=(double *)calloc(lein,sizeof(double));
  for(i=0;i<lein;i++){
    ow=i-ww2;  if(ow<0) ow=0;
    up=i+ww2p; if(up>lein) up=lein;
    for(u0=0.,j=ow;j<up;j++){ u0+=indd[j]; }
    rz[i]=u0/(up-ow);
  }
  return rz;
}
double *ayoo2(double *indd,int lein,int ww)
{
  static double *rz=NULL;

  int i,j,k,n0,bu1,le1,le2,*lisa=NULL,*lile=NULL,*limw=NULL;
  double u0,*su0=NULL,*v0=NULL,*v1=NULL;
  
  rz=(double *)calloc(4*lein,sizeof(double)); su0=rz+3*lein;
  lisa=(int *)calloc(3*lein,sizeof(int)); lile=lisa+lein; limw=lisa+2*lein;

  if(ww==1){
    memcpy(rz+lein,indd+lein,2*lein*sizeof(double)); v0=rz+lein; v1=v0+lein;
    for(i=0;i<lein;i++){
      for(u0=0.,k=0,j=v1[i];k<v0[i];k++){ u0+=indd[j+k]; } rz[i]=u0;
  } }
  else{
    for(j=i=0;i<lein;i++){ if(fabs(indd[i])>1.e-10){ lisa[j]=i; j++; }} n0=j;
    for(u0=0.,k=j=i=0;i<n0;i++){
      u0+=indd[lisa[i]]; j++;
      if(j==ww){
        su0[k]=u0; lile[k]=lisa[i]-lisa[k]; limw[k]=lisa[i]+lisa[k];
        u0-=indd[lisa[k]]; k++; j--;
    } }
    bu1=k;

    v0=rz+lein; v1=v0+lein;
    for(i=0;i<lisa[0];i++){
      rz[i]=su0[0]; v0[i]=lisa[0]+lile[0]-i; v1[i]=i;
    }
    for(j=0;i<lisa[bu1];i++){
      for(;j<bu1-1;j++){
        le1=abs(2*i-limw[j]); le2=abs(2*i-limw[j+1]);
        if(le2>le1){ break; }
      }
      rz[i]=su0[j]; v0[i]=lile[j]; v1[i]=lisa[j];
    }
    for(   ;i<lein;i++){
      rz[i]=su0[j]; v0[i]=lile[j]; v1[i]=lisa[j];
      if(i-lisa[j]>v0[i]){ v0[i]=i-lisa[j]; }
  } }

  return rz;
}
char *ldfi(FILE *in,int *oule)
{
  static char *wlcch;
  int i,le,fid,utsz=3000000;
  char bf[3000000];

  fid=fileno(in); lseek(fid,0,SEEK_SET);
  for(le=0,i=1;i;le+=i) { i=read(fid,bf,utsz); }
  wlcch=(char *)calloc(le+1,sizeof(char));
  lseek(fid,0,SEEK_SET); read(fid,wlcch,le); le=strlen(wlcch);
  *oule=le; return wlcch;
}
char *nxst(char **bf)
{
  static char rz[1000],*s0=NULL;
  int i=0,j=0;

  for(s0=*bf;isspace(s0[i]);i++);
  for(j=0;!isspace(s0[i])||isblank(s0[i]);i++,j++){ rz[j]=s0[i]; }
  rz[j]='\0'; *bf+=i;

  return rz;
}
char *nxwd(char **bf)
{
  static char rz[80],*s0=NULL;
  int i=0,j=0;

  for(s0=*bf;isspace(s0[i]);i++);
  for(j=0;s0[i]&&!isspace(s0[i]);i++,j++){ rz[j]=s0[i]; }
  rz[j]='\0'; *bf+=i;

  return rz;
}
int main(int argc, char *argv[])
{
  double u0,u1,*msk=NULL,*df=NULL,*rz2=NULL,*v0=NULL,*prp=NULL;
  int i,j,le,ttle,nseq,lseq,lRef,*igs0=NULL,*l0=NULL;
  char *dd0=NULL,*nx0=NULL,*alseq=NULL,*s0=NULL,s1[100],ab[]="CTGA";
  FILE *inp=NULL;

  if(argc!=2){ WR"= ERROR ! Usage: mk1prb clustal_ouput.aln !\n"); exit(1); }
  inp=fopen(argv[1],"rb");
  if(!inp){ WR"= ERROR ! File \"%s\" not found !\n",argv[1]); exit(1); }

  dd0=ldfi(inp,&le); alseq=(char *)calloc(le,sizeof(char));

  nx0=dd0; nxst(&nx0); s0=nx0; lRef=strlen(Ref);
  for(j=1,i=-1;*s0!='*';i++){
    s0=nxst(&nx0); if(!strncmp(s0,Ref,lRef)){ j=0; }
  }
  if(j){ WR"= ERROR ! sequence \"%s\" not found !\n",Ref); exit(1); }
  nseq=i;

  nx0=dd0;
  while((nx0=strstr(nx0,Ref))){ nxwd(&nx0); strcat(alseq,nxwd(&nx0)); }
  lseq=strlen(alseq); ttle=nseq*lseq;

  nx0=dd0; nxst(&nx0);
  for(i=0;i<nseq;i++){
    s0=nxst(&nx0); if(!strncmp(s0,Ref,lRef)) continue;
    strcpy(s1,nxwd(&s0)); strcat(s1," ");
    strcat(alseq,nxwd(&s0)); s0=nx0;
    while((s0=strstr(s0,s1))){
      nxwd(&s0); strcat(alseq,nxwd(&s0));
  } }
  igs0=(int *)calloc(ttle,sizeof(int));
  msk=(double *)calloc(5*lseq,sizeof(double));
    prp=msk+lseq; df=msk+2*lseq;

  for(i=0;i<ttle;i++){ igs0[i]=strchr(ab,alseq[i])-ab; }
  for(i=1;i<nseq;i++){
    l0=igs0+i*lseq;
    for(j=0;j<lseq;j++){ if(igs0[j]!=l0[j]){ df[j]++; } }
  }

  for(i=0;i<lseq;i++){
    if(!strncmp(alseq+i,"CA",2)){ msk[i]=1.; }
    if(!strncmp(alseq+i,"TG",2)){ msk[i+1]=1.; }
  }

  for(u0=1./(nseq-1),i=0;i<lseq;i++){ df[i]*=u0*msk[i]; } 

  rz2=ayoo2(msk,lseq,Lww1); u1=Lww1; u1=(u1-1)/u1;
  memcpy(df+lseq,rz2+lseq,2*lseq*sizeof(double));

  for(v0=rz2+lseq,i=0;i<lseq;i++){ v0[i]=1./(v0[i]+1.e-30); }
  for(i=0;i<lseq;i++){ rz2[i]*=u1*v0[i]; };
  memcpy(prp,rz2,lseq*sizeof(double));

  rz2=ayoo2(df,lseq,1);
  for(v0=rz2+lseq,i=0;i<lseq;i++){ v0[i]=1./(v0[i]+1.e-30); }
  for(i=0;i<lseq;i++){ rz2[i]*=u1*v0[i]; };
  for(i=0;i<lseq;i++){ prp[i]=rz2[i]/prp[i]; };

  if(Lww2>1){ prp=ayoo1(prp,lseq,Lww2); }

  for(i=0;i<lseq;i++){ WR"%7d, %9.2f,\n",(i+1),prp[i]*100); }

  return 0;
}


