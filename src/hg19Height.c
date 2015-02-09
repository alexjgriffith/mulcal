/* This file is part of peakAnalysis,
   http://github.com/alexjgriffith/alpha-score/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include "hg19Height.h"


int getChromosome(char value[10])
{
  int i;
  for(i=0;i<25;i++)
    if(!strcmp(value,hg19Chrom[i].value))
      return i;
  return -1;
}

int getChromosomeShort(char value[2])
{
  int i;
  for(i=0;i<25;i++)
    if(value[0]==hg19Chrom[i].shortValue[0] && value[1]==hg19Chrom[i].shortValue[1])
      return i;
  return -1;
}

void getChromosomeValue(int order,char * value)
{
  int i;
  for(i=0;i<25;i++)
    if(order==hg19Chrom[i].order)
      {
	strcpy(value,hg19Chrom[i].value);
	return ;
      }
  strcpy(value,hg19Chrom[25].value);
}

void rankChromosomes(char ** chroms,int *length, int *out)
{
  int i;
  char string[256];
  for(i=0;i<*length;i++)
    {
      sscanf(chroms[i],"%s",string);
      out[i]=getChromosome(string);
    }
}

void valueChromosomes(char ** chroms,int *length, int *out)
{
  int i;
  for(i=0;i<*length;i++)
    {
      getChromosomeValue(out[i],chroms[i]);
    }
}

int compare(int start1, int end1,int start2, int end2)
{
  int width = (end1-start1)*2;
  int difference = (start1+end1-start2-end2);
  if(width-difference<0)
    return 0 ; //move reads forward
  else if(width+difference<0)
    return 1 ; //move peaks forward
  else
    return 2 ; //add to score and move reads forward
}

void pileup(char ** filename,int * chro,int *start,
	    int *end,int *peaknum,int *scores)
{
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char string[1024];
  int inStart,inEnd,inChro;
  int i=0;
  int decision;
  while(fgets(buffer,1024, f))
    {      
      sscanf(buffer,"chr%s\t%d\t%d", string,&inStart,&inEnd);
      inChro=getChromosomeShort(string);
      if(inChro==chro[i])
	{
	  decision=compare(start[i],end[i],inStart,inEnd);
	    if(decision==1)
	      {
	      while(1)
		{
		  i=i+1;
		  if(i>=*peaknum-1) return;
		  if(inChro!=chro[i]){break;}
		  decision=compare(start[i],end[i],inStart,inEnd);
		  if(decision!=1) break;
		}	      		    
	      }
	    if(decision==2) scores[i]++;
	}
      else if(inChro>chro[i])
	{	  
	while(1)
	  {
	    i=i+1;	    
	    if( i>=*peaknum-1){ return;}
	    if(inChro<=chro[i]){ 
	      break;}
	  }
	}
    }
  fclose(f);
}

void file_length(char ** filename,int * i)
{
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  while(fgets(buffer,1024, f))
    {      
      *i=*i+1;
    }
}

void read_bed(char ** filename,char ** chrom,int *start, int *end)
{
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char  *s;
  char string[256];
  int t;
  int i=0;
  while(fgets(buffer,1024, f))
    {      
      s = strtok (buffer,"\n");
      sscanf(s,"%s\t%d\t%d",string,&start[i],&end[i]);
      t=getChromosome(string);
      getChromosomeValue(t,chrom[i]);
      i=i+1;
    }
  fclose(f);
}


void R_init_myLib(DllInfo *info)
{
R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
