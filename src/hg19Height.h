/* This file is part of peakAnalysis,
   http://github.com/alexjgriffith/alpha-score/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int getChromosome(char value[10]);
int getChromosomeShort(char value[2]);
void getChromosomeValue(int order,char * value);
void rankChromosomes(char ** chroms,int *length, int *out);
void valueChromosomes(char ** chroms,int *length, int *out);
int compare(int start1, int end1,int start2, int end2);
void pileup(char ** filename,int * chro,int *start,
	    int *end,int *peaknum,int *scores);
void file_length(char ** filename,int * i);
void read_bed(char ** filename,char ** chrom,int *start, int *end);

R_CMethodDef cMethods[]={
{"read_bed",(DL_FUNC) &read_bed,4,{STRSXP,STRSXP,INTSXP,INTSXP}},
{"file_length",(DL_FUNC) &file_length,2,{STRSXP,INTSXP} },
{"pileup",(DL_FUNC) &pileup,6, {STRSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP}},
{"rankChromosomes",(DL_FUNC) &rankChromosomes,6,{STRSXP,INTSXP,INTSXP}},
{"valueChromosomes",(DL_FUNC) &valueChromosomes,6,{STRSXP,INTSXP,INTSXP}},
{NULL,NULL,0}
  };

void R_init_myLib(DllInfo *info);


typedef struct chromosomes {
  int order;
  char value[10];
  char shortValue[2];} chromosomes;

static const chromosomes hg19Chrom[] = {
{0,"chr1","1"},
{1,"chr2","2"},
{2,"chr3","3"},
{3,"chr4","4"},
{4,"chr5","5"},
{5,"chr6","6"},
{6,"chr7","7"},
{7,"chrY","Y"},
{8,"chr8","8"},
{9,"chr9","9"},
{10,"chr10","10"},
{11,"chr11","11"},
{12,"chr12","12"},
{13,"chr13","13"},
{14,"chr14","14"},
{15,"chr15","15"},
{16,"chr16","16"},
{17,"chr17","17"},
{18,"chr18","18"},
{19,"chr19","19"},
{20,"chr20","20"},
{21,"chrX","X"},
{22,"chr21","21"},
{23,"chr22","22"},
{24,"chrM","M"},
{25,"chr",""}};

