//
//  main.c
//  Genomics Project
//
//  Created by Ahmed  Elshetany  on 11/3/19.
//  Copyright © 2019 Ahmed  Elshetany . All rights reserved.
//

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "genomicsIO.h"
typedef char* string;
const string alu = "CGAGACCACGGTGAAACCCCGTC";


int aluSearch(string str, int numberOfCharsToRead,  int * pPos, int indexArray[40]);
void binCalculator(int *bin1, int *bin2, int *bin3, int *bin4, int *bin5, int w, int r , int * indexArray, int c);
void processChromosome(string fileName, int fileSize, int* numberOfAluPerChro ,int* indexArray);
int getMedianIndex(int* indexArray, int c);
void reset(int * i, int * numberOfAluPerChro, int * indexArray );


int i = 0;
int numberOfAluPerChro = 0;
int indexArray[40] = {0};


int main()
{
    
    int sizeOfChr1 = 230481012;
    int sizeOfChr2 = 240548228;
    int sizeOfChr3 = 198100135;
    int sizeOfChr4 = 189752667;
    int sizeOfChr5 = 181265378;
    int sizeOfChr6 = 170078522;
    int sizeOfChr7 = 158970131;
    int sizeOfChr8 = 144768136;
    int sizeOfChr9 = 121790550;
    int sizeOfChr10= 133262962;
    int sizeOfChr11= 134533742;
    int sizeOfChr12= 133137816;
    int sizeOfChr13=  97983125;
    int sizeOfChr14=  90568149;
    int sizeOfChr15=  84641325;
    int sizeOfChr16=  81805943;
    int sizeOfChr17=  82920204;
    int sizeOfChr18=  80089605;
    int sizeOfChr19=  58440758;
    int sizeOfChr20=  63944257;
    int sizeOfChr21=  40088619;
    int sizeOfChr22=  39159777;
    int sizeOfChrX = 154893029;
    int sizeOfChrY =  26415043;
    
    
    
    InitializeGenomeIO();
    
    processChromosome("chr1.fa", sizeOfChr1, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr2.fa", sizeOfChr2, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr3.fa", sizeOfChr3, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr4.fa", sizeOfChr4, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr5.fa", sizeOfChr5, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr6.fa", sizeOfChr6, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr7.fa", sizeOfChr7, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr8.fa", sizeOfChr8, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr9.fa", sizeOfChr9, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr10.fa", sizeOfChr10, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr11.fa", sizeOfChr11, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr12.fa", sizeOfChr12, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr13.fa", sizeOfChr13, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr14.fa", sizeOfChr14, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr15.fa", sizeOfChr15, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr16.fa", sizeOfChr16, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr17.fa", sizeOfChr17, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr18.fa", sizeOfChr18, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr19.fa", sizeOfChr19, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr20.fa", sizeOfChr20, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr21.fa", sizeOfChr21, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chr22.fa", sizeOfChr22, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chrX.fa", sizeOfChrX, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    processChromosome("chrY.fa", sizeOfChrY, &numberOfAluPerChro,indexArray);
    reset(&i, &numberOfAluPerChro, indexArray );
    
    
    CleanUpGenomeIO();
    return 0;
}

int aluSearch(string dna, int numberOfCharsToRead, int * pPos, int indexArray[40]) {
    
    int count = 0;
    const char *tmp = dna;
    while((tmp = strstr(tmp, alu)))
    {
        if(tmp!=NULL){
            *pPos = (int) tmp - (int) dna + GetTotalDnaRead() - numberOfCharsToRead;
            printf("alu found @ index =  %d \n",*pPos );
            indexArray[i] = *pPos;
            i++;
        }
       count++;
       tmp++;
    }
    
    return count;
}

void binCalculator(int *bin1, int *bin2, int *bin3, int *bin4, int *bin5, int w, int r , int * indexArray, int c) {
    for(int i = 0; i < c ; i ++) {
        if (indexArray[i] >= 0 && indexArray[i] <= w - 1){
            (*bin1) ++;
        }
        if (indexArray[i] >= w && indexArray[i] <= 2*w - 1){
            (*bin2) ++;
        }
        if (indexArray[i] >= 2*w && indexArray[i] <= 3*w - 1){
            (*bin3) ++;
        }
        if (indexArray[i] >= 3*w && indexArray[i] <= 4*w - 1){
            (*bin4) ++;
        }
        if (indexArray[i] >= 4*w && indexArray[i] <= 5*w - 1){
            (*bin5) ++;
        }
    }
    
}

int getMedianIndex(int* indexArray, int c) {
    
    int n = (c + 1) / 2 - 1;
    return indexArray[n];
}

void processChromosome(string fileName, int fileSize, int* numberOfAluPerChro, int* indexArray) {
    int range = fileSize;
    int width = range / 5;
    int bin1 = 0;
    int bin2 = 0;
    int bin3 = 0;
    int bin4 = 0;
    int bin5 = 0;
    int numberOfCharsToRead = 2000000;
    * numberOfAluPerChro = 0;
    int DnaActuallyRead;
    int* pDnaActuallyRead = &DnaActuallyRead;
    *pDnaActuallyRead = 0;
    OpenChromosomeFile(fileName);
    while(GetTotalDnaRead() < fileSize ){
        string dna = (string) calloc (numberOfCharsToRead,sizeof(char));
        dna = FetchDna(numberOfCharsToRead, &DnaActuallyRead);
        int position;
        *numberOfAluPerChro = (*numberOfAluPerChro) + aluSearch(dna, numberOfCharsToRead, &position, indexArray);
        CleanUpFetchedDna(dna);
    }
    printf("number of times found : %d \n", *numberOfAluPerChro);
    int medianIndex = getMedianIndex(indexArray, *numberOfAluPerChro);
    printf("Median Index = %d \n",medianIndex);
    binCalculator(&bin1, &bin2, &bin3, &bin4, &bin5, width, range, indexArray, *numberOfAluPerChro );
    printf(" bin1 = %d \n bin2 = %d \n bin3 = %d \n bin4 = %d \n bin5 = %d \n  ", bin1 , bin2, bin3, bin4, bin5 );
    
}

void reset(int * i, int * numberOfAluPerChro, int * indexArray ){
    *i = 0;
    *numberOfAluPerChro = 0;
    memset(indexArray, 0, 160);
    
    
}
