//
//  genomicsIO.c
//  Genomics Project
//
//  Created by Ahmed  Elshetany  on 11/3/19.
//  Copyright © 2019 Ahmed  Elshetany . All rights reserved.
//

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
# include <ctype.h>
#include "genomicsIO.h"


#define POINTER_HISTORY_SIZE 1000
#define EXPECTED_LINE_SIZE 52
#define MB_PRINT_AT (30 * 1000000)

#define IS_GENE(x) x == 'a' || x == 'A' || \
x == 'g' ||  x == 'G' || \
x == 'c' || x == 'C' || \
x == 't' || x == 'T'

typedef struct _PointerData
{
    int pointer;
    int allocated;
} PointerData;

FILE* fp = NULL;
int bytesRead = 0;
PointerData pointerHistory[POINTER_HISTORY_SIZE];
int bPointerCycled = 0;
int lastPointerIndex = 0;
char chrName[10];
int chrNum=0;
int gCh;
int lineCnt = 0;
int totalByteCnt;

int SkipJunkData(int *pLineCnt);

void InitializeGenomeIO()
{
    memset(pointerHistory, 0, sizeof(pointerHistory));
    bytesRead = 0;
    bPointerCycled = 0;
}

void OpenChromosomeFile(char* sFileName)
{
    char tmp[10];
    if (fp != NULL)
    {
        printf("ERROR: Chromosome file already open\n");
        return;
    }
    if (NULL == sFileName || sFileName[0] == '\0')
    {
        printf("ERROR: File name was empty\n");
        return;
    }
    if (sFileName[strlen(sFileName) - 2] != 'f' || sFileName[strlen(sFileName) - 1] != 'a')
    {
        printf("ERROR: File should have a .fa extension after being unzipped\n");
        return;
    }
    if (NULL == (fp = fopen(sFileName, "r")))
    {
        printf("ERROR: Could not open file %s. Check path/file location and try again\n",sFileName);
        return;
    }
    fgets(tmp, 10, fp);
    snprintf(chrName, 10, "%s", &tmp[1]);
    chrNum = strtol(&tmp[4], NULL, 10);
    totalByteCnt = 0;
    bytesRead = 0;
    lineCnt = 0;
    if (-1 == (gCh = SkipJunkData(&lineCnt)))
    {
        printf("ERROR: Closing %s\n", sFileName);
        fclose(fp);
        fp = NULL;
        return;
    }
#ifdef PRINT_INFO
    printf("INFO: Ready to process Chromosome %d\n", chrNum);
#endif
}
char* FetchDna(unsigned int dnaToRead, int* pDnaActuallyRead)
{
    char* pRv;
    
    if (NULL == pDnaActuallyRead)
    {
        printf("ERROR: pBytesActuallyRead parameter needs to be a valid pointer\n");
        return NULL;
    }
    if (dnaToRead > 2000000)
    {
        printf("ERROR: Will not read more than 2MB at a time (You entered %d)\n", dnaToRead);
        return NULL;
    }
    if (NULL == (pRv = (char *) malloc(dnaToRead +1)))
    {
        printf("ERROR: Failed to allocate space for %d bytes\n", dnaToRead);
        return NULL;
    }
    
    *pDnaActuallyRead = 0;

    if (pointerHistory[lastPointerIndex].allocated == 1)
    {
        printf("WARNING: Overwriting allocated memory. Have you been freeing data?\n");
        free((void *) pointerHistory[lastPointerIndex].pointer);
    }
    pointerHistory[lastPointerIndex].pointer = (int)pRv;
    pointerHistory[lastPointerIndex].allocated = 1;
    lastPointerIndex++;
    if (lastPointerIndex >= POINTER_HISTORY_SIZE)
    {
        lastPointerIndex = 0;
        bPointerCycled = 1;
    }

    while (*pDnaActuallyRead != dnaToRead)
    {
        int ch;
        lineCnt++;
        if ((totalByteCnt++) % MB_PRINT_AT == 0)
        {
#ifdef PRINT_INFO
            printf("INFO: %d MB read\n", totalByteCnt/1000000);
#endif
        }
        if (gCh != -1)
        {
            ch = gCh;
            gCh = -1;
        }
        else if (EOF == (ch = fgetc(fp)))
        {
            printf("WARNING: Unexpected end of %s\n", chrName);
            printf("WARNING: File may be corrupted. Did you unzip it?\n");
            fclose(fp);
            fp = NULL;
            return NULL;
        }
        else if (ch == '\n' && lineCnt != (EXPECTED_LINE_SIZE - 1))
        {
            printf("ERROR: File did not match expected format. It may be corrupted\n");
            return NULL;
        }
        else if (ch == '\n')
        {
            lineCnt = 0;
            continue;
        }
        else if (ch == '\0' && lineCnt != EXPECTED_LINE_SIZE)
        {
            printf("ERROR: File did not match expected format. It may be corrupted\n");
            return NULL;
        }
        else if(ch == '\0')
        {
            lineCnt = 0;
            continue;
        }

        if (IS_GENE(ch))
        {
            bytesRead++;
            pRv[(*pDnaActuallyRead)] = toupper(ch);
            (*pDnaActuallyRead)++;
            if (*pDnaActuallyRead == (dnaToRead))
            {
                pRv[*pDnaActuallyRead] = '\0';
                if (lineCnt == (EXPECTED_LINE_SIZE - 2))
                {
                    fgetc(fp);
                    if ((totalByteCnt++) % MB_PRINT_AT == 0)
                    {
#ifdef PRINT_INFO
                        printf("INFO: %d MB read\n", totalByteCnt / 1000000);
#endif
                    }
                    lineCnt = 0;
                }
                break;
            }
        }
        else if (ch == 'N' || ch == 'n')
        {
            int ch2;
            if ((ch2 = SkipJunkData(&lineCnt)) < 0)
            {
                pRv[*pDnaActuallyRead] = '\0';
                return pRv;
            }
        }
        else
        {
            printf("ERROR: File did not match expected format. It may be corrupted\n");
            return NULL;
        }
    }
    return pRv;
}
void CleanUpFetchedDna(char* sGeneData)
{
    int endPt = lastPointerIndex;
    if (NULL == sGeneData)
    {
        printf("ERROR: Cannot clean up a NULL pointer\n");
        return;
    }
    if (bPointerCycled)
    {
        endPt = POINTER_HISTORY_SIZE;
    }
    for (int i = 0; i < endPt; i++)
    {
        if (pointerHistory[i].pointer == (int)sGeneData)
        {
            pointerHistory[i].allocated = 0;
            pointerHistory[i].pointer = NULL;
            free(sGeneData);
            return;
        }
    }

    printf("WARNING: Pointer address was not found in history...This might cause code to crash\n");
    free(sGeneData);
    return;
}

int GetTotalDnaRead()
{
    return bytesRead;
}

int SkipJunkData(int *pLineCnt)
{
    int i;
    int lineCnt;
    int ch;

    while (!feof(fp))
    {
        (*pLineCnt)++;
        if ((totalByteCnt++) % MB_PRINT_AT == 0)
        {
#ifdef PRINT_INFO
            printf("INFO: %d MB read\n", totalByteCnt / 1000000);
#endif
        }
        if (EOF == (ch = fgetc(fp)))
        {
            break;
        }
        if (ch != 'N' && ch != 'n' && ch != '\n' && ch != '\0')
        {
            if (IS_GENE(ch))
            {
                gCh = ch;
                (*pLineCnt)--;
                return ch;
            }
            printf("ERROR: File did not match expected format. It may be corrupted\n");
            return -10;
        }
        else if (ch == '\n')
        {
            *pLineCnt = 0;
        }
    }
#ifdef PRINT_INFO
    printf("INFO: Closing %sINFO: %d total DNA read\n",chrName,bytesRead);
#endif
    fclose(fp);
    fp = NULL;
    return EOF;
}

void CleanUpGenomeIO()
{
    int endPt = lastPointerIndex;
    int cnt = 0;
    if (bPointerCycled)
    {
        endPt = POINTER_HISTORY_SIZE;
    }
    for (int i = 0; i < endPt; i++)
    {
        if (pointerHistory[i].allocated == 1 && NULL != pointerHistory[i].pointer)
        {
            free((void *)pointerHistory[i].pointer);
            pointerHistory[i].allocated = 0;
            cnt++;
        }
    }
    if (cnt > 0)
    {
        printf("SCOLDING: Had to clean up %d pointers after you!", cnt);
        return;
    }
}

