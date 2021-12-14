//
//  genomicsIO.h
//  Genomics Project
//
//  Created by Ahmed  Elshetany  on 11/3/19.
//  Copyright © 2019 Ahmed  Elshetany . All rights reserved.
//

#ifndef genomicsIO_h
#define genomicsIO_h

#pragma once

//Delete this line if you don't want to see INFO statements
#define PRINT_INFO

//Startup for IO system
void InitializeGenomeIO();

//Give the file name as a string " "
//If you place the files in the same directory as your source code you will
//not have to include any additional path information.
void OpenChromosomeFile(char* sFileName);

//Only 2MB can be read in one call.
//INPUT: Pass the number of DNA characters you want to read.
//OUTPUT: pDnaActuallyRead will return the number of DNA characters actually read.
//        dnaToRead != pDnaActuallyRead when the file ends or an error occurs.
//RETURN: Returns a pointer to the string this function creates.
//        Pointer is NULL if there was an error.
char* FetchDna(unsigned int dnaToRead, int* pDnaActuallyRead);

//FetchDna creates memory for you. Clean it up with this when you no longer need it.
void CleanUpFetchedDna(char* sGeneData);

//A count of every DNA byte read. This does not include the 'N' and 'n' characters,
//tabs or newline characters.
int GetTotalDnaRead();

//When you are completely done make sure the memory is managed well by calling this.
void CleanUpGenomeIO();

/*Normal use of these functions would be:
   1) Decide how large of a string you will extract, 1K? 2K? 100K?
   2) For each chromosome file, 1 to 22, X and Y....
      a) Open a chromosome file with OpenChromosomeFile function.
      b) Call FetchData to return a string of the size you chose above
      c) Search that string to see if it contains the ALU gene.
         i) NOTE: If the string you are extracting is large enough it may appear multiple times
            within the same string.
         ii) If found add this data to your statistics.
         iii) Call CleanUpFetchedDNa() with the string so that it can be released into memory.
      d) Repeat until the number of DNA characters actually read does not match the requested
             amount. This means the file has ended (or an error has occurred).
      e) Output the statistics for that chromosome
*/



#endif /* genomicsIO_h */
