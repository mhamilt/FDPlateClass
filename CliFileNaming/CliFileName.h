//
//  CliFileName.cpp
//  FDPlateClass
//
//  Created by mhamilt7 on 08/05/2018.
//  Copyright © 2018 mhamilt7. All rights reserved.
//

#include <iostream>
#include <string>

/**
  Create a file name from Command Line INterface arguement

 @param argv The argv[] character pointer array from `int main (int argc, const char *argv[])`
 @param defaultFileName For debugging declare a default file path in the user home folder
 @return a pointer to a character array containing the file name
 */
char *CliSetFilename (const char *argv[], const char *defaultFileName)
{
    char *outputfname = nullptr;
    
    if(!argv[1])
    {
        const char *homedir = getenv("HOME");
        if(!homedir)
        {
            printf("Couldn't find Home Directory. Enter a filename\n");
            return nullptr;
        }
        printf("NO FILE NAME ENTERED\nUSING DEFAULT FILE DESTINATION\n");
//        const char *def_fname = "/Downloads/Plate.wav";
        const int length = int(strlen(homedir)) + int(strlen(defaultFileName));
        outputfname = new char[length+1]();
        strncpy(outputfname,homedir, int(strlen(homedir)));
        strcat(outputfname, defaultFileName);
    }
    else
    {
        const int length = int(strlen(argv[1]));
        outputfname = new char[length+1]();
        strncpy(outputfname, argv[1], length);
    }
    
    return outputfname;
}
