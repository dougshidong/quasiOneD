#ifndef output_h
#define output_h
#include <string>
#include <stdio.h>
#include "globals.h"
template <typename vec>

void outVec(std::string fn2, std::string ra, vec data)
{
    // wa = "w" for read
    // wa = "a" for append
    FILE *Results;
    std::string fn = "./Results/temp" + fn2;
    Results = fopen(fn.c_str(), ra.c_str());
    for(int i = 0; i < (int)data.size(); i++)
    {
        fprintf(Results, "%.15f\n", data[i]);
    }

    fclose(Results);
}
#endif
