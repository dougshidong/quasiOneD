#ifndef output_h
#define output_h
#include <string>
template <typename vec>

void outVec(std::string fn2, std::string ra, vec data)
{
    // ra = "r" for read
    // ra = "a" for append
    FILE *Results;
    std::string fn = fname() + fn2;
    Results = fopen(fn.c_str(), ra.c_str());
    for(int i = 0; i < data.size(); i++)
        fprintf(Results, "%.15f\n", data[i]);

    fclose(Results);
}
#endif
