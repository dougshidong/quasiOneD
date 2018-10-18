#ifndef OUTPUT_H
#define OUTPUT_H
#include <string>
#include <stdio.h>
template <typename vec>

void outVec(std::string case_name, std::string filename, std::string ra, vec data) {
    // wa = "w" for read
    // wa = "a" for append
    FILE *Results;
    std::string fn = "./Results/"+ case_name + filename;
    Results = fopen(fn.c_str(), ra.c_str());
    for (int i = 0; i < (int)data.size(); i++) {
        fprintf(Results, "%1.15e\n", data[i]);
    }

    fclose(Results);
}
#endif
