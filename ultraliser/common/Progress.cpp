#include "Progress.h"

void printProgressBar(const size_t& current,
                      const size_t& total,
                      const size_t barLength)
{
    float percentage = ((100.f * current) / (1.f * total));
    float starts = std::floor((percentage * barLength) / 100.f);
    float spaces = std::floor(barLength - starts);
    // std::string bar = "* │";
    std::string bar = "* Progress │";
    for(int _i_ = 0; _i_ < int(starts); _i_++) bar += "▒";
    for(int _i_ = 0; _i_ < int(spaces); _i_++) bar += " "; bar+= "│";
    printf("\r\t%s (%2.2f %%)", bar.c_str(), percentage);
    fflush(stdout);
}

void printFractionProgressBar(const size_t& current,
                              const size_t& total,
                              const size_t barLength)
{
    float percentage = ((100.f * current) / (1.f * total));
    if (static_cast< size_t >(percentage) % 10 == 0)
    {
        float starts = std::floor((percentage * barLength) / 100.f);
        float spaces = std::floor(barLength - starts);
        std::string bar = "* Progress │";
        for(int i = 0; i < int(starts); i++) bar += "▒";
        for(int i = 0; i < int(spaces); i++) bar += " "; bar+= "│";
        printf("\r\t%s (%2.2f %%)", bar.c_str(), percentage);
        fflush(stdout);
    }
}

void progressUpdate(size_t& progressValue)
{
#ifdef ULTRALISER_USE_OPENMP
#pragma omp atomic
#endif
    ++progressValue;
}
