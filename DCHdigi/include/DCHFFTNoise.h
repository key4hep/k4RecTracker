// @auther: Muhammad Saiel
#ifndef DCH_FFT_NOISE_H
#define DCH_FFT_NOISE_H

#include <string>
#include <vector>

class TRandom3;

namespace dchfft {

bool loadFFTNoiseTemplate(
    const std::string& rootFile,
    const std::string& treeName,
    std::vector<double>& magTemplate,
    double& ampNorm,
    double& maxFreq,
    int& size,
    std::string* errMsg = nullptr
);

// Generate colored noise samples (length = nSamples) using FFT magnitude template.
// Units: dt_ns is in ns, so frequency is in 1/ns.
std::vector<double> makeFFTNoiseSamples(
    int nSamples,
    double dt_ns,
    const std::vector<double>& magTemplate,
    double maxFreqTemplate,
    TRandom3& rnd,
    bool removeDC
);

} // namespace dchfft

#endif

