// @auther: Muhammad Saiel
#pragma once

#include <string>
#include <cmath>

// forward declares (avoid ROOT headers in .h)
class TFile;
class TGraph2D;
class TRandom3;

class DCHXT2DLUT {
public:
  DCHXT2DLUT() = default;
  ~DCHXT2DLUT();

  // Load xt_mean and xt_error from ROOT file (par.root)
  bool load(const std::string& rootFile,
            const std::string& meanName  = "xt_mean",
            const std::string& sigmaName = "xt_error");

  bool isLoaded() const { return m_loaded; }

  // Query mean and sigma at (x,y) in cm; outputs in ns
  bool meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const;

  // Sample drift time using provided RNG
  double sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const;

private:
  bool m_loaded {false};

  TFile*    m_file   {nullptr};
  TGraph2D* m_gMean  {nullptr}; // xt_mean
  TGraph2D* m_gSigma {nullptr}; // xt_error

  // small helper: recognize "bad outside-hull" results from TGraph2D
  static bool isInvalid(double mu, double sig, double x, double y);
};

