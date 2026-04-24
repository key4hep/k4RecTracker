#include "DCHXT2DLUT.h"

#include "TFile.h"
#include "TGraph2D.h"
#include "TRandom3.h"

#include <iostream>

DCHXT2DLUT::~DCHXT2DLUT() {
  // Close and delete ROOT file if we own it
  if (m_file) {
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  m_gMean = nullptr;
  m_gSigma = nullptr;
  m_loaded = false;
}

/// Helper function: checks if the LUT query result looks invalid.
bool DCHXT2DLUT::isInvalid(double mu, double sig, double x, double y) {
  if (!std::isfinite(mu) || !std::isfinite(sig)) return true;
  if (mu < 0.0 || sig < 0.0) return true;

  // Common ROOT behavior: outside hull returns ~0
  const double r = std::hypot(x, y);
  if (r > 1e-4 && std::abs(mu) < 1e-12 && std::abs(sig) < 1e-12) return true;

  return false;
}

/// Load the XT LUT from a ROOT file.
/// - Fetches two TGraph2D objects:
///     meanName  -> "xt_mean"  (mean drift time surface)
///     sigmaName -> "xt_error" (time spread surface)
bool DCHXT2DLUT::load(const std::string& rootFile,
                      const std::string& meanName,
                      const std::string& sigmaName) {
  // cleanup if re-loading
  if (m_file) {
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  m_gMean = nullptr;
  m_gSigma = nullptr;
  m_loaded = false;

  m_file = TFile::Open(rootFile.c_str(), "READ");
  if (!m_file || m_file->IsZombie()) {
    std::cerr << "[DCHXT2DLUT] ERROR: cannot open file " << rootFile << "\n";
    if (m_file) { delete m_file; m_file = nullptr; }
    return false;
  }

  // Retrieve objects from the file:
  m_file->GetObject(meanName.c_str(),  m_gMean);
  m_file->GetObject(sigmaName.c_str(), m_gSigma);

  if (!m_gMean || !m_gSigma) {
    std::cerr << "[DCHXT2DLUT] ERROR: missing TGraph2D " << meanName
              << " or " << sigmaName << " in " << rootFile << "\n";
    m_file->Close();
    delete m_file;
    m_file = nullptr;
    m_gMean = nullptr;
    m_gSigma = nullptr;
    return false;
  }

  // Mark that everything is ready to use
  m_loaded = true;
  std::cout << "[DCHXT2DLUT] Loaded " << meanName << " and " << sigmaName
            << " from " << rootFile << "\n";
  return true;
}

/// Query mean and sigma drift-time at a given (x,y) point.
bool DCHXT2DLUT::meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const {
  if (!m_loaded || !m_gMean || !m_gSigma) {
    mu_ns = 0.0;
    sig_ns = 0.0;
    return false;
  }

  // ROOT interpolates the 2D scattered data surface at (x,y)
  mu_ns  = m_gMean->Interpolate(x_cm, y_cm);
  sig_ns = m_gSigma->Interpolate(x_cm, y_cm);

  // Filter out NaN/inf/negative/outside-hull zero behavior
  if (isInvalid(mu_ns, sig_ns, x_cm, y_cm)) {
    mu_ns = 0.0;
    sig_ns = 0.0;
    return false;
  }
  return true;
}

/// Sample drift time in ns using the digitizer RNG.
/// 1) Get (mu, sigma) from meanSigma()
/// 2) If valid and sigma>0 -> Gaussian smear: t ~ N(mu, sigma)
/// 3) If sigma==0 -> use t = mu
/// 4) Clamp to >=0 (no negative drift times)
double DCHXT2DLUT::sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const {
  double mu = 0.0, sig = 0.0;
  const bool ok = meanSigma(x_cm, y_cm, mu, sig);

  // If invalid/outside, return 0 (simple behavior)
  if (!ok) return 0.0;

  // Apply diffusion/time-smearing
  double t = (sig > 0.0) ? rng.Gaus(mu, sig) : mu;
  if (t < 0.0) t = 0.0;
  return t;
}

