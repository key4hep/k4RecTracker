// @auther: Muhammad saiel
//
#ifndef DCH_PHYSICS_TOOLS_H
#define DCH_PHYSICS_TOOLS_H
//#include "XTRELTIME.h"
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
//include EDM4hep
#include <edm4hep/MCParticle.h>
#include "TF1.h"

// Function compute betaGamma:
inline double compute_beta_gamma(const edm4hep::MCParticle& mc)
{
  const auto& p = mc.getMomentum();        // GeV
  const double p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  const double m     = mc.getMass();       // GeV
  if (m <= 0.0) return 0.0;
  return p_mag / m;
}

// This function is used to calculate the dNdx using BB.h file:
// when we have a beam-test TGraph gNcldx(bg)->(clusters/cm),
inline double DCHdigi_v03::get_dNcldx_per_cm(double betagamma, const edm4hep::MCParticle& mc) const
{
  // 1) reconstruct momentum from β
  const double m_GeV = mc.getMass();
  const double p_GeV = betagamma * m_GeV;

  // 2) call your header: bethe_bloch(xp, par)
  double xp[1]  = { p_GeV };                           // GeV (header multiplies by 1e3 internally)
  double par[2] = { m_GeV * 1.0e3,                     // mass in MeV
                    m_MeanExcEnergy_eV.value() * 1e-6  // I in MeV
                  };
  const double dEdx_MeVcm2_per_g = BB::bethe_bloch(xp, par);
  const double dEdx_MeV_per_cm = dEdx_MeVcm2_per_g * m_GasDensity_g_cm3.value();
  const double lambda_per_cm = (dEdx_MeV_per_cm * 1.0e6) / m_W_eff_eV.value();

  return (lambda_per_cm > 0.0) ? lambda_per_cm : 0.0;
}

namespace
{

inline int ZT_poisson(double mu, std::mt19937_64& rng)
{
  if (!(mu > 0) || !std::isfinite(mu)) return 0;
  std::poisson_distribution<int> pois(mu);
  int k = 0;
  do {k = pois(rng); }  while (k == 0);
  return k;
}

// Generate exponentially-distributed cluster positions along a step of length l_mm.
// Output: vector of positions in mm, measured from the start of the step.
inline std::vector<double>
generate_cluster_positions_mm(double l_mm, double lambda_per_cm, std::mt19937_64& rng)
{
  std::vector<double> pos_mm;
  if (l_mm <= 0.0 || lambda_per_cm <= 0.0) return pos_mm;

  const double l_cm = 0.1 * l_mm;
  std::exponential_distribution<double> expo(lambda_per_cm);

  double s_cm = 0.0;
  while (true) {
    s_cm += expo(rng);           // draw next spacing
    if (s_cm >= l_cm) break;     // stop if we exceed the step length
    pos_mm.push_back(10.0 * s_cm); // store position in mm
  }
  return pos_mm;

}

// --- Experimental cluster size probabilities w(n) for He-iSobutane (90:10) ---
// Based on Fischle et al., NIM A301 (1991)
inline const std::vector<double>& w_cluster_He_iC4H10()
{
  constexpr int KMAX = 35;
  static std::vector<double> w;
  static bool init = false;

  if (!init) {
    w.assign(KMAX, 0.0);
    w[0]  = 76.60;
    w[1]  = 12.50;
    w[2]  =  4.60;
    w[3]  =  2.00;
    w[4]  =  1.20;
    w[5]  =  0.75;
    w[6]  =  0.50;
    w[7]  =  0.36;
    w[8]  =  0.25;
    w[9]  =  0.19;
    w[10] =  0.14;
    w[11] =  0.10;
    w[12] =  0.08;
    w[13] =  0.06;
    w[14] =  0.048;
    w[15] =  0.043;
    w[16] =  0.038;
    w[17] =  0.034;
    w[18] =  0.030;

    //for electrons more than 19 use the function: 20/k^2
    for (int k = 20; k <= KMAX; ++k) {
      w[k - 1] = 10.9 / (double(k) * double(k));
    }
    for (auto& x : w) x *= 0.01;

    const double sum = std::accumulate(w.begin(), w.end(), 0.0);
    if (sum > 0.0) {
      for (auto& x : w) x /= sum;
    }
    init = true;
  }
  return w;
}

inline double mean_cluster_size_He()
{
  const auto& w = w_cluster_He_iC4H10();
  double Ne_mean = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    Ne_mean += (double(i) + 1.0) * w[i];
  }
  return Ne_mean;
}

// --- Draw a cluster size n_e according to experimental probabilities ---
inline int sample_cluster_size(std::mt19937_64& rng)
{
  const auto& w = w_cluster_He_iC4H10();
  std::discrete_distribution<int> dist(w.begin(), w.end());
  int n = dist(rng) + 1; // +1 because bins start at 1
  return n;
}

} //end namespace

/***************************************************************************
 * This block of code is needed when we are using 1D LUD
 * *************************************************************************
inline double DCHdigi_v03::electronDriftTime(double r_cm, TRandom3& myRandom) const
{
  if (!m_xtHelper) {
    return 0.0;
  }

  // call XTREL::spaceToTime(dist, version) inherited by XTRELTIME:
  // dist: distance in cm; version=0 uses table with interpolation
  Float_t* x2time = m_xtHelper->spaceToTime(static_cast<Float_t>(r_cm), 0);

  // x2time[0] = mean drift time
  // x2time[1] = longitudinal diffusion sigma (time units)
  const double t0    = static_cast<double>(x2time[0]);
  const double sigma = static_cast<double>(x2time[1]);

  // smear with Gaussian diffusion
  const double t = (t0 + myRandom.Gaus(0.0, sigma))*1000.0;
  const double t_ns = std::max(0.0, t);

  return t_ns;
}
****************************************************************************/

// This block is sampling avalache charge for 
// one electron using polya distibution
inline double DCHdigi_v03::avalancheCharge(TRandom3& myRandom) const
{
	if (!m_polya) return 0.0;

	double q = m_polya->GetRandom(&myRandom);
	return (q > 0.0) ? q : 0.0;

}

// ----------------------------------------------------------------------
//  Single-electron pulse shape (Point 7: rise and fall times)
// ----------------------------------------------------------------------
// t_ns   : observation time (ns)
// t0_ns  : electron arrival time (ns)
// q      : avalanche charge (arbitrary units, e.g. from Polya)
//
// Shape: double exponential CR-RC-like pulse
// ----------------------------------------------------------------------
inline double DCHdigi_v03::singleElectronPulse(double t_ns,
                                               double t0_ns,
                                               double q) const
{
  const double dt = t_ns - t0_ns;
  if (dt <= 0.0) {
    // No signal before the electron arrival
    return 0.0;
  }

  const double tau_r = m_pulseRiseTime_ns.value();  // ns
  const double tau_f = m_pulseFallTime_ns.value();  // ns

  // Safety: if parameters are not set, just return a delta-like pulse
  if (tau_r <= 0.0 || tau_f <= 0.0 || tau_r == tau_f) {
    return q;
  }

  // Double exponential shape (not normalized, we keep q as overall scale)
  const double term_f = std::exp(-dt / tau_f);
  const double term_r = std::exp(-dt / tau_r);

  // Optional normalization factor so peak ~ q (roughly)
  const double norm = 1.0 / (tau_f - tau_r);

  const double q_eff = q * m_pulseAmplitudeScale.value();

  return q_eff * norm * (term_f - term_r);
}
#endif
