// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file kinkHistManager.h
/// \brief histogram manager for kink histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace kinkhistmanager
{
// enum for kink histograms
enum KinkHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kKinkAngle,
  kDcaMothToPV,
  kDcaDaugToPV,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsKinkAngle,
  kPtVsDecayRadius,
  kKinkHistLast
};

#define KINK_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                         \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

template <const char* Prefix>
struct ConfSigmaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_BINNING(1.1, 1.3)
};
#undef KINK_DEFAULT_BINNING

constexpr const char PrefixSigmaBinning1[] = "SigmaBinning1";
using ConfSigmaBinning1 = ConfSigmaBinning<PrefixSigmaBinning1>;

template <const char* Prefix>
struct ConfKinkQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis kinkAngle{"kinkAngle", {{100, 0, 3.15}}, "Kink Angle (rad)"};
  o2::framework::ConfigurableAxis dcaMothToPV{"dcaMothToPV", {{150, 0, 1.5}}, "Mother DCA to PV (cm)"};
  o2::framework::ConfigurableAxis dcaDaugToPV{"dcaDaugToPV", {{150, 0, 1.5}}, "Daughter DCA to PV (cm)"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex position (cm)"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius (cm)"};
};

constexpr const char PrefixSigmaQaBinning1[] = "SigmaQaBinning1";
using ConfSigmaQaBinning1 = ConfKinkQaBinning<PrefixSigmaQaBinning1>;

// must be in sync with enum KinkHist
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<KinkHist>, kKinkHistLast> KinkHistTable = {{
  {kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
  {kEta, o2::framework::kTH1F, "hEta", "Pseudorapidity; #eta; Entries"},
  {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
  {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
  {kSign, o2::framework::kTH1F, "hSign", "Sign; sign; Entries"},
  {kKinkAngle, o2::framework::kTH1F, "hKinkAngle", "Kink Angle; Angle (rad); Entries"},
  {kDcaMothToPV, o2::framework::kTH1F, "hDcaMothToPV", "Mother DCA to PV; DCA (cm); Entries"},
  {kDcaDaugToPV, o2::framework::kTH1F, "hDcaDaugToPV", "Daughter DCA to PV; DCA (cm); Entries"},
  {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "Decay Vertex X; x (cm); Entries"},
  {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Decay Vertex Y; y (cm); Entries"},
  {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Decay Vertex Z; z (cm); Entries"},
  {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Decay Distance from PV; r (cm); Entries"},
  {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse Decay Radius; r_{xy} (cm); Entries"},
  {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}); #eta"},
  {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}); #varphi"},
  {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi; #eta"},
  {kPtVsKinkAngle, o2::framework::kTH2F, "hPtVsKinkAngle", "p_{T} vs kink angle; p_{T} (GeV/#it{c}); kink angle (rad)"},
  {kPtVsDecayRadius, o2::framework::kTH2F, "hPtVsDecayRadius", "p_{T} vs transverse decay radius; p_{T} (GeV/#it{c}); r_{xy} (cm)"}
}};

template <typename T>
auto makeKinkAnalysisHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<KinkHist, std::vector<framework::AxisSpec>> makeKinkQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kKinkAngle, {confBinningQa.kinkAngle}},
    {kDcaMothToPV, {confBinningQa.dcaMothToPV}},
    {kDcaDaugToPV, {confBinningQa.dcaDaugToPV}},
    {kDecayVtxX, {confBinningQa.decayVertex}},
    {kDecayVtxY, {confBinningQa.decayVertex}},
    {kDecayVtxZ, {confBinningQa.decayVertex}},
    {kDecayVtx, {confBinningQa.decayVertex}},
    {kTransRadius, {confBinningQa.transRadius}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsKinkAngle, {confBinningAnalysis.pt, confBinningQa.kinkAngle}},
    {kPtVsDecayRadius, {confBinningAnalysis.pt, confBinningQa.transRadius}}};
}

constexpr char PrefixSigmaQa[] = "SigmaQA/";
constexpr char PrefixSigma1[] = "Sigma1/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class KinkHistManager
/// \brief Class for histogramming kink properties
template <const char* prefix, modes::Mode mode, modes::Kink kink>
class KinkHistManager
{
 public:
  /// Destructor
  virtual ~KinkHistManager() = default;
  
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  /// \param Specs Map of histogram specifications
  void init(o2::framework::HistogramRegistry* registry, std::map<KinkHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, KinkHistTable), GetHistDesc(kPt, KinkHistTable), GetHistType(kPt, KinkHistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, KinkHistTable), GetHistDesc(kEta, KinkHistTable), GetHistType(kEta, KinkHistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, KinkHistTable), GetHistDesc(kPhi, KinkHistTable), GetHistType(kPhi, KinkHistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, KinkHistTable), GetHistDesc(kMass, KinkHistTable), GetHistType(kMass, KinkHistTable), {Specs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSign, KinkHistTable), GetHistDesc(kSign, KinkHistTable), GetHistType(kSign, KinkHistTable), {Specs[kSign]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(prefix) + std::string(QaDir);

      // Basic kinematic histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPt, KinkHistTable), GetHistDesc(kPt, KinkHistTable), GetHistType(kPt, KinkHistTable), {Specs[kPt]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kEta, KinkHistTable), GetHistDesc(kEta, KinkHistTable), GetHistType(kEta, KinkHistTable), {Specs[kEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhi, KinkHistTable), GetHistDesc(kPhi, KinkHistTable), GetHistType(kPhi, KinkHistTable), {Specs[kPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMass, KinkHistTable), GetHistDesc(kMass, KinkHistTable), GetHistType(kMass, KinkHistTable), {Specs[kMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kSign, KinkHistTable), GetHistDesc(kSign, KinkHistTable), GetHistType(kSign, KinkHistTable), {Specs[kSign]});

      // Kink-specific QA histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kKinkAngle, KinkHistTable), GetHistDesc(kKinkAngle, KinkHistTable), GetHistType(kKinkAngle, KinkHistTable), {Specs[kKinkAngle]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDcaMothToPV, KinkHistTable), GetHistDesc(kDcaMothToPV, KinkHistTable), GetHistType(kDcaMothToPV, KinkHistTable), {Specs[kDcaMothToPV]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDcaDaugToPV, KinkHistTable), GetHistDesc(kDcaDaugToPV, KinkHistTable), GetHistType(kDcaDaugToPV, KinkHistTable), {Specs[kDcaDaugToPV]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxX, KinkHistTable), GetHistDesc(kDecayVtxX, KinkHistTable), GetHistType(kDecayVtxX, KinkHistTable), {Specs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxY, KinkHistTable), GetHistDesc(kDecayVtxY, KinkHistTable), GetHistType(kDecayVtxY, KinkHistTable), {Specs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxZ, KinkHistTable), GetHistDesc(kDecayVtxZ, KinkHistTable), GetHistType(kDecayVtxZ, KinkHistTable), {Specs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtx, KinkHistTable), GetHistDesc(kDecayVtx, KinkHistTable), GetHistType(kDecayVtx, KinkHistTable), {Specs[kDecayVtx]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTransRadius, KinkHistTable), GetHistDesc(kTransRadius, KinkHistTable), GetHistType(kTransRadius, KinkHistTable), {Specs[kTransRadius]});

      // 2D QA histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, KinkHistTable), GetHistDesc(kPtVsEta, KinkHistTable), GetHistType(kPtVsEta, KinkHistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, KinkHistTable), GetHistDesc(kPtVsPhi, KinkHistTable), GetHistType(kPtVsPhi, KinkHistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, KinkHistTable), GetHistDesc(kPhiVsEta, KinkHistTable), GetHistType(kPhiVsEta, KinkHistTable), {Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsKinkAngle, KinkHistTable), GetHistDesc(kPtVsKinkAngle, KinkHistTable), GetHistType(kPtVsKinkAngle, KinkHistTable), {Specs[kPtVsKinkAngle]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsDecayRadius, KinkHistTable), GetHistDesc(kPtVsDecayRadius, KinkHistTable), GetHistType(kPtVsDecayRadius, KinkHistTable), {Specs[kPtVsDecayRadius]});
    }
  }

  /// Fill histograms for kink candidates
  /// \param kinkcandidate Kink candidate to fill histograms for
  template <typename T>
  void fill(T const& kinkcandidate)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, KinkHistTable)), kinkcandidate.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, KinkHistTable)), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, KinkHistTable)), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, KinkHistTable)), kinkcandidate.mass());

      if constexpr (isEqual(kink, modes::Kink::kSigma)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, KinkHistTable)), kinkcandidate.sign());
      }
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      // Basic kinematic histograms
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPt, KinkHistTable)), kinkcandidate.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kEta, KinkHistTable)), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhi, KinkHistTable)), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMass, KinkHistTable)), kinkcandidate.mass());
      
      if constexpr (isEqual(kink, modes::Kink::kSigma)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kSign, KinkHistTable)), kinkcandidate.sign());
      }

      // Kink-specific QA histograms
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kKinkAngle, KinkHistTable)), kinkcandidate.kinkAngle());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDcaMothToPV, KinkHistTable)), kinkcandidate.dcaMothToPV());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDcaDaugToPV, KinkHistTable)), kinkcandidate.dcaDaugToPV());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, KinkHistTable)), kinkcandidate.decayVtxX());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, KinkHistTable)), kinkcandidate.decayVtxY());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, KinkHistTable)), kinkcandidate.decayVtxZ());
      
      // Calculate decay distance from PV
      float decayDistance = std::sqrt(kinkcandidate.decayVtxX() * kinkcandidate.decayVtxX() + 
                                     kinkcandidate.decayVtxY() * kinkcandidate.decayVtxY() + 
                                     kinkcandidate.decayVtxZ() * kinkcandidate.decayVtxZ());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, KinkHistTable)), decayDistance);
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, KinkHistTable)), kinkcandidate.transRadius());

      // 2D QA histograms
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, KinkHistTable)), kinkcandidate.pt(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, KinkHistTable)), kinkcandidate.pt(), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, KinkHistTable)), kinkcandidate.phi(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsKinkAngle, KinkHistTable)), kinkcandidate.pt(), kinkcandidate.kinkAngle());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsDecayRadius, KinkHistTable)), kinkcandidate.pt(), kinkcandidate.transRadius());
    }
  }

  /// Fill histograms for kink candidates - overload with track table argument
  /// \param kinkcandidate Kink candidate to fill histograms for  
  /// \param tracks Track table (ignored for kinks)
  template <typename T1, typename T2>
  void fill(T1 const& kinkcandidate, T2 const& /*tracks*/)
  {
    fill(kinkcandidate);
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace kinkhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_