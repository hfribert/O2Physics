// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file LambdaHistManager.h
/// \brief histogram manager for lambda histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCEHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCEHISTMANAGER_H_

#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femtounited
{
namespace twotrackresonancehistmanager
{
// enum for track histograms
enum TwoTrackResonanceHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsMass,
  kTwoTrackResonanceHistLast
};

template <const char* Prefix>
struct ConfTwoTrackResonanceBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{200, 1.0, 1.2}}, "Mass"};
};

constexpr const char PrefixTwoTrackResonanceBinning1[] = "TwoTrackResonanceBinning1";
constexpr const char PrefixTwoTrackResonanceBinning2[] = "TwoTrackResonanceBinning2";
constexpr const char PrefixTwoTrackResonanceBinning3[] = "TwoTrackResonanceBinning3";
using ConfTwoTrackResonanceBinning1 = ConfTwoTrackResonanceBinning<PrefixTwoTrackResonanceBinning1>;
using ConfTwoTrackResonanceBinning2 = ConfTwoTrackResonanceBinning<PrefixTwoTrackResonanceBinning2>;
using ConfTwoTrackResonanceBinning3 = ConfTwoTrackResonanceBinning<PrefixTwoTrackResonanceBinning3>;

constexpr std::array<histmanager::HistInfo<TwoTrackResonanceHist>, kTwoTrackResonanceHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hResonanceMass", "Invariant mass resonance; m (GeV/#it{c}^{2}; Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi;p_{T} (GeV/#it{c});#varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsMass, o2::framework::kTH2F, "hPtVsMass", "p_{T} vs invariant mass; p_{T} (GeV/#it{c}); M (GeV/#it{c}^{2}"}}};

template <typename T>
std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>> makeTwoTrackResonanceHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}}};
};

template <typename T>
auto makeTwoTrackResonanceQaHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsMass, {confBinningAnalysis.pt, confBinningAnalysis.mass}}};
};

constexpr char PrefixTwoTrackResonanceQa[] = "TwoTrackResonanceQA/";
constexpr char PrefixRho1[] = "Rho1/";
constexpr char PrefixPhi1[] = "Phi1/";
constexpr char PrefixKstar1[] = "Kstar1/";
constexpr char PrefixAntiKstar1[] = "AntiKstar1/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class TwoTrackResonanceHistManager
{
 public:
  /// Destructor
  virtual ~TwoTrackResonanceHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {Specs[kMass]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(prefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsMass, HistTable), GetHistDesc(kPtVsMass, HistTable), GetHistType(kPtVsMass, HistTable), {Specs[kPtVsMass]});
    }
  }

  template <modes::Mode mode, typename T>
  void fill(T const& resonance)
  {

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), resonance.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), resonance.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), resonance.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), resonance.resonanceMass());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), resonance.pt(), resonance.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), resonance.pt(), resonance.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), resonance.phi(), resonance.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsMass, HistTable)), resonance.pt(), resonance.resonanceMass());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace twotrackresonancehistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCEHISTMANAGER_H_
