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

/// \file v0HistManager.h
/// \brief histogram manager for vzero histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_

#include "PWGCF/FemtoUnited/Core/histManager.h"
#include "PWGCF/FemtoUnited/Core/modes.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femtounited
{
namespace v0histmanager
{
// enum for track histograms
enum V0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kMassLambda,
  kMassAntiLambda,
  kMassK0short,
  kCosPa,
  kDecayDauDca,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsLambdaMass,
  kPtVsAntiLambdaMass,
  kPtVsK0shortMass,
  kLambdaMassVsAntiLambdaMass,
  kK0shortMassVsLambdaMass,
  kK0shortMassVsAntiLambdaMass,
  kV0HistLast
};

#define V0_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                         \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

template <const char* Prefix>
struct ConfLambdaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(1.0, 1.2)
};
template <const char* Prefix>
struct ConfK0shortBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(0.475, 0.515)
};
#undef V0_DEFAULT_BINNING

constexpr const char PrefixLambdaBinning1[] = "LambdaBinning1";
using ConfLambdaBinning1 = ConfLambdaBinning<PrefixLambdaBinning1>;
constexpr const char PrefixK0shortBinning1[] = "K0shortBinning1";
using ConfK0shortBinning1 = ConfK0shortBinning<PrefixK0shortBinning1>;

template <const char* Prefix>
struct ConfV0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis massLambda{"massLambda", {{200, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massAntiLambda{"massAntiLambda", {{100, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massK0short{"massK0short", {{200, 0.45, 0.55}}, "Mass for k0short hypothesis"};
};

constexpr const char PrefixLambdaQaBinning1[] = "LambdaQaBinning1";
using ConfLambdaQaBinning1 = ConfV0QaBinning<PrefixLambdaQaBinning1>;

constexpr const char PrefixK0shortQaBinning1[] = "K0shortQaBinning1";
using ConfK0shortQaBinning1 = ConfV0QaBinning<PrefixK0shortQaBinning1>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<V0Hist>, kV0HistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv}} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kMassLambda, o2::framework::kTH1F, "hMassLambda", "#Lambda mass; m_{p#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kMassAntiLambda, o2::framework::kTH1F, "hMassAntiLambda", "#bar{#Lambda} mass; m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2}); Entries"},
   {kMassK0short, o2::framework::kTH1F, "hMassK0short", "K^{0}_{s} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; coa(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse radius ; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsLambdaMass, o2::framework::kTH2F, "hPtVsLambdaMass", "p_{T} vs #Lambda mass; p_{T} (GeV/#it{c}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kPtVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiLambdaMass", "p_{T} vs #bar{#Lambda} mass; p_{T} (GeV/#it{c}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kPtVsK0shortMass, o2::framework::kTH2F, "hPtVsK0shortMass", "p_{T} vs K^{0}_{S} mass; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsLambdaMass, o2::framework::kTH2F, "hK0shortMassVsLambdaMass", " K^{0}_{S} mass vs #Lambda mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsAntiLambdaMass, o2::framework::kTH2F, "hK0shortMassVsAntiLambdaMass", "K^{0}_{S} mass vs #bar{#Lambda} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kLambdaMassVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiMass", "p_{T} vs K^{0}_{S} mass; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}"}}};

template <typename T>
auto makeV0AnalysisHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0QaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kCosPa, {confBinningQa.cosPa}},
    {kDecayDauDca, {confBinningQa.dauDcaAtDecay}},
    {kDecayVtxX, {confBinningQa.decayVertex}},
    {kDecayVtxY, {confBinningQa.decayVertex}},
    {kDecayVtxZ, {confBinningQa.decayVertex}},
    {kDecayVtx, {confBinningQa.decayVertex}},
    {kTransRadius, {confBinningQa.transRadius}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsCosPa, {confBinningAnalysis.pt, confBinningQa.cosPa}},
    {kMassLambda, {confBinningQa.massLambda}},
    {kMassAntiLambda, {confBinningQa.massAntiLambda}},
    {kMassK0short, {confBinningQa.massK0short}},
    {kPtVsLambdaMass, {confBinningAnalysis.pt, confBinningQa.massLambda}},
    {kPtVsAntiLambdaMass, {confBinningAnalysis.pt, confBinningQa.massAntiLambda}},
    {kPtVsK0shortMass, {confBinningAnalysis.pt, confBinningQa.massK0short}},
    {kLambdaMassVsAntiLambdaMass, {confBinningQa.massLambda, confBinningQa.massAntiLambda}},
    {kK0shortMassVsLambdaMass, {confBinningQa.massK0short, confBinningQa.massLambda}},
    {kK0shortMassVsAntiLambdaMass, {confBinningQa.massK0short, confBinningQa.massAntiLambda}}};
};

constexpr char PrefixLambdaQa[] = "LambdaQA/";
constexpr char PrefixLambda1[] = "Lambda1/";
constexpr char PrefixK0shortQa[] = "K0shortQa/";
constexpr char PrefixK0short[] = "K0short1/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class V0HistManager
{
 public:
  /// Destructor
  virtual ~V0HistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<V0Hist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {Specs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {Specs[kSign]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(prefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kCosPa, HistTable), GetHistDesc(kCosPa, HistTable), GetHistType(kCosPa, HistTable), {Specs[kCosPa]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayDauDca, HistTable), GetHistDesc(kDecayDauDca, HistTable), GetHistType(kDecayDauDca, HistTable), {Specs[kDecayDauDca]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxX, HistTable), GetHistDesc(kDecayVtxX, HistTable), GetHistType(kDecayVtxX, HistTable), {Specs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxY, HistTable), GetHistDesc(kDecayVtxY, HistTable), GetHistType(kDecayVtxY, HistTable), {Specs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxZ, HistTable), GetHistDesc(kDecayVtxZ, HistTable), GetHistType(kDecayVtxZ, HistTable), {Specs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtx, HistTable), GetHistDesc(kDecayVtx, HistTable), GetHistType(kDecayVtx, HistTable), {Specs[kDecayVtx]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTransRadius, HistTable), GetHistDesc(kTransRadius, HistTable), GetHistType(kTransRadius, HistTable), {Specs[kTransRadius]});

      // qa 2d
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsCosPa, HistTable), GetHistDesc(kPtVsCosPa, HistTable), GetHistType(kPtVsCosPa, HistTable), {Specs[kPtVsCosPa]});

      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassLambda, HistTable), GetHistDesc(kMassLambda, HistTable), GetHistType(kMassLambda, HistTable), {Specs[kMassLambda]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassAntiLambda, HistTable), GetHistDesc(kMassAntiLambda, HistTable), GetHistType(kMassAntiLambda, HistTable), {Specs[kMassAntiLambda]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassK0short, HistTable), GetHistDesc(kMassK0short, HistTable), GetHistType(kMassK0short, HistTable), {Specs[kMassK0short]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsLambdaMass, HistTable), GetHistDesc(kPtVsLambdaMass, HistTable), GetHistType(kPtVsLambdaMass, HistTable), {Specs[kPtVsLambdaMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsAntiLambdaMass, HistTable), GetHistDesc(kPtVsAntiLambdaMass, HistTable), GetHistType(kPtVsAntiLambdaMass, HistTable), {Specs[kPtVsAntiLambdaMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsK0shortMass, HistTable), GetHistDesc(kPtVsK0shortMass, HistTable), GetHistType(kPtVsK0shortMass, HistTable), {Specs[kPtVsK0shortMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kLambdaMassVsAntiLambdaMass, HistTable), GetHistDesc(kLambdaMassVsAntiLambdaMass, HistTable), GetHistType(kLambdaMassVsAntiLambdaMass, HistTable), {Specs[kLambdaMassVsAntiLambdaMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kK0shortMassVsLambdaMass, HistTable), GetHistDesc(kK0shortMassVsLambdaMass, HistTable), GetHistType(kK0shortMassVsLambdaMass, HistTable), {Specs[kK0shortMassVsLambdaMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kK0shortMassVsAntiLambdaMass, HistTable), GetHistDesc(kK0shortMassVsAntiLambdaMass, HistTable), GetHistType(kK0shortMassVsAntiLambdaMass, HistTable), {Specs[kK0shortMassVsAntiLambdaMass]});
    }
  }

  template <modes::Mode mode, modes::V0 v0, typename T>
  void fill(T const& v0candidate)
  {

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), v0candidate.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), v0candidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), v0candidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), v0candidate.mass());

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), v0candidate.sign());
      }
      if constexpr (isFlagSet(v0, modes::V0::kK0short)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), 0);
      }
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kCosPa, HistTable)), v0candidate.cosPa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), v0candidate.dauDca());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), v0candidate.decayVtxX());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), v0candidate.decayVtxY());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), v0candidate.decayVtxZ());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), v0candidate.decayVtx());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), v0candidate.transRadius());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), v0candidate.pt(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), v0candidate.pt(), v0candidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), v0candidate.phi(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsCosPa, HistTable)), v0candidate.pt(), v0candidate.cosPa());

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        float massLambda, massAntiLambda;
        if (v0candidate.sign() > 0) {
          massLambda = v0candidate.mass();
          massAntiLambda = v0candidate.massAnti();
        } else {
          massLambda = v0candidate.massAnti();
          massAntiLambda = v0candidate.mass();
        }
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassLambda, HistTable)), massLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassAntiLambda, HistTable)), massAntiLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassK0short, HistTable)), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), massLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), massAntiLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kLambdaMassVsAntiLambdaMass, HistTable)), massLambda, massAntiLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.massK0short(), massLambda);
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.massK0short(), massAntiLambda);
      }
      if constexpr (isFlagSet(v0, modes::V0::kK0short)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassLambda, HistTable)), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassAntiLambda, HistTable)), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassK0short, HistTable)), v0candidate.mass());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.mass());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kLambdaMassVsAntiLambdaMass, HistTable)), v0candidate.massLambda(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massAntiLambda());
      }
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace v0histmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_
