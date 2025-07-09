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
/// \brief histogram manager for vzero histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_

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
namespace v0histmanager
{
// enum for track histograms
enum V0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kLambdaMass,
  kAntiLambdaMass,
  kK0shortMass,
  // qa variables
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

#define V0_BINNING                                                       \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};         \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"}; \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};

template <const char* Prefix>
struct ConfLambdaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis lambdaMass{"lambdaMass", {{200, 1.0, 1.2}}, "Mass"};
  o2::framework::ConfigurableAxis antiLambdaMass{"antiLambdaMass", {{200, 1.0, 1.2}}, "Mass"};
  V0_BINNING
};
constexpr const char PrefixLambdaBinning1[] = "LambdaBinning1";
using ConfLambdaBinning1 = ConfLambdaBinning<PrefixLambdaBinning1>;

template <const char* Prefix>
struct ConfK0shortBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis k0shortMass{"k0shortMass", {{200, 0.47, 0.51}}, "Mass"};
  V0_BINNING
};
constexpr const char PrefixK0shortBinning1[] = "K0shortBinning1";
using ConfK0shortBinning1 = ConfK0shortBinning<PrefixK0shortBinning1>;

#undef V0_BINNING

#define V0_QA_BINNING                                                                                              \
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};                      \
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"}; \
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};                     \
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};

template <const char* Prefix>
struct ConfLambdaQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_QA_BINNING
  o2::framework::ConfigurableAxis k0shortMass{"k0shortMass", {{100, 0.45, 0.55}}, "Mass for kaon hypothesis"};
};
constexpr const char PrefixLambdaQaBinning1[] = "LambdaQaBinning1";
using ConfLambdaQaBinning1 = ConfLambdaQaBinning<PrefixLambdaQaBinning1>;

template <const char* Prefix>
struct ConfK0shortQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_QA_BINNING
};
constexpr const char PrefixK0shortQaBinning1[] = "K0shortQaBinning1";
using ConfK0shortQaBinning1 = ConfK0shortQaBinning<PrefixK0shortQaBinning1>;

#undef V0_QA_BINNING

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<V0Hist>, kV0HistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kLambdaMass, o2::framework::kTH1F, "hLambdaMass", "#{Lambda} mass; m_{p{#pi}^{-}} (GeV/#it{c}^{2}; Entries"},
   {kAntiLambdaMass, o2::framework::kTH1F, "hAntiLambdaMass", "#bar{#{Lambda}} mass;m_{#bar{p}{#pi}^{+}} (GeV/#it{c}^{2}; Entries"},
   {kK0shortMass, o2::framework::kTH1F, "hK0shortMass", "K^{0}_{S} mass hypothesis ; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; coa(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Tranverse radius ; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsLambdaMass, o2::framework::kTH2F, "hPtVsLambdaMass", "p_{T} vs mass #{Lambda}; p_{T} (GeV/#it{c}); m_{p{#pi}^{-}} (GeV/#it{c}^{2}"},
   {kPtVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiLambdaMass", "p_{T} vs mass #bar{#{Lambda}}; p_{T} (GeV/#it{c}); m_{#bar{p}{#pi}^{+}} (GeV/#it{c}^{2}"},
   {kPtVsK0shortMass, o2::framework::kTH2F, "hPtVsK0shortMass", "p_{T} vs mass K^{0}_{S}; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}"},
   {kK0shortMassVsLambdaMass, o2::framework::kTH2F, "hK0shortMassVsLambdaMass", " K^{0}_{S} mass vs #{Lambda} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}; m_{p{#pi}^{-}} (GeV/#it{c}^{2}"},
   {kK0shortMassVsAntiLambdaMass, o2::framework::kTH2F, "hK0shortMassVsAntiLambdaMass", "K^{0}_{S} mass vs #bar{#{Lambda}} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}; m_{#bar{p}{#pi}^{+}} (GeV/#it{c}^{2}"},
   {kLambdaMassVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiMass", "p_{T} vs K^{0}_{S} mass; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}"}}};

template <typename T>
auto makeV0AnalysisHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}}};
};

template <typename T>
auto makeLambdaHistSpecMap(const T& confBinningAnalysis)
{
  auto map = makeV0AnalysisHistSpecMap(confBinningAnalysis);
  map.emplace(kLambdaMass, {confBinningAnalysis.lambdaMass});
  map.emplace(kAntiLambdaMass, {confBinningAnalysis.antiLambdaMass});
  return map;
};

template <typename T>
auto makeKaonHistSpecMap(const T& confBinningAnalysis)
{
  auto map = makeV0AnalysisHistSpecMap(confBinningAnalysis);
  map.emplace(kK0shortMass, {confBinningAnalysis.k0shortMass});
  return map;
};

template <typename T1, typename T2>
auto makeV0QaHistSpecMap(const T1& confBinningAnalysis, const T2 confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
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
    {kPtVsCosPa, {confBinningAnalysis.pt, confBinningQa.cosPa}}};
};

template <typename T1, typename T2>
auto makeLambdaQaHistSpecMap(const T1& confBinningAnalysis, const T2 confBinningQa)
{
  // get default map
  auto map = makeV0QaHistSpecMap(confBinningAnalysis, confBinningQa);
  // add lambda specific histograms
  map.emplace(kK0shortMass, confBinningQa.k0shortMass);
  map.emplace(kPtVsLambdaMass, {confBinningAnalysis.pt, confBinningAnalysis.lambdaMass});
  map.emplace(kPtVsAntiLambdaMass, {confBinningAnalysis.pt, confBinningAnalysis.antiLambdaMass});
  map.emplace(kLambdaMassVsAntiLambdaMass, {confBinningAnalysis.lambdaMass, confBinningAnalysis.antiLambdaMass});
  map.emplace(kK0shortMassVsLambdaMass, {confBinningQa.k0shortMass, confBinningAnalysis.lambdaMass});
  map.emplace(kK0shortMassVsAntiLambdaMass, {confBinningQa.k0shortMass, confBinningAnalysis.antiLambdaMass});
  return map;
};

template <typename T1, typename T2>
auto makeKaonQaHistSpecMap(const T1& confBinningAnalysis, const T2 confBinningQa)
{
  // get default map
  auto map = makeV0QaHistSpecMap(confBinningAnalysis, confBinningQa);
  // add k0short specific histograms
  map.emplace(kPtVsK0shortMass, {confBinningAnalysis.pt, confBinningAnalysis.k0shortMass});
  return map;
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
class v0histmanager
{
 public:
  /// Destructor
  virtual ~v0histmanager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode, modes::V0 v0>
  void init(o2::framework::HistogramRegistry* registry, std::map<V0Hist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        mHistogramRegistry->add(analysisDir + GetHistNamev2(kLambdaMass, HistTable), GetHistDesc(kLambdaMass, HistTable), GetHistType(kLambdaMass, HistTable), {Specs[kLambdaMass]});
        mHistogramRegistry->add(analysisDir + GetHistNamev2(kAntiLambdaMass, HistTable), GetHistDesc(kAntiLambdaMass, HistTable), GetHistType(kAntiLambdaMass, HistTable), {Specs[kAntiLambdaMass]});
      }
      if constexpr (modes::isFlagSet(v0, modes::V0::kK0short)) {
        mHistogramRegistry->add(analysisDir + GetHistNamev2(kK0shortMass, HistTable), GetHistDesc(kK0shortMass, HistTable), GetHistType(kK0shortMass, HistTable), {Specs[kK0shortMass]});
      }
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

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        mHistogramRegistry->add(qaDir + GetHistNamev2(kK0shortMass, HistTable), GetHistDesc(kK0shortMass, HistTable), GetHistType(kK0shortMass, HistTable), {Specs[kK0shortMass]});

        mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsLambdaMass, HistTable), GetHistDesc(kPtVsLambdaMass, HistTable), GetHistType(kPtVsLambdaMass, HistTable), {Specs[kPtVsLambdaMass]});
        mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsAntiLambdaMass, HistTable), GetHistDesc(kPtVsAntiLambdaMass, HistTable), GetHistType(kPtVsAntiLambdaMass, HistTable), {Specs[kPtVsAntiLambdaMass]});
        mHistogramRegistry->add(qaDir + GetHistNamev2(kLambdaMassVsAntiLambdaMass, HistTable), GetHistDesc(kLambdaMassVsAntiLambdaMass, HistTable), GetHistType(kLambdaMassVsAntiLambdaMass, HistTable), {Specs[kLambdaMassVsAntiLambdaMass]});
        mHistogramRegistry->add(qaDir + GetHistNamev2(kK0shortMassVsLambdaMass, HistTable), GetHistDesc(kK0shortMassVsLambdaMass, HistTable), GetHistType(kK0shortMassVsLambdaMass, HistTable), {Specs[kK0shortMassVsLambdaMass]});
        mHistogramRegistry->add(qaDir + GetHistNamev2(kK0shortMassVsAntiLambdaMass, HistTable), GetHistDesc(kK0shortMassVsAntiLambdaMass, HistTable), GetHistType(kK0shortMassVsAntiLambdaMass, HistTable), {Specs[kK0shortMassVsAntiLambdaMass]});
      }

      if constexpr (isFlagSet(v0, modes::V0::kK0short)) {
        mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsK0shortMass, HistTable), GetHistDesc(kPtVsK0shortMass, HistTable), GetHistType(kPtVsK0shortMass, HistTable), {Specs[kPtVsK0shortMass]});
      }
    }
  }

  template <modes::Mode mode, modes::V0 v0, typename T>
  void fill(T const& lambda)
  {

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), lambda.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), lambda.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), lambda.phi());

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kLambdaMass, HistTable)), lambda.lambdaMass());
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kAntiLambdaMass, HistTable)), lambda.antiLambdaMass());
      }
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kCosPa, HistTable)), lambda.cosPa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), lambda.dauDca());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), lambda.decayVtxX());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), lambda.decayVtxY());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), lambda.decayVtxZ());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), std::hypot(lambda.decayVtxX(), lambda.decayVtxY(), lambda.decayVtxZ()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), lambda.transRadius());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), lambda.pt(), lambda.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), lambda.pt(), lambda.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), lambda.phi(), lambda.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsCosPa, HistTable)), lambda.pt(), lambda.cosPa());

      if constexpr (isFlagSet(v0, modes::V0::kLambda)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kK0shortMass, HistTable)), lambda.kaonMass());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsLambdaMass, HistTable)), lambda.pt(), lambda.lambdaMass());
        mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsAntiLambdaMass, HistTable)), lambda.pt(), lambda.antiLambdaMass());
      }
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace v0histmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_V0HISTMANAGER_H_
