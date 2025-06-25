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

#ifndef PWGCF_FEMTOUNITED_CORE_LAMBDAHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_LAMBDAHISTMANAGER_H_

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
namespace lambdahistmanager
{
// enum for track histograms
enum LambdaHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kAntiMass,
  // qa variables
  kCosinePointingAngle,
  kDecayDauDca,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  kKaonMass,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kCosinePointingAngleVsPt,
  // qa for daughters
  kLambdaHistLast
};

template <const char* Prefix>
struct ConfLambdaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{200, 1.0, 1.2}}, "Mass"};
};

constexpr const char PrefixLambdaBinning1[] = "LambdaBinning1";
using ConfLambdaBinning1 = ConfLambdaBinning<PrefixLambdaBinning1>;

template <const char* Prefix>
struct ConfLambdaQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis cosinePointingAngle{"cosinePointingAngle", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis kaonMass{"kaonMass", {{100, 0.45, 0.55}}, "Mass for kaon hypothesis"};
};

constexpr const char PrefixLambdaQaBinning1[] = "LambdaQaBinning1";
using ConfLambdaQaBinning1 = ConfLambdaQaBinning<PrefixLambdaQaBinning1>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<LambdaHist>, kLambdaHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hLambdaMass", "Invariant mass #{Lambda}; m_{p{#pi}^{-}} (GeV/#it{c}^{2}; Entries"},
   {kAntiMass, o2::framework::kTH1F, "hAntiLambdaMass", "hAntiLambdaMass #bar{#{Lambda}};Invariant mass; m_{#bar{p}{#pi}^{+}} (GeV/#it{c}^{2}; Entries"},
   {kCosinePointingAngle, o2::framework::kTH1F, "hCosinePointingAngle", "Cosine of pointing angle; coa(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Tranverse radius ; r_{xy} (cm); Entries"},
   {kKaonMass, o2::framework::kTH1F, "hKaonMass", "Kaon mass hypothesis ; m_{K} (GeV/#it{c}^{2}); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kCosinePointingAngleVsPt, o2::framework::kTH2F, "hCosinePointingAngleVsPt", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"}}};

template <typename T>
std::map<LambdaHist, std::vector<framework::AxisSpec>> makeLambdaHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<LambdaHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kAntiMass, {confBinningAnalysis.mass}}};
};

template <typename T1, typename T2>
std::map<LambdaHist, std::vector<framework::AxisSpec>> makeLambdaQaHistSpecMap(const T1& confBinningAnalysis, const T2 confiBinningQa)
{
  return std::map<LambdaHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kAntiMass, {confBinningAnalysis.mass}},
    {lambdahistmanager::kCosinePointingAngle, {confiBinningQa.cosinePointingAngle}},
    {lambdahistmanager::kDecayDauDca, {confiBinningQa.dauDcaAtDecay}},
    {lambdahistmanager::kDecayVtxX, {confiBinningQa.decayVertex}},
    {lambdahistmanager::kDecayVtxY, {confiBinningQa.decayVertex}},
    {lambdahistmanager::kDecayVtxZ, {confiBinningQa.decayVertex}},
    {lambdahistmanager::kDecayVtx, {confiBinningQa.decayVertex}},
    {lambdahistmanager::kTransRadius, {confiBinningQa.transRadius}},
    {lambdahistmanager::kKaonMass, {confiBinningQa.kaonMass}},
    {lambdahistmanager::kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {lambdahistmanager::kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {lambdahistmanager::kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {lambdahistmanager::kCosinePointingAngleVsPt, {confiBinningQa.cosinePointingAngle, confBinningAnalysis.pt}}};
};

constexpr char PrefixLambdaQa[] = "Lambdas/LambdaQA/";
constexpr char PrefixLambda1[] = "Lambdas/Lambda1/";
constexpr char PrefixLambda2[] = "Lambdas/Lambda2/";
constexpr char PrefixLambda3[] = "Lambdas/Lambda3/";

constexpr char PrefixLambdaCascade[] = "LambdaPosDaughQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class LambdaHistManager
{
 public:
  /// Destructor
  virtual ~LambdaHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<LambdaHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {Specs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kAntiMass, HistTable), GetHistDesc(kAntiMass, HistTable), GetHistType(kAntiMass, HistTable), {Specs[kAntiMass]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(prefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayDauDca, HistTable), GetHistDesc(kDecayDauDca, HistTable), GetHistType(kDecayDauDca, HistTable), {Specs[kDecayDauDca]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxX, HistTable), GetHistDesc(kDecayVtxX, HistTable), GetHistType(kDecayVtxX, HistTable), {Specs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxY, HistTable), GetHistDesc(kDecayVtxY, HistTable), GetHistType(kDecayVtxY, HistTable), {Specs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxZ, HistTable), GetHistDesc(kDecayVtxZ, HistTable), GetHistType(kDecayVtxZ, HistTable), {Specs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtx, HistTable), GetHistDesc(kDecayVtx, HistTable), GetHistType(kDecayVtx, HistTable), {Specs[kDecayVtx]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTransRadius, HistTable), GetHistDesc(kTransRadius, HistTable), GetHistType(kTransRadius, HistTable), {Specs[kTransRadius]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kKaonMass, HistTable), GetHistDesc(kKaonMass, HistTable), GetHistType(kKaonMass, HistTable), {Specs[kKaonMass]});

      // qa 2d
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});
    }
  }

  template <modes::Mode mode, typename T>
  void fill(T const& v0)
  {

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), v0.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), v0.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), v0.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), v0.lambdaMass());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kAntiMass, HistTable)), v0.antiLambdaMass());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), v0.dauDCA());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), v0.decayVtxX());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), v0.decayVtxY());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), v0.decayVtxZ());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), std::hypot(v0.decayVtxX(), v0.decayVtxY(), v0.decayVtxZ()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), v0.transRadius());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kKaonMass, HistTable)), v0.kaonMass());

      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), v0.pt(), v0.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), v0.pt(), v0.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), v0.phi(), v0.eta());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace lambdahistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_LAMBDAHISTMANAGER_H_
