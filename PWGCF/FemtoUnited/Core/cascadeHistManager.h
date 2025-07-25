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

/// \file cascadeHistManager.h
/// \brief histogram manager for cascade histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_CASCADEHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_CASCADEHISTMANAGER_H_

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
namespace cascadehistmanager
{
// enum for track histograms
enum CascadeHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kMassXi,
  kMassOmega,
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
  kPtVsMassXi,
  kPtVsMassOmega,
  kMassXiVsMassOmega,
  kCascadeHistLast
};

template <const char* Prefix>
struct ConfCascadeBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{500, 1.5f, 2.f}}, "Mass"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
};

constexpr const char PrefixXiBinning[] = "XiBinning";
using ConfXiBinning = ConfCascadeBinning<PrefixXiBinning>;
constexpr const char PrefixOmegaBinning[] = "OmegaBinning";
using ConfOmegaBinning = ConfCascadeBinning<PrefixOmegaBinning>;

template <const char* Prefix>
struct ConfCascadeQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis massXi{"massXi", {{400, 1.2f, 1.6f}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massOmega{"massOmega", {{400, 1.4f, 1.8f}}, "mass for antiparticle hypothesis"};
};

constexpr const char PrefixXiQaBinning[] = "XiQaBinning";
using ConfXiQaBinning = ConfCascadeQaBinning<PrefixXiQaBinning>;

constexpr const char PrefixOmegatQaBinning[] = "OmegaQaBinning";
using ConfOmegaQaBinning = ConfCascadeQaBinning<PrefixOmegatQaBinning>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CascadeHist>, kCascadeHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv}} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kMassXi, o2::framework::kTH1F, "hMassXi", "Mass #Xi; m_{#Lambda#pi} (GeV/#it{c}^{2}); Entries"},
   {kMassOmega, o2::framework::kTH1F, "hMassOmega", "mass #Omega; m_{#LambdaK} (GeV/#it{c}^{2}); Entries"},
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
   {kPtVsMassXi, o2::framework::kTH2F, "hPtVsMassXi", "p_{T} vs mass #Xi; p_{T} (GeV/#it{c}); m_{#Lambda#pi} (GeV/#it{c}^{2})"},
   {kPtVsMassOmega, o2::framework::kTH2F, "hPtVsMassOmega", "p_{T} vs mass #Omega; p_{T} (GeV/#it{c}); m_{#LambdaK} (GeV/#it{c}^{2})"},
   {kMassXiVsMassOmega, o2::framework::kTH2F, "hMassXiVsMassOmega", "mass #Xi vs mass #Omega; m_{#Lambda#pi} (GeV/#it{c}^{2}); m_{#LambdaK} (GeV/#it{c}^{2})"}}};

template <typename T>
auto makeCascadeAnalysisHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CascadeHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<CascadeHist, std::vector<framework::AxisSpec>> makeCascadeQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<CascadeHist, std::vector<framework::AxisSpec>>{
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
    {kMassXi, {confBinningQa.massXi}},
    {kMassOmega, {confBinningQa.massOmega}},
    {kPtVsMassXi, {confBinningAnalysis.pt, confBinningQa.massXi}},
    {kPtVsMassOmega, {confBinningAnalysis.pt, confBinningQa.massOmega}},
    {kMassXiVsMassOmega, {confBinningQa.massXi, confBinningQa.massOmega}}};
};

constexpr char PrefixXiQa[] = "XiQA/";
constexpr char PrefixXi[] = "Xi/";
constexpr char PrefixOmegaQa[] = "OmegaQa/";
constexpr char PrefixOmega[] = "Omega/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class CascadeHistManager
{
 public:
  /// Destructor
  virtual ~CascadeHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> Specs)
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

      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassXi, HistTable), GetHistDesc(kMassXi, HistTable), GetHistType(kMassXi, HistTable), {Specs[kMassXi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassOmega, HistTable), GetHistDesc(kMassOmega, HistTable), GetHistType(kMassOmega, HistTable), {Specs[kMassOmega]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsMassXi, HistTable), GetHistDesc(kPtVsMassXi, HistTable), GetHistType(kPtVsMassXi, HistTable), {Specs[kPtVsMassXi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsMassOmega, HistTable), GetHistDesc(kPtVsMassOmega, HistTable), GetHistType(kPtVsMassOmega, HistTable), {Specs[kPtVsMassOmega]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMassXiVsMassOmega, HistTable), GetHistDesc(kMassXiVsMassOmega, HistTable), GetHistType(kMassXiVsMassOmega, HistTable), {Specs[kMassXiVsMassOmega]});
    }
  }

  template <modes::Mode mode, modes::Cascade cascade, typename T>
  void fill(T const& cascadeCandidate)
  {

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), cascadeCandidate.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), cascadeCandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), cascadeCandidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), cascadeCandidate.mass());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), cascadeCandidate.sign());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQA)) {

      float massXi, massOmega;
      if constexpr (isFlagSet(cascade, modes::Cascade::kXi)) {
        massXi = cascadeCandidate.mass();
        massOmega = cascadeCandidate.massOmega();
      }
      if constexpr (isFlagSet(cascade, modes::Cascade::kOmega)) {
        massXi = cascadeCandidate.massXi();
        massOmega = cascadeCandidate.mass();
      }

      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kCosPa, HistTable)), cascadeCandidate.cascadeCosPa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), cascadeCandidate.cascadeDauDca());
      // mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), cascadeCandidate.decayVtxX());
      // mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), cascadeCandidate.decayVtxY());
      // mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), cascadeCandidate.decayVtxZ());
      // mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), cascadeCandidate.decayVtx());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), cascadeCandidate.cascadeTransRadius());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), cascadeCandidate.pt(), cascadeCandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), cascadeCandidate.pt(), cascadeCandidate.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), cascadeCandidate.phi(), cascadeCandidate.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsCosPa, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassXi, HistTable)), massXi);
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassOmega, HistTable)), massOmega);
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsMassXi, HistTable)), cascadeCandidate.pt(), massXi);
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsMassOmega, HistTable)), cascadeCandidate.pt(), massOmega);
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kMassXiVsMassOmega, HistTable)), massXi, massOmega);
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace cascadehistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CASCADEHISTMANAGER_H_
