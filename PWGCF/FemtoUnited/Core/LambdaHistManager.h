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
  kPosDauPt,
  kPosDauEta,
  kPosDauPhi,
  kNegDauPt,
  kNegDauEta,
  kNegDauPhi,
  // qa variables
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
  // qa for daughters
  kPosDauTpcCluster,
  kPosDauPtVsDcaxy,
  kPosDauPtVsDcaz,
  kPosDauPtVsDca,
  kPosDauProtonTpcNsigma,
  kPosDauPionTpcNsigma,
  kNegDauTpcCluster,
  kNegDauPtVsDcaxy,
  kNegDauPtVsDcaz,
  kNegDauPtVsDca,
  kNegDauProtonTpcNsigma,
  kNegDauPionTpcNsigma,
  kLambdaHistLast
};

constexpr std::string_view AnalysisDir = "LambdaHistograms/Analysis/";
constexpr std::string_view QaDir = "LambdaHistograms/QA/";
constexpr std::string_view PidDir = "LambdaHistograms/Daughters/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<LambdaHist>, kLambdaHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hLambdaMass", "Invariant mass #{Lambda}; m_{p{#pi}^{-}} (GeV/#it{c}^{2}; Entries"},
   {kAntiMass, o2::framework::kTH1F, "hAntiLambdaMass", "hAntiLambdaMass #bar{#{Lambda}};Invariant mass; m_{#bar{p}{#pi}^{+}} (GeV/#it{c}^{2}; Entries"},
   {kPosDauPt, o2::framework::kTH1F, "hPosDauPt", "Transverse Momentum of positive daughter; p_{T} (GeV/#it{c}); Entries"},
   {kPosDauEta, o2::framework::kTH1F, "hPosDauEta", "Pseudorapdity of positive daughter; #eta; Entries"},
   {kPosDauPhi, o2::framework::kTH1F, "hPosDauPhi", "Azimuthal angle of positive daughter; #varphi; Entries"},
   {kNegDauPt, o2::framework::kTH1F, "hNegDauPt", "Transverse Momentum of negative daughter; p_{T} (GeV/#it{c}); Entries"},
   {kNegDauEta, o2::framework::kTH1F, "hNegDauEta", "Pseudorapdity of negative daughter; #eta; Entries"},
   {kNegDauPhi, o2::framework::kTH1F, "hNegDauPhi", "Azimuthal angle of negative daughter; #varphi; Entries"},
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
   {kPosDauTpcCluster, o2::framework::kTH1F, "hPosDauTpcCluster", "TPC cluster found (daughter^{+}) ; TPC cluster found, Entries"},
   {kPosDauPtVsDcaxy, o2::framework::kTH2F, "hPosDauPtVsDcaxy", "p_{T} vs DCA_{XY} (daughter +) ; p_{T} (GeV/#it{c}) ; DCA_{XY} (cm)"},
   {kPosDauPtVsDcaz, o2::framework::kTH2F, "hPosDauPtVsDcaz", "p_{T} vs DCA_{Z} (daughter +) ; p_{T} (GeV/#it{c}) ; DCA_{Z}"},
   {kPosDauPtVsDca, o2::framework::kTH2F, "hPosDauPtVsDca", "p_{T} vs DCA (daughter +) ; p_{T} (GeV/#it{c}) ; DCA"},
   {kPosDauProtonTpcNsigma, o2::framework::kTH2F, "hPosDauProtonTpcNsigma", "TPC Proton PID (daughter +) ; p (GeV/#it{c}) ; n#sigma_{TPC}"},
   {kPosDauPionTpcNsigma, o2::framework::kTH2F, "hPosDauPionTpcNsigma", "TPC Pion PID (daughter +) ; p (GeV/#it{c}) ; n#sigma_{TPC}"},
   {kNegDauTpcCluster, o2::framework::kTH1F, "hNegDauTpcCluster", "TPC cluster found (daughter -) ; TPC cluster found; Entries"},
   {kNegDauPtVsDcaxy, o2::framework::kTH2F, "hNegDauPtVsDcaxy", "p_{T} vs DCA_{XY} (daughter -) ; p_{T} (GeV/#it{c}) ; DCA_{XY} (cm)"},
   {kNegDauPtVsDcaz, o2::framework::kTH2F, "hNegDauPtVsDcaz", "p_{T} vs DCA_{Z} (daughter -) ; p_{T} (GeV/#it{c}) ; DCA_{Z}"},
   {kNegDauPtVsDca, o2::framework::kTH2F, "hNegDauPtVsDca", "p_{T} vs DCA (daughter -) ; p_{T} (GeV/#it{c}) ; DCA"},
   {kNegDauProtonTpcNsigma, o2::framework::kTH2F, "hNegDauProtonTpcNsigma", "TPC Proton PID (daughter -) ; p (GeV/#it{c}) ; n#sigma_{TPC}"},
   {kNegDauPionTpcNsigma, o2::framework::kTH2F, "hNegDauPionTpcNsigma", "TPC Pion PID (daughter -) ; p (GeV/#it{c}) ; n#sigma_{TPC}"}}};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
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

    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {Specs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kAntiMass, HistTable), GetHistDesc(kAntiMass, HistTable), GetHistType(kAntiMass, HistTable), {Specs[kAntiMass]});

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauPt, HistTable), GetHistDesc(kPosDauPt, HistTable), GetHistType(kPosDauPt, HistTable), {Specs[kPosDauPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauEta, HistTable), GetHistDesc(kPosDauEta, HistTable), GetHistType(kPosDauEta, HistTable), {Specs[kPosDauEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauPhi, HistTable), GetHistDesc(kPosDauPhi, HistTable), GetHistType(kPosDauPhi, HistTable), {Specs[kPosDauPhi]});

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauPt, HistTable), GetHistDesc(kNegDauPt, HistTable), GetHistType(kNegDauPt, HistTable), {Specs[kNegDauPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauEta, HistTable), GetHistDesc(kNegDauEta, HistTable), GetHistType(kNegDauEta, HistTable), {Specs[kNegDauEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauPhi, HistTable), GetHistDesc(kNegDauPhi, HistTable), GetHistType(kNegDauPhi, HistTable), {Specs[kNegDauPhi]});
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(QaDir);

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

      // qa daughters
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauTpcCluster, HistTable), GetHistDesc(kPosDauTpcCluster, HistTable), GetHistType(kPosDauTpcCluster, HistTable), {Specs[kPosDauTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauProtonTpcNsigma, HistTable), GetHistDesc(kPosDauProtonTpcNsigma, HistTable), GetHistType(kPosDauProtonTpcNsigma, HistTable), {Specs[kPosDauProtonTpcNsigma]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauPionTpcNsigma, HistTable), GetHistDesc(kPosDauPionTpcNsigma, HistTable), GetHistType(kPosDauPionTpcNsigma, HistTable), {Specs[kPosDauPionTpcNsigma]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauPtVsDcaxy, HistTable), GetHistDesc(kPosDauPtVsDcaxy, HistTable), GetHistType(kPosDauPtVsDcaxy, HistTable), {Specs[kPosDauPtVsDcaxy]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauPtVsDcaz, HistTable), GetHistDesc(kPosDauPtVsDcaz, HistTable), GetHistType(kPosDauPtVsDcaz, HistTable), {Specs[kPosDauPtVsDcaz]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDauPtVsDca, HistTable), GetHistDesc(kPosDauPtVsDca, HistTable), GetHistType(kPosDauPtVsDca, HistTable), {Specs[kPosDauPtVsDca]});

      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauTpcCluster, HistTable), GetHistDesc(kNegDauTpcCluster, HistTable), GetHistType(kNegDauTpcCluster, HistTable), {Specs[kNegDauTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauPtVsDcaxy, HistTable), GetHistDesc(kNegDauPtVsDcaxy, HistTable), GetHistType(kNegDauPtVsDcaxy, HistTable), {Specs[kNegDauPtVsDcaxy]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauPtVsDcaz, HistTable), GetHistDesc(kNegDauPtVsDcaz, HistTable), GetHistType(kNegDauPtVsDcaz, HistTable), {Specs[kNegDauPtVsDcaz]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauPtVsDca, HistTable), GetHistDesc(kNegDauPtVsDca, HistTable), GetHistType(kNegDauPtVsDca, HistTable), {Specs[kNegDauPtVsDca]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauProtonTpcNsigma, HistTable), GetHistDesc(kNegDauProtonTpcNsigma, HistTable), GetHistType(kNegDauProtonTpcNsigma, HistTable), {Specs[kNegDauProtonTpcNsigma]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDauPionTpcNsigma, HistTable), GetHistDesc(kNegDauPionTpcNsigma, HistTable), GetHistType(kNegDauPionTpcNsigma, HistTable), {Specs[kNegDauPionTpcNsigma]});
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& v0, T2 /*tracks*/)
  {

    auto posDaughter = v0.template posDauLambda_as<T2>();
    auto negDaughter = v0.template negDauLambda_as<T2>();

    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), v0.pt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), v0.eta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), v0.phi());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), v0.lambdaMass());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kAntiMass, HistTable)), v0.antiLambdaMass());

      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauPt, HistTable)), posDaughter.pt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauEta, HistTable)), posDaughter.eta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauPhi, HistTable)), posDaughter.phi());

      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauPt, HistTable)), negDaughter.pt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauEta, HistTable)), negDaughter.eta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauPhi, HistTable)), negDaughter.phi());
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), v0.dauDCA());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), v0.decayVtxX());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), v0.decayVtxY());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), v0.decayVtxZ());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), std::hypot(v0.decayVtxX(), v0.decayVtxY(), v0.decayVtxZ()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), v0.transRadius());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kKaonMass, HistTable)), v0.kaonMass());

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), v0.pt(), v0.eta());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), v0.pt(), v0.phi());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), v0.phi(), v0.eta());

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauTpcCluster, HistTable)), posDaughter.tpcNClsFound());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauPtVsDcaxy, HistTable)), posDaughter.pt(), posDaughter.dcaXY());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauPtVsDcaz, HistTable)), posDaughter.pt(), posDaughter.dcaZ());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauPtVsDca, HistTable)), posDaughter.pt(), std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauPionTpcNsigma, HistTable)), posDaughter.p(), posDaughter.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDauProtonTpcNsigma, HistTable)), posDaughter.p(), posDaughter.tpcNSigmaPr());

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauTpcCluster, HistTable)), negDaughter.tpcNClsFound());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauPtVsDcaxy, HistTable)), negDaughter.pt(), negDaughter.dcaXY());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauPtVsDcaz, HistTable)), negDaughter.pt(), negDaughter.dcaZ());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauPtVsDca, HistTable)), negDaughter.pt(), std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauPionTpcNsigma, HistTable)), negDaughter.p(), negDaughter.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDauProtonTpcNsigma, HistTable)), negDaughter.p(), negDaughter.tpcNSigmaPr());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace lambdahistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_LAMBDAHISTMANAGER_H_
