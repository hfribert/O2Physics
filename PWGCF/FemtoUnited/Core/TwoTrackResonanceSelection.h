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

/// \file TwoTrackResonanceSelection.h
/// \brief Vzero selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_

#include "RecoDecay.h"

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/Configurable.h"
#include <CommonConstants/PhysicsConstants.h>

#include "Math/Vector4D.h"

#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace twotrackresonanceselection
{

struct ConfTwoTrackResonanceDaughterFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoTrackResonanceDaughterFilter");
  o2::framework::Configurable<float> posDaughPtMin{"posDaughPtMin", 0.f, "Minimum pT for positive Daughter"};
  o2::framework::Configurable<float> posDaughPtMax{"posDaughPtMax", 6.f, "Maximum pT for positive Daughter"};
  o2::framework::Configurable<float> negDaughPtMin{"negDaughPtMin", 0.f, "Minimum pT for negative Daughter"};
  o2::framework::Configurable<float> negDaughPtMax{"negDaughPtMax", 6.f, "Maximum pT for negative Daughter"};

  o2::framework::Configurable<float> posDaughEtaMin{"posDaughEtaMin", -0.9f, "Minimum eta for positive Daughter"};
  o2::framework::Configurable<float> posDaughEtaMax{"posDaughEtaMax", 0.9f, "Maximum eta for positive Daughter"};
  o2::framework::Configurable<float> negDaughEtaMin{"negDaughEtaMin", -0.9f, "Minimum eta for negative Daughter"};
  o2::framework::Configurable<float> negDaughEtaMax{"negDaughEtaMax", 0.9f, "Maximum eta for negative Daughter"};
  o2::framework::Configurable<float> posDaughPhiMin{"posDaughPhiMin", 0.f, "Minimum phi for positive Daughter"};
  o2::framework::Configurable<float> posDaughPhiMax{"posDaughPhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi for positive Daughter"};
  o2::framework::Configurable<float> negDaughPhiMin{"negDaughPhiMin", 0.f, "Minimum phi for negative Daughter"};
  o2::framework::Configurable<float> negDaughPhiMax{"negDaughPhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi for negative Daughter"};
};

template <const char* Prefix>
struct ConfTwoTrackResonanceFilters : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for Resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 6.f, "Maximum invariant mass for Resonance"};
};
constexpr const char PrefixRhoFilters[] = "RhoFilters1";
constexpr const char PrefixPhiFilters[] = "PhiFilters1";
constexpr const char PrefixKstarFilters[] = "KstarFilters1";
constexpr const char PrefixAntiKstarFilters[] = "AntiKstarFilters1";
using ConfRhoFilters = ConfTwoTrackResonanceFilters<PrefixRhoFilters>;
using ConfPhiFilters = ConfTwoTrackResonanceFilters<PrefixPhiFilters>;
using ConfKstarFilters = ConfTwoTrackResonanceFilters<PrefixKstarFilters>;
using ConfAntiKstarFilters = ConfTwoTrackResonanceFilters<PrefixAntiKstarFilters>;

template <const char* Prefix>
struct ConfTwoTrackResonanceBits : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  // track quality cuts for resonance daughters
  o2::framework::Configurable<std::vector<float>> posDaughTpcClustersMin{"posDaughTpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};
  o2::framework::Configurable<std::vector<std::string>> posDaughDcaxyMax{"posDaughDcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"};
  o2::framework::Configurable<std::vector<std::string>> posDaughDcazMax{"posDaughDcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcClustersMin{"negDaughTpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};
  o2::framework::Configurable<std::vector<std::string>> negDaughDcaxyMax{"negDaughDcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"};
  o2::framework::Configurable<std::vector<std::string>> negDaughDcazMax{"negDaughDcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};

  // momentum threshold for TOF PID
  // this is set and fixed in the producer, not in the pair task
  // like this we just need to store 2 additional bits instead of two floats containg the momentum of the daughters
  // for initalization of the partition we cannnot access the information which is linked in another table
  o2::framework::Configurable<float> posDaughMinMomentumForTof{"posDaughMinMomentumForTof", 0.5f, "Minimum momentum to required TOF PID (positive daughers)"};
  o2::framework::Configurable<float> negDaughMinMomentumForTof{"negDaughMinMomentumForTof", 0.5f, "Minimum momentum to required TOF PID (negative daughers)"};

  // track pid cuts for positive daughters
  o2::framework::Configurable<std::vector<float>> posDaughItsPion{"posDaughItsPion", {}, "Maximum |nsimga_Pion| ITS for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcPion{"posDaughTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTofPion{"posDaughTofPion", {}, "Maximum |nsimga_Pion| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpctofPion{"posDaughTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughItsKaon{"posDaughItsKaon", {}, "Maximum |nsimga_Kaon| ITS for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcKaon{"posDaughTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTofKaon{"posDaughTofKaon", {}, "Maximum |nsimga_Kaon| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpctofKaon{"posDaughTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for positive daughter tracks"};

  // track pid cuts for negative daughters
  o2::framework::Configurable<std::vector<float>> negDaughItsPion{"negDaughItsPion", {}, "Maximum |nsimga_Pion| ITS for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcPion{"negDaughTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTofPion{"negDaughTofPion", {}, "Maximum |nsimga_Pion| TOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpctofPion{"negDaughTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughItsKaon{"negDaughItsKaon", {}, "Maximum |nsimga_Kaon| ITS for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcKaon{"negDaughTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTofKaon{"negDaughTofKaon", {}, "Maximum |nsimga_Kaon| TOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpctofKaon{"negDaughTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for negative daughter tracks"};
};
constexpr const char PrefixRhoBits[] = "RhoBits";
constexpr const char PrefixPhiBits[] = "PhiBits";
constexpr const char PrefixKstarBits[] = "KstarBits";
constexpr const char PrefixAntiKstarBits[] = "AntiKstarBits";
using ConfRhoBits = ConfTwoTrackResonanceBits<PrefixRhoBits>;
using ConfPhiBits = ConfTwoTrackResonanceBits<PrefixPhiBits>;
using ConfKstarBits = ConfTwoTrackResonanceBits<PrefixKstarBits>;
using ConfAntiKstarBits = ConfTwoTrackResonanceBits<PrefixAntiKstarBits>;

template <const char* Prefix>
struct ConfTwoTrackResonaceSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for Resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 6.f, "Maximum invariant mass for Resonance"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> maskPosDauBelowThres{"maskPosDauBelowThres", 0u, "Bitmask for resonance selection"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> maskPosDauAboveThres{"maskPosDauAboveThres", 0u, "Bitmask for resonance selection"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> maskNegDauBelowThres{"maskNegDauBelowThres", 0u, "Bitmask for resonance selection"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> maskNegDauAboveThres{"maskNegDauAboveThres", 0u, "Bitmask for resonance selection"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceType> type{"type", 4u, "Resonance type (Rho=0, Phi=1, Kstar=2, AntiKstar=3)"};
};

constexpr const char PrefixTwoTrackResonanceSelection1[] = "TwoTrackResonanceSelection1";
constexpr const char PrefixRhoSelection1[] = "RhoSelection1";
constexpr const char PrefixPhiSelection1[] = "PhiSelection1";
constexpr const char PrefixKstarSelection1[] = "KstarSelection1";
constexpr const char PrefixAntiKstarSelection1[] = "AntiKstarSelection1";
using ConfTwoTrackResonaceSelection1 = ConfTwoTrackResonaceSelection<PrefixTwoTrackResonanceSelection1>;
using ConfRhoSelection1 = ConfTwoTrackResonaceSelection<PrefixRhoSelection1>;
using ConfPhiSelection1 = ConfTwoTrackResonaceSelection<PrefixPhiSelection1>;
using ConfKstarSelection1 = ConfTwoTrackResonaceSelection<PrefixKstarSelection1>;
using ConfAntiKstarSelection1 = ConfTwoTrackResonaceSelection<PrefixAntiKstarSelection1>;

/// The different selections this task is capable of doing
enum TwoTrackResonanceSels {
  // selection for positive daughter
  kPosDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kPosDauDcaxyMax,  ///< Min. DCA of the positive daughers at primary vertex
  kPosDauDcazMax,   ///< Min. DCA of the positive daughers at primary vertex
  kPosDaughTpcPion,
  kPosDaughTofPion,
  kPosDaughTpctofPion,
  kPosDaughTpcKaon,
  kPosDaughTofKaon,
  kPosDaughTpctofKaon,

  // selection for negative daughter
  kNegDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kNegDauDcaxyMax,  ///< Min. DCA of the positive daughers at primary vertex
  kNegDauDcazMax,   ///< Min. DCA of the positive daughers at primary vertex
  kNegDaughTpcPion,
  kNegDaughTofPion,
  kNegDaughTpctofPion,
  kNegDaughTpcKaon,
  kNegDaughTofKaon,
  kNegDaughTpctofKaon,

  kResonanceSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TwoTrackResonanceSelection : public BaseSelection<float, o2::aod::femtodatatypes::TwoTrackResonaceMaskType, twotrackresonanceselection::kResonanceSelsMax>
{
 public:
  TwoTrackResonanceSelection() {}
  virtual ~TwoTrackResonanceSelection() = default;

  template <o2::analysis::femtounited::modes::TwoTrackResonace reso, typename T1, typename T2, typename T3>
  void configure(T1 const& config, T2 const& filter, T3 const& daughterFilter)
  {
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(reso, o2::analysis::femtounited::modes::TwoTrackResonace::kPhi)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;

      this->addSelection(config.posDaughTpcKaon.value, kPosDaughTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTofKaon.value, kPosDaughTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTpctofKaon.value, kPosDaughTpctofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpcKaon.value, kNegDaughTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTofKaon.value, kNegDaughTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpctofKaon.value, kNegDaughTpctofKaon, limits::kAbsUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(reso, o2::analysis::femtounited::modes::TwoTrackResonace::kRho)) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;

      this->addSelection(config.posDaughTpcPion.value, kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTofPion.value, kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTpctofPion.value, kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpcPion.value, kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTofPion.value, kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpctofPion.value, kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(reso, o2::analysis::femtounited::modes::TwoTrackResonace::kKstar)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;

      this->addSelection(config.posDaughTpcKaon.value, kPosDaughTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTofKaon.value, kPosDaughTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTpctofKaon.value, kPosDaughTpctofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpcPion.value, kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTofPion.value, kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpctofPion.value, kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(reso, o2::analysis::femtounited::modes::TwoTrackResonace::kKstarBar)) {
      mPosDaughterMass = o2::constants::physics::MassKMinus;
      mNegDaughterMass = o2::constants::physics::MassPiPlus;

      this->addSelection(config.posDaughTpcPion.value, kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTofPion.value, kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDaughTpctofPion.value, kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpcKaon.value, kNegDaughTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTofKaon.value, kNegDaughTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDaughTpctofKaon.value, kNegDaughTpctofKaon, limits::kAbsUpperLimit, false, false);
    }
    mMassMin = filter.massMin.value;
    mMassMax = filter.massMax.value;
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    mPosDaughMinimalMomentumForTof = config.posDaughMinMomentumForTof.value;
    mNegDaughMinimalMomentumForTof = config.negDaughMinMomentumForTof.value;

    this->addSelection(config.posDaughTpcClustersMin.value, kPosDauTpcClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.posDaughDcaxyMax.name, daughterFilter.posDaughPtMin.value, daughterFilter.posDaughPtMax.value, config.posDaughDcaxyMax.value, kPosDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.posDaughDcazMax.name, daughterFilter.posDaughPtMin.value, daughterFilter.posDaughPtMax.value, config.posDaughDcazMax.value, kPosDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    // rho negative daughter selections
    this->addSelection(config.negDaughTpcClustersMin.value, kNegDauTpcClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.negDaughDcaxyMax.name, daughterFilter.negDaughPtMin.value, daughterFilter.negDaughPtMax.value, config.negDaughDcaxyMax.value, kNegDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.negDaughDcazMax.name, daughterFilter.negDaughPtMin.value, daughterFilter.negDaughPtMax.value, config.negDaughDcazMax.value, kNegDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
  };

  template <typename Tracks>
  void reconstructResonance(Tracks const& posDaughter, Tracks const& negDaughter)
  {

    ROOT::Math::PtEtaPhiMVector vecPosDaughter{posDaughter.pt(), posDaughter.eta(), posDaughter.phi(), mPosDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecNegDaughter{negDaughter.pt(), negDaughter.eta(), negDaughter.phi(), mNegDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecResonance = vecPosDaughter + vecNegDaughter;

    // cache kinematics
    mMass = vecResonance.M();
    mPt = vecResonance.Pt();
    mEta = vecResonance.Eta();
    mPhi = RecoDecay::constrainAngle(vecResonance.Phi());
  }

  bool checkFilters() const
  {
    if (mMass < mMassMin || mMass > mMassMax ||
        mPt < mPtMin || mPt > mPtMax ||
        mEta < mEtaMin || mEta > mEtaMax ||
        mPhi < mPhiMin || mPhi > mPhiMax) {
      return false;
    }
    return true;
  }

  float getPt() const { return mPt; }
  float getEta() const { return mEta; }
  float getPhi() const { return mPhi; }
  float getMass() const { return mMass; }

  template <typename T>
  bool hasTofAboveThreshold(T const& positiveDaughter, T const& negativeDaughter) const
  {
    // If track momentum exceeds threshold, we require valid TOF info
    return !(positiveDaughter.p() > mPosDaughMinimalMomentumForTof && !positiveDaughter.hasTOF()) &&
           !(negativeDaughter.p() > mNegDaughMinimalMomentumForTof && !negativeDaughter.hasTOF());
  }

  template <typename Tracks>
  void applySelections(Tracks const& posDaughter, Tracks const& negDaughter)
  {
    this->reset();
    // for resoanace topological selectsion are in general not possible, so only selections on the daughters are performed

    // positive daughter selections
    this->updateLimits(kPosDauDcaxyMax, posDaughter.pt());
    this->evaluateObservable(kPosDauDcaxyMax, posDaughter.dcaXY());
    this->updateLimits(kPosDauDcazMax, posDaughter.pt());
    this->evaluateObservable(kPosDauDcazMax, posDaughter.dcaZ());
    this->evaluateObservable(kPosDauTpcClsMin, posDaughter.tpcNClsFound());
    this->evaluateObservable(kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDaughTofPion, posDaughter.tofNSigmaPi());
    this->evaluateObservable(kPosDaughTpctofPion, std::hypot(posDaughter.tpcNSigmaPi(), posDaughter.tofNSigmaPi()));
    this->evaluateObservable(kPosDaughTpcKaon, posDaughter.tpcNSigmaKa());
    this->evaluateObservable(kPosDaughTofKaon, posDaughter.tofNSigmaKa());
    this->evaluateObservable(kPosDaughTpctofKaon, std::hypot(posDaughter.tpcNSigmaKa(), posDaughter.tofNSigmaKa()));

    // negative daughter selections
    this->updateLimits(kNegDauDcaxyMax, negDaughter.pt());
    this->evaluateObservable(kNegDauDcaxyMax, negDaughter.dcaXY());
    this->updateLimits(kNegDauDcazMax, negDaughter.pt());
    this->evaluateObservable(kNegDauDcazMax, negDaughter.dcaZ());
    this->evaluateObservable(kNegDauTpcClsMin, negDaughter.tpcNClsFound());
    this->evaluateObservable(kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->evaluateObservable(kNegDaughTofPion, negDaughter.tofNSigmaPi());
    this->evaluateObservable(kNegDaughTpctofPion, std::hypot(negDaughter.tpcNSigmaPi(), negDaughter.tofNSigmaPi()));

    this->evaluateObservable(kNegDaughTpcKaon, negDaughter.tpcNSigmaKa());
    this->evaluateObservable(kNegDaughTofKaon, negDaughter.tofNSigmaKa());
    this->evaluateObservable(kNegDaughTpctofKaon, std::hypot(negDaughter.tpcNSigmaKa(), negDaughter.tofNSigmaKa()));

    this->assembleBitmask();
  };

  bool checkRhoHypothesis()
  {
    return (this->passesOptionalCut(kPosDaughTpcPion) || this->passesOptionalCut(kPosDaughTofPion) || this->passesOptionalCut(kPosDaughTofPion)) &&
           (this->passesOptionalCut(kNegDaughTpcPion) || this->passesOptionalCut(kNegDaughTofPion) || this->passesOptionalCut(kNegDaughTofPion));
  }

  bool checkPhiHypothesis()
  {
    return (this->passesOptionalCut(kPosDaughTpcKaon) || this->passesOptionalCut(kPosDaughTofKaon) || this->passesOptionalCut(kPosDaughTofKaon)) &&
           (this->passesOptionalCut(kNegDaughTpcKaon) || this->passesOptionalCut(kNegDaughTofKaon) || this->passesOptionalCut(kNegDaughTofKaon));
  }

  bool checkKstarHypothesis()
  {
    return (this->passesOptionalCut(kPosDaughTpcKaon) || this->passesOptionalCut(kPosDaughTofKaon) || this->passesOptionalCut(kPosDaughTofKaon)) &&
           (this->passesOptionalCut(kNegDaughTpcPion) || this->passesOptionalCut(kNegDaughTofPion) || this->passesOptionalCut(kNegDaughTofPion));
  }
  bool checkKstarBarHypothesis()
  {
    return (this->passesOptionalCut(kPosDaughTpcPion) || this->passesOptionalCut(kPosDaughTofPion) || this->passesOptionalCut(kPosDaughTofPion)) &&
           (this->passesOptionalCut(kNegDaughTpcKaon) || this->passesOptionalCut(kNegDaughTofKaon) || this->passesOptionalCut(kNegDaughTofKaon));
  }

 protected:
  // kinematic variables of the resonance
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;

  // kinematic selections of the resonance
  float mMassMin = 0.f;
  float mMassMax = 6.f;
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  // minimum momentum of the daughers to ask for tof information
  float mPosDaughMinimalMomentumForTof = 99.f;
  float mNegDaughMinimalMomentumForTof = 99.f;

  float mPosDaughterMass = 0.f;
  float mNegDaughterMass = 0.f;
};

} // namespace twotrackresonanceselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_
