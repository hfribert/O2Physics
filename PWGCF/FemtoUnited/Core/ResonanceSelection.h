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

/// \file ResonanceSelection.h
/// \brief Vzero selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_RESONANCESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_RESONANCESELECTION_H_

#include "RecoDecay.h"

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

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

// enum of supported resonances
enum class TwoTrackResonace : o2::aod::femtodatatypes::TwoTrackResonaceType {
  kRho,
  kPhi,
  kKstar,
  kAntiKstar,
  kTwoTrackResonaceMax
};

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
using ConfRhoFilters = ConfTwoTrackResonanceFilters<PrefixRhoFilters>;
using ConfPhiFilters = ConfTwoTrackResonanceFilters<PrefixPhiFilters>;
using ConfKstarFilters = ConfTwoTrackResonanceFilters<PrefixKstarFilters>;

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

  o2::framework::Configurable<float> minMomentumForTof{"minMomentumForTof", 0.5f, "Minimum momentum to required TOF PID (all daughers)"};

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
using ConfRhoBits = ConfTwoTrackResonanceBits<PrefixRhoBits>;
using ConfPhiBits = ConfTwoTrackResonanceBits<PrefixPhiBits>;
using ConfKstarBits = ConfTwoTrackResonanceBits<PrefixKstarBits>;

template <const char* Prefix>
struct ConfResonanceSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for Resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 6.f, "Maximum invariant mass for Resonance"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> mask{"mask", 0, "Bitmask for resonance selection"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonaceMaskType> type{"type", 1, "Resonance type (Rho=1, Phi=2, Kstar=4)"};
};

constexpr const char PrefixRhoSelection1[] = "RhoSelection1";
constexpr const char PrefixPhiSelection1[] = "PhiSelection1";
constexpr const char PrefixKstarSelection1[] = "KstarSelection1";
using ConfRhoSelection1 = ConfResonanceSelection<PrefixRhoSelection1>;
using ConfPhiSelection1 = ConfResonanceSelection<PrefixPhiSelection1>;
using ConfKstarSelection1 = ConfResonanceSelection<PrefixKstarSelection1>;

/// The different selections this task is capable of doing
enum TwoTrackResonanceSels {
  // selection for positive daughter
  kPosDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kPosDauDcaxyMax,  ///< Min. DCA of the positive daughers at primary vertex
  kPosDauDcazMax,   ///< Min. DCA of the positive daughers at primary vertex
  kPosDaughItsPion,
  kPosDaughTpcPion,
  kPosDaughTofPion,
  kPosDaughTpctofPion,
  kPosDaughItsKaon,
  kPosDaughTpcKaon,
  kPosDaughTofKaon,
  kPosDaughTpctofKaon,

  // selection for negative daughter
  kNegDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kNegDauDcaxyMax,  ///< Min. DCA of the positive daughers at primary vertex
  kNegDauDcazMax,   ///< Min. DCA of the positive daughers at primary vertex
  kNegDaughItsPion,
  kNegDaughTpcPion,
  kNegDaughTofPion,
  kNegDaughTpctofPion,
  kNegDaughItsKaon,
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

  void setResonanceType(TwoTrackResonace type)
  {
    mType = type;
    if (mType == TwoTrackResonace::kRho) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mPosDaughterMass = o2::constants::physics::MassPiMinus;
    } else if (mType == TwoTrackResonace::kPhi) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mPosDaughterMass = o2::constants::physics::MassKMinus;
    } else if (mType == TwoTrackResonace::kKstar) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mPosDaughterMass = o2::constants::physics::MassPiMinus;
    } else if (mType == TwoTrackResonace::kAntiKstar) {
      mPosDaughterMass = o2::constants::physics::MassKMinus;
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
    } else {
      LOG(warn) << "Resonance type is not supported";
    }
  }

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

  void setMassLimits(float massMin, float massMax)
  {
    mMassMin = massMin;
    mMassMax = massMax;
  }

  void setPtLimits(float ptMin, float ptMax)
  {
    mPtMin = ptMin;
    mPtMax = ptMax;
  }

  void setEtaLimits(float etaMin, float etaMax)
  {
    mEtaMin = etaMin;
    mEtaMax = etaMax;
  }

  void setPhiLimits(float phiMin, float phiMax)
  {
    mPhiMin = phiMin;
    mPhiMax = phiMax;
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

  float getPt() const { return mPt; };
  float getEta() const { return mEta; };
  float getPhi() const { return mPhi; };
  float getMass() const { return mMass; };
  TwoTrackResonace getType() const { return mType; };

  void setMinimalMomentumForTof(float minimalMomentumForTof)
  {
    mMinimalMomentumForTof = minimalMomentumForTof;
  }

  template <typename T>
  bool hasTofAboveThreshold(T const& positiveDaughter, T const& negativeDaughter) const
  {
    // If track momentum exceeds threshold, we require valid TOF info
    return !(positiveDaughter.p() > mMinimalMomentumForTof && !positiveDaughter.hasTOF()) &&
           !(negativeDaughter.p() > mMinimalMomentumForTof && !negativeDaughter.hasTOF());
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
    this->evaluateObservable(kPosDaughItsPion, posDaughter.itsNsigmaPi());
    this->evaluateObservable(kPosDaughTpcPion, posDaughter.tpcNsigmaPi());
    this->evaluateObservable(kPosDaughTofPion, posDaughter.tofNsigmaPi());
    this->evaluateObservable(kPosDaughTpctofPion, std::hypot(posDaughter.tpcNsigmaPi(), posDaughter.tofNsigmaPi()));
    this->evaluateObservable(kPosDaughItsKaon, posDaughter.itsNsigmaKa());
    this->evaluateObservable(kPosDaughTpcKaon, posDaughter.tpcNsigmaKa());
    this->evaluateObservable(kPosDaughTofKaon, posDaughter.tofNsigmaKa());
    this->evaluateObservable(kPosDaughTpctofKaon, std::hypot(posDaughter.tpcNsigmaKa(), posDaughter.tofNsigmaKa()));

    // negative daughter selections
    this->updateLimits(kNegDauDcaxyMax, posDaughter.pt());
    this->evaluateObservable(kNegDauDcaxyMax, negDaughter.dcaXY());
    this->updateLimits(kNegDauDcazMax, posDaughter.pt());
    this->evaluateObservable(kNegDauDcazMax, negDaughter.dcaZ());
    this->evaluateObservable(kNegDauTpcClsMin, negDaughter.tpcNClsFound());
    this->evaluateObservable(kNegDaughItsPion, negDaughter.itsNsigmaPi());
    this->evaluateObservable(kNegDaughTpcPion, negDaughter.tpcNsigmaPi());
    this->evaluateObservable(kNegDaughTofPion, negDaughter.tofNsigmaPi());
    this->evaluateObservable(kNegDaughTpctofPion, std::hypot(negDaughter.tpcNsigmaPi(), negDaughter.tofNsigmaPi()));

    this->evaluateObservable(kNegDaughItsKaon, negDaughter.itsNsigmaKa());
    this->evaluateObservable(kNegDaughTpcKaon, negDaughter.tpcNsigmaKa());
    this->evaluateObservable(kNegDaughTofKaon, negDaughter.tofNsigmaKa());
    this->evaluateObservable(kNegDaughTpctofKaon, std::hypot(negDaughter.tpcNsigmaKa(), negDaughter.tofNsigmaKa()));

    this->assembleBitmask();
  };

  bool checkDaughterPids()
  {
    if (mType == TwoTrackResonace::kRho) {
      return (this->passesOptionalCut(kPosDaughTpcPion) || this->passesOptionalCut(kPosDaughTofPion) || this->passesOptionalCut(kPosDaughTofPion)) &&
             (this->passesOptionalCut(kNegDaughTpcPion) || this->passesOptionalCut(kNegDaughTofPion) || this->passesOptionalCut(kNegDaughTofPion));
    } else if (mType == TwoTrackResonace::kPhi) {
      return (this->passesOptionalCut(kPosDaughTpcKaon) || this->passesOptionalCut(kPosDaughTofKaon) || this->passesOptionalCut(kPosDaughTofKaon)) &&
             (this->passesOptionalCut(kNegDaughTpcKaon) || this->passesOptionalCut(kNegDaughTofKaon) || this->passesOptionalCut(kNegDaughTofKaon));
    } else if (mType == TwoTrackResonace::kKstar) {
      return (this->passesOptionalCut(kPosDaughTpcKaon) || this->passesOptionalCut(kPosDaughTofKaon) || this->passesOptionalCut(kPosDaughTofKaon)) &&
             (this->passesOptionalCut(kNegDaughTpcPion) || this->passesOptionalCut(kNegDaughTofPion) || this->passesOptionalCut(kNegDaughTofPion));
    } else if (mType == TwoTrackResonace::kAntiKstar) {
      return (this->passesOptionalCut(kPosDaughTpcPion) || this->passesOptionalCut(kPosDaughTofPion) || this->passesOptionalCut(kPosDaughTofPion)) &&
             (this->passesOptionalCut(kNegDaughTpcKaon) || this->passesOptionalCut(kNegDaughTofKaon) || this->passesOptionalCut(kNegDaughTofKaon));
    } else {
      LOG(warn) << "Resonance type is not supported";
      return false;
    }
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
  float mMinimalMomentumForTof = 99.f;

  TwoTrackResonace mType;
  float mPosDaughterMass = 0.f;
  float mNegDaughterMass = 0.f;
};

} // namespace twotrackresonanceselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_RESONANCESELECTION_H_
