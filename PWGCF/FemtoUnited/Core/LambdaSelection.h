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

/// \file LambdaSelection.h
/// \brief Lambda selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_LAMBDASELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_LAMBDASELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include "Framework/Configurable.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace lambdaselection
{

struct ConfLambdaFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LambdaFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
};

struct ConfLambdaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LambdaBits");
  // lambda selections
  o2::framework::Configurable<std::vector<float>> dcaDaughMax{"dcaDaughMax", {1.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.99f}, "Minimum cosine of pointing angle"};
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};
  o2::framework::Configurable<std::vector<float>> decayVtxMax{"decayVtxMax", {100.f}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"};
  o2::framework::Configurable<float> kaonMassRejectionLow{"kaonMassRejectionLow", 0.48f, "Lower limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
  o2::framework::Configurable<float> kaonMassRejectionHigh{"kaonMassRejectionHigh", 0.515f, "Upper limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
  // positive daughter selections
  o2::framework::Configurable<std::vector<float>> posDaughDcaMin{"posDaughDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcClustersMin{"posDaughTpcClustersMin", {70.f}, "Minimum number of TPC clusters for daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcPion{"posDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcProton{"posDaughTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC for positive daughter tracks"};
  // negative daughter selections
  o2::framework::Configurable<std::vector<float>> negDaughDcaMin{"negDaughDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcClustersMin{"negDaughTpcClustersMin", {70.f}, "Minimum number of TPC clusters for daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcPion{"negDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcProton{"negDaughTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC negative for daughter tracks"};
};

template <const char* Prefix>
struct ConfLambdaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for Lambda"};
  o2::framework::Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for Lambda"};
  o2::framework::Configurable<float> antiMassMin{"antiMassMin", 0.f, "Minimum invariant mass for AntiLambda"};
  o2::framework::Configurable<float> antiMassMax{"antiMassMax", 999.f, "Maximum invariant mass for AntiLambda"};
  o2::framework::Configurable<o2::aod::femtodatatypes::LambdaMaskType> mask{"mask", 6, "Bitmask for lambda selection"};
};

// Define unique prefixes as constexpr string literals
constexpr const char PrefixLambdaSelection1[] = "LambdaSelection1";
constexpr const char PrefixLambdaSelection2[] = "LambdaSelection2";
constexpr const char PrefixLambdaSelection3[] = "LambdaSelection3";

// Instantiate different instances with unique prefixes
using ConfLambdaSelection1 = ConfLambdaSelection<PrefixLambdaSelection1>;
using ConfLambdaSelection2 = ConfLambdaSelection<PrefixLambdaSelection2>;
using ConfLambdaSelection3 = ConfLambdaSelection<PrefixLambdaSelection3>;

/// The different selections this task is capable of doing
enum LambdaSels {
  // selections for lambdas
  kDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kCpaMin,      ///< Min. CPA (cosine pointing angle)
  kTransRadMin, ///< Min. transverse radius
  kTransRadMax, ///< max. transverse radius
  kDecayVtxMax, ///< Max. distance of decay vertex in x,y,z

  // selection for positive daugther
  kPosDauDcaMin,      ///< Min. DCA of the positive daughers at primary vertex
  kPosDauTpcClsMin,   ///< Min. number of TPC clusters of positive daughter
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter

  // selection for negative daugther
  kNegDauDcaMin,      ///< Min. DCA of the positive daughers at primary vertex
  kNegDauTpcClsMin,   ///< Min. number of TPC clusters of positive daughter
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  kLambdaSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class LambdaSelection : public BaseSelection<float, o2::aod::femtodatatypes::LambdaMaskType, kLambdaSelsMax>
{
 public:
  LambdaSelection() {}
  virtual ~LambdaSelection() = default;

  void setK0ShortMassLimits(float lower, float upper)
  {
    if (lower < 0 || upper < 0) {
      mMassKaonLowerLimit = -1.f;
      mMassKaonUpperLimit = -1.f;
      return;
    }
    mMassKaonLowerLimit = lower;
    mMassKaonUpperLimit = upper;
  }

  template <typename T>
  bool checkK0ShortMass(T const& v0)
  {
    if (mMassKaonLowerLimit < 0 || mMassKaonUpperLimit < 0) {
      return true;
    }
    if (v0.mK0Short() > mMassKaonLowerLimit && v0.mK0Short() < mMassKaonUpperLimit) {
      return false;
    } else {
      return true;
    }
  }

  template <typename V0, typename Tracks>
  void applySelections(V0 const& v0, Tracks const& /*tracks*/)
  {
    this->reset();
    // vzero selections
    this->evaluateObservable(kDcaDaughMax, v0.dcaV0daughters());
    this->evaluateObservable(kCpaMin, v0.v0cosPA());
    this->evaluateObservable(kTransRadMin, v0.v0radius());
    this->evaluateObservable(kTransRadMax, v0.v0radius());
    // for decay vertex, the x,y and z coordinate have to be below a certain threshold
    // compare the largest of the 3 to the limit set by the bit
    std::array<float, 3> decayCoordinates = {v0.x(), v0.y(), v0.z()};
    this->evaluateObservable(kDecayVtxMax, *std::max_element(decayCoordinates.begin(), decayCoordinates.end()));

    // positive daughter selections
    auto posDaughter = v0.template posTrack_as<Tracks>();
    this->evaluateObservable(kPosDauDcaMin, std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()));
    this->evaluateObservable(kPosDauTpcClsMin, posDaughter.tpcNClsFound());
    this->evaluateObservable(kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDaughTpcProton, posDaughter.tpcNSigmaPr());

    // negative daughter selections
    auto negDaughter = v0.template negTrack_as<Tracks>();
    this->evaluateObservable(kNegDauDcaMin, std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ()));
    this->evaluateObservable(kNegDauTpcClsMin, negDaughter.tpcNClsFound());
    this->evaluateObservable(kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->evaluateObservable(kNegDaughTpcProton, negDaughter.tpcNSigmaPr());

    this->assembleBitmask();
  };

  bool checkDaughterPidsLambda()
  {
    return this->passesOptionalCut(kNegDaughTpcPion) && this->passesOptionalCut(kPosDaughTpcProton);
  }

  bool checkDaughterPidsAntiLambda()
  {
    return this->passesOptionalCut(kPosDaughTpcPion) && this->passesOptionalCut(kNegDaughTpcPion);
  }

 protected:
  float mMassKaonLowerLimit = -1.f;
  float mMassKaonUpperLimit = -1.f;
};
} // namespace lambdaselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_LAMBDASELECTION_H_
