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

/// \file V0Selecion.h
/// \brief v0 selections
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/Configurable.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace v0selection
{

// filters applied in the producer task
struct ConfV0Filters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("V0Filters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMinLambda{"massMinLambda", 1.f, "Minimum mass for Lambda hypothesis"};
  o2::framework::Configurable<float> massMaxLambda{"massMaxLambda", 1.2f, "Maximum mass for Lambda hypothesis"};
  o2::framework::Configurable<float> rejectMassMinK0short{"rejectMassMinK0short", 0.48f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxK0short{"rejectMassMaxK0short", 0.515f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> massMinK0short{"massMinK0short", 0.48f, "Minimum mass for K0Short hypothesis"};
  o2::framework::Configurable<float> massMaxK0short{"massMaxK0short", 0.515f, "Maximum mass for K0Short hypothesis"};
  o2::framework::Configurable<bool> fillLambdas{"fillLambdas", true, "Fill Lambda candidates"};
  o2::framework::Configurable<bool> fillK0shorts{"fillK0shorts", true, "Fill K0Short candidates"};
};

// selections bits for all v0s
#define V0_DEFAULT_BITS                                                                                                                                          \
  o2::framework::Configurable<std::vector<float>> dcaDauMax{"dcaDaughMax", {1.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};                    \
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.99f}, "Minimum cosine of pointing angle"};                                                 \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                                          \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                                         \
  o2::framework::Configurable<std::vector<float>> decayVtxMax{"decayVtxMax", {100.f}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"DauAbsEtaMax", {0.8f}, "Minimum DCA of the daughters from primary vertex (cm)"};                 \
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"DauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"DauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};

// derived selection bits for lambda
struct ConfLambdaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LambdaBits");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTpcProton{"posDaughTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcProton{"negDaughTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC negative for daughter tracks"};
};

// derived selection bits for K0Short
struct ConfK0shortBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("K0short");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDaughTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
};

#undef V0_DEFAULT_BITS

// base selection for analysis task for v0s
#define V0_DEFAULT_SELECTIONS                                                                           \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                 \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                               \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                            \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                             \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                              \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::V0MaskType> mask{"mask", 6, "Bitmask for v0 selection"};

// base selection for analysis task for lambdas
template <const char* Prefix>
struct ConfLambdaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS
  o2::framework::Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for Lambda"};
  o2::framework::Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for Lambda"};
  o2::framework::Configurable<float> antiMassMin{"antiMassMin", 0.f, "Minimum invariant mass for AntiLambda"};
  o2::framework::Configurable<float> antiMassMax{"antiMassMax", 999.f, "Maximum invariant mass for AntiLambda"};
};

// Define unique prefixes as constexpr string literals
constexpr const char PrefixLambdaSelection1[] = "LambdaSelection1";
using ConfLambdaSelection1 = ConfLambdaSelection<PrefixLambdaSelection1>;

// base selection for analysis task for k0Short
template <const char* Prefix>
struct ConfK0shortSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS
  o2::framework::Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for K0Short"};
  o2::framework::Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for K0Short"};
};

#undef V0_DEFAULT_SELECTIONS

// Define unique prefixes as constexpr string literals
constexpr const char PrefixK0ShortSelection1[] = "K0ShortSelection1";
using ConfK0shortSelection1 = ConfK0shortSelection<PrefixK0ShortSelection1>;

/// The different selections for v0s
enum V0Seles {
  // selections for lambdas
  kCpaMin,      ///< Min. CPA (cosine pointing angle)
  kDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kDecayVtxMax, ///< Max. distance of decay vertex in x,y,z
  kTransRadMin, ///< Min. transverse radius
  kTransRadMax, ///< max. transverse radius

  // selection for daugther
  kDauAbsEtaMax, ///< Max. absolute pseudo rapidity
  kDauDcaMin,    ///< Min. DCA of the positive daughers at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter

  // pid selection for daughters
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  kV0SelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class V0Selection : public BaseSelection<float, o2::aod::femtodatatypes::V0MaskType, kV0SelsMax>
{
 public:
  V0Selection() {}
  virtual ~V0Selection() = default;

  template <o2::analysis::femtounited::modes::V0 v0, typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(v0, o2::analysis::femtounited::modes::V0::kLambda)) {
      mMassLambdaLowerLimit = filter.massMinLambda.value;
      mMassLambdaUpperLimit = filter.massMaxLambda.value;
      mMassK0shortLowerLimit = filter.rejectMassMinK0short.value;
      mMassK0shortUpperLimit = filter.rejectMassMaxK0short.value;
    }
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(v0, o2::analysis::femtounited::modes::V0::kK0short)) {
      mMassK0shortLowerLimit = filter.massMinK0short.value;
      mMassK0shortUpperLimit = filter.massMaxK0short.value;
    }

    this->addSelection(config.dcaDauMax.value, v0selection::kDcaDaughMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.cpaMin.value, v0selection::kCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMin.value, v0selection::kTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMax.value, v0selection::kTransRadMax, limits::kUpperLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, v0selection::kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaMin.value, v0selection::kDauDcaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, v0selection::kDauTpcClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.posDauTpcPion.value, v0selection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.posDauTpcProton.value, v0selection::kPosDaughTpcProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.negDauTpcPion.value, v0selection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.negDauTpcProton.value, v0selection::kNegDaughTpcProton, limits::kAbsUpperLimit, false, false);
  };

  template <typename V0, typename Tracks>
  void applySelections(V0 const& v0, Tracks const& /*tracks*/)
  {
    this->reset();
    // v0 selections
    this->evaluateObservable(kCpaMin, v0.v0cosPA());
    this->evaluateObservable(kDcaDaughMax, v0.dcaV0daughters());
    // for decay vertex, the x,y and z coordinate have to be below a certain threshold
    // compare the largest of the 3 to the limit set by the bit
    std::array<float, 3> decayCoordinates = {v0.x(), v0.y(), v0.z()};
    this->evaluateObservable(kDecayVtxMax, *std::max_element(decayCoordinates.begin(), decayCoordinates.end()));
    this->evaluateObservable(kTransRadMin, v0.v0radius());
    this->evaluateObservable(kTransRadMax, v0.v0radius());

    // daughter selection
    // for daughter selections, both have to fit the same track quality selection, so we store only one bit for both
    // take largest/smallest from both daughters and evalute the observable with this value
    auto posDaughter = v0.template posTrack_as<Tracks>();
    auto negDaughter = v0.template negTrack_as<Tracks>();

    std::array<float, 2> etaDaughters = {std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 2> dcaDaughters = {std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()), std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ())};
    this->evaluateObservable(kDauDcaMin, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 2> clustersDaughters = {1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // daughter pid selections
    this->evaluateObservable(kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDaughTpcProton, posDaughter.tpcNSigmaPr());
    this->evaluateObservable(kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->evaluateObservable(kNegDaughTpcProton, negDaughter.tpcNSigmaPr());

    this->assembleBitmask();
  };

  template <typename T>
  bool checkLambdaHypothesis(T const& v0)
  {
    return (this->passesOptionalCut(kNegDaughTpcPion) && this->passesOptionalCut(kPosDaughTpcProton)) && // check PID for daughters
           (v0.mLambda() > mMassLambdaLowerLimit && v0.mLambda() < mMassLambdaUpperLimit) &&             // inside Lambda window
           (v0.mK0Short() < mMassK0shortLowerLimit || v0.mK0Short() > mMassK0shortUpperLimit);           // outside K0short window
  }

  template <typename T>
  bool checkAntiLambdaHypothesis(T const& v0)
  {
    return (this->passesOptionalCut(kPosDaughTpcPion) && this->passesOptionalCut(kNegDaughTpcProton)) && // check PID for daughters
           (v0.mAntiLambda() > mMassLambdaLowerLimit && v0.mAntiLambda() < mMassLambdaUpperLimit) &&     // inside AntiLambda window
           (v0.mK0Short() < mMassK0shortLowerLimit || v0.mK0Short() > mMassK0shortUpperLimit);           // outside K0short window
  }

  template <typename T>
  bool checkK0shortHypothesis(T const& v0)
  {
    return (this->passesOptionalCut(kPosDaughTpcPion) && this->passesOptionalCut(kNegDaughTpcPion)) && // check PID for daughters
           (v0.mK0Short() > mMassK0shortLowerLimit && v0.mK0Short() < mMassK0shortUpperLimit);         // inside K0short window
  }

 protected:
  float mMassK0shortLowerLimit = 0.f;
  float mMassK0shortUpperLimit = 99.f;

  float mMassLambdaLowerLimit = 0.f;
  float mMassLambdaUpperLimit = 99.f;
};
} // namespace v0selection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_
