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

#ifndef PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include "Framework/Configurable.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace cascadeselection
{

struct ConfCascadeFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massXiMin{"massXiMin", 1.25f, "Minimum Xi mass"};
  o2::framework::Configurable<float> massXiMax{"massXiMax", 1.39f, "Maximum Xi mass"};
  o2::framework::Configurable<float> rejectMassXiMin{"rejectmassXiMin", 1.25f, "Reject Minimum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> rejectMassXiMax{"rejectMassXiMax", 1.39f, "Rejection Maximum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> massOmegaMin{"massOmegaMin", 1.6f, "Minimum Omega mass"};
  o2::framework::Configurable<float> massOmegaMax{"massOmegaMax", 1.73f, "Maximum Omega mass"};
  o2::framework::Configurable<float> rejectMassOmegaMin{"rejectMassOmegaMin", 1.6f, "Reject minimum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> rejectMassOmegaMax{"rejectMassOmegaMax", 1.73f, "Reject maximum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> massLambdaMin{"massLambdaMin", 1.6f, "Minimum Lambda mass"};
  o2::framework::Configurable<float> massLambdaMax{"massLambdaMax", 1.73f, "Maximum Lambda mass"};
};

struct ConfCascadeBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeBits");
  // cascade selections
  o2::framework::Configurable<std::vector<float>> cascadeCpaMin{"cascadeCpaMin", {0.99f}, "Minimum cosine of pointing angle"};
  o2::framework::Configurable<std::vector<float>> cascadeTransRadMin{"cascadeTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};
  o2::framework::Configurable<std::vector<float>> cascadeDcaDauMax{"cascadeDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"};
  // lambda selections
  o2::framework::Configurable<std::vector<float>> lambdaCpaMin{"lambdaCpaMin", {0.99f}, "Minimum cosine of pointing angle"};
  o2::framework::Configurable<std::vector<float>> lambdaTransRadMin{"lambdaTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};
  o2::framework::Configurable<std::vector<float>> lambdaDcaDauMax{"lambdaDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> lambdaDcaToPvMax{"lambdaDcaToPvMin", {0.3f}, "Maximum DCA between the daughters at decay vertex (cm)"};
  // bachelor/daughter selections
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"DauAbsEtaMax", {0.8f}, "Minimum DCA of the daughters from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"DauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"DauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};
  // bachelor/daughter pid selections
  o2::framework::Configurable<std::vector<float>> bachelorTpcPion{"bachelorTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for bachelor tracks"};
  o2::framework::Configurable<std::vector<float>> bachelorTpcKaon{"bachelorTpcKaon", {5.f}, "Maximum |nsimga_Kaon| TPC for bachelor tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTpc{"posDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpc{"negDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for negative daughter tracks"};
};

template <const char* Prefix>
struct ConfCascadeSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for Cascade"};
  o2::framework::Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for Cascade"};
  o2::framework::Configurable<float> antiMassMin{"antiMassMin", 0.f, "Minimum invariant mass for AntiCascade"};
  o2::framework::Configurable<float> antiMassMax{"antiMassMax", 999.f, "Maximum invariant mass for AntiCascade"};
  o2::framework::Configurable<o2::aod::femtodatatypes::CascadeMaskType> mask{"mask", 6, "Bitmask for cascade selection"};
};

// Define unique prefixes as constexpr string literals
constexpr const char PrefixCascadeSelection1[] = "CascadeSelection1";
constexpr const char PrefixCascadeSelection2[] = "CascadeSelection2";
constexpr const char PrefixCascadeSelection3[] = "CascadeSelection3";

// Instantiate different instances with unique prefixes
using ConfCascadeSelection1 = ConfCascadeSelection<PrefixCascadeSelection1>;
using ConfCascadeSelection2 = ConfCascadeSelection<PrefixCascadeSelection2>;
using ConfCascadeSelection3 = ConfCascadeSelection<PrefixCascadeSelection3>;

/// The different selections this task is capable of doing
enum CascadeSels {
  // selections for cascades
  kCascadeCpaMin,      ///< Min. CPA (cosine pointing angle)
  kCascadeDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kCascadeTransRadMin, ///< max. transverse radius

  // selection for lambda daughter
  kLambdaCpaMin,      ///< Min. DCA of the positive daughers at primary vertex
  kLambdaDcaDauMax,   ///< TPC Pion PID for positive daughter
  kLambdaTransRadMin, ///< Min. number of TPC clusters of positive daughter
  kLambdaDcaToPvMin,  ///< TPC Pion PID for positive daughter

  // selection for bachelor/daugthers
  kDauAbsEtaMax, ///< Min. DCA of the positive daughers at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kDauDcaMin,    ///< TPC Pion PID for negative daughter

  // PID selection for cascade bachelor
  kBachelorTpcPion, ///< TPC Proton PID for negative daughter
  kBachelorTpcKaon, ///< TPC Proton PID for negative daughter
                    ///
  // PID selection for lambda daughers
  kPosDauTpc, ///< TPC Proton PID for negative daughter
  kNegDauTpc, ///< TPC Proton PID for negative daughter

  kCascadeSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class CascadeSelection : public BaseSelection<float, o2::aod::femtodatatypes::CascadeMaskType, kCascadeSelsMax>
{
 public:
  CascadeSelection() {}
  virtual ~CascadeSelection() = default;

  template <typename T1, typename T2, typename T3>
  void applySelections(T1 const& cascade, T2 const& /*tracks*/, T3 const& col)
  {
    this->reset();
    // cascade selections
    this->evaluateObservable(kCascadeCpaMin, cascade.casccosPA(col.posX(), col.posY(), col.posZ()));
    this->evaluateObservable(kCascadeDcaDaughMax, cascade.dcacascdaughters());
    this->evaluateObservable(kCascadeTransRadMin, cascade.cascradius());

    // lambda selection
    this->evaluateObservable(kLambdaCpaMin, cascade.v0cosPA(col.posY(), col.posY(), col.posZ()));
    this->evaluateObservable(kLambdaDcaDauMax, cascade.dcaV0daughters());
    this->evaluateObservable(kLambdaTransRadMin, cascade.v0radius());
    this->evaluateObservable(kLambdaDcaToPvMin, cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));

    auto bachelor = cascade.template bachelor_as<T2>();
    auto posDaughter = cascade.template posTrack_as<T2>();
    auto negDaughter = cascade.template negTrack_as<T2>();

    // daughter selections
    std::array<float, 3> etaDaughters = {std::fabs(bachelor.eta()), std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 3> dcaDaughters = {std::hypot(bachelor.dcaXY(), bachelor.dcaZ()),
                                         std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()),
                                         std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ())};
    this->evaluateObservable(kDauAbsEtaMax, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 3> clustersDaughters = {1.f * bachelor.tpcNClsFound(), 1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauAbsEtaMax, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // bachelor pid selection
    // check both pion and kaon PID for xi and omega
    this->evaluateObservable(kBachelorTpcPion, bachelor.tpcNSigmaPi());
    this->evaluateObservable(kBachelorTpcKaon, bachelor.tpcNSigmaKa());

    // depending on the charge, we check lambda and antilambda hypothesis
    if (cascade.sign() < 0) {
      this->evaluateObservable(kPosDauTpc, bachelor.tpcNSigmaPr());
      this->evaluateObservable(kNegDauTpc, bachelor.tpcNSigmaPi());
    } else if (cascade.sign() > 0) {
      this->evaluateObservable(kPosDauTpc, bachelor.tpcNSigmaPi());
      this->evaluateObservable(kNegDauTpc, bachelor.tpcNSigmaPr());
    } else {
      LOG(warn) << "Encountered Cascade candidate with 0 charge";
    }

    this->assembleBitmask();
  };

  void setXiMassLimits(float lower, float upper)
  {
    mXiMassLowerLimit = lower;
    mXiMassUpperLimit = upper;
  };

  void setOmegaMassLimits(float lower, float upper)
  {
    mOmegaMassLowerLimit = lower;
    mOmegaMassUpperLimit = upper;
  };

  template <typename T>
  bool checkXiHypothesis(T const& cascade)
  {
    return (this->passesOptionalCut(kBachelorTpcPion) && this->passesOptionalCut(kPosDauTpc) && this->passesOptionalCut(kNegDauTpc)) &&
           (mXiMassLowerLimit < cascade.mXi() && mXiMassUpperLimit > cascade.mXi());
  };

  template <typename T>
  bool checkOmegaHypothesis(T const& cascade)
  {
    return (this->passesOptionalCut(kBachelorTpcKaon) && this->passesOptionalCut(kPosDauTpc) && this->passesOptionalCut(kNegDauTpc)) &&
           (mXiMassLowerLimit < cascade.mOmega() && mXiMassUpperLimit > cascade.mOmega());
  };

 protected:
  float mXiMassLowerLimit = 0.f;
  float mXiMassUpperLimit = 999.f;

  float mOmegaMassLowerLimit = 0.f;
  float mOmegaMassUpperLimit = 999.f;
};
} // namespace cascadeselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
