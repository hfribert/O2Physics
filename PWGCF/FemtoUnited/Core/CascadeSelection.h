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

/// \file CascadeSelection.h
/// \brief cascade selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

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
  o2::framework::Configurable<float> rejectMassXiMin{"rejectMassXiMin", 1.25f, "Reject Minimum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> rejectMassXiMax{"rejectMassXiMax", 1.39f, "Rejection Maximum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> massOmegaMin{"massOmegaMin", 1.6f, "Minimum Omega mass"};
  o2::framework::Configurable<float> massOmegaMax{"massOmegaMax", 1.73f, "Maximum Omega mass"};
  o2::framework::Configurable<float> rejectMassOmegaMin{"rejectMassOmegaMin", 1.6f, "Reject minimum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> rejectMassOmegaMax{"rejectMassOmegaMax", 1.73f, "Reject maximum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> massLambdaMin{"massLambdaMin", 1.6f, "Minimum Lambda mass"};
  o2::framework::Configurable<float> massLambdaMax{"massLambdaMax", 1.73f, "Maximum Lambda mass"};
};

#define CASCADE_DEFAULT_BITS                                                                                                                               \
  o2::framework::Configurable<std::vector<float>> cascadeCpaMin{"cascadeCpaMin", {0.99f}, "Minimum cosine of pointing angle"};                             \
  o2::framework::Configurable<std::vector<float>> cascadeTransRadMin{"cascadeTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> cascadeDcaDauMax{"cascadeDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> lambdaCpaMin{"lambdaCpaMin", {0.99f}, "Minimum cosine of pointing angle"};                               \
  o2::framework::Configurable<std::vector<float>> lambdaTransRadMin{"lambdaTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                        \
  o2::framework::Configurable<std::vector<float>> lambdaDcaDauMax{"lambdaDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"};   \
  o2::framework::Configurable<std::vector<float>> lambdaDcaToPvMin{"lambdaDcaToPvMin", {0.3f}, "Maximum DCA between the daughters at decay vertex (cm)"};  \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Minimum DCA of the daughters from primary vertex (cm)"};           \
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"dauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};                \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};    \
  o2::framework::Configurable<std::vector<float>> posDauTpc{"posDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTpc{"negDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for negative daughter tracks"};

struct ConfXiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcPion{"bachelorTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for bachelor tracks"};
};

struct ConfOmegaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcKaon{"bachelorTpcKaon", {5.f}, "Maximum |nsimga_Kaon| TPC for bachelor tracks"};
};

#undef CASCADE_DEFAULT_BITS

template <const char* Prefix>
struct ConfCascadeSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> pdgCode{"pdgCode", 3312, "Track PDG code"};
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Lambda (+1 for Lambda and -1 for Antilambda"};
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for Cascade"};
  o2::framework::Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for Cascade"};
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
  kLambdaCpaMin,      ///< Min. DCA of the lambda daughers at primary vertex
  kLambdaDcaDauMax,   ///< TPC PID for daughters (Pion/Proton)
  kLambdaTransRadMin, ///< Min. number of TPC clusters of daughter
  kLambdaDcaToPvMin,  ///< Min. DCA to primary vertex of daughter lambda

  // selection for bachelor/daugthers
  kDauAbsEtaMax, ///< Min. DCA of the daughers/bachelor at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of daughters/bachelor
  kDauDcaMin,    ///< TPC Pion PID for negative daughter

  // PID selection for cascade bachelor
  kBachelorTpcPion, ///< TPC Pion PID for bachelor
  kBachelorTpcKaon, ///< TPC Kaon PID for bachelor
                    ///
  // PID selection for lambda daughers
  kPosDauTpc, ///< TPC PID for positive daughter
  kNegDauTpc, ///< TPC PID for negative daughter

  kCascadeSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class CascadeSelection : public BaseSelection<float, o2::aod::femtodatatypes::CascadeMaskType, kCascadeSelsMax>
{
 public:
  CascadeSelection() {}
  virtual ~CascadeSelection() = default;

  template <o2::analysis::femtounited::modes::Cascade cascade, typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(cascade, o2::analysis::femtounited::modes::Cascade::kXi)) {
      mXiMassLowerLimit = filter.massXiMin.value;
      mXiMassLowerLimit = filter.massXiMax.value;
      mOmegaMassLowerLimit = filter.rejectMassOmegaMin.value;
      mOmegaMassUpperLimit = filter.rejectMassOmegaMax.value;
      mType = o2::analysis::femtounited::modes::Cascade::kXi;
      this->addSelection(config.bachelorTpcPion.value, kBachelorTpcPion, limits::kAbsUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(cascade, o2::analysis::femtounited::modes::Cascade::kOmega)) {
      mOmegaMassLowerLimit = filter.massOmegaMin.value;
      mOmegaMassUpperLimit = filter.massOmegaMax.value;
      mXiMassLowerLimit = filter.rejectMassXiMin.value;
      mXiMassLowerLimit = filter.rejectMassXiMax.value;
      mType = o2::analysis::femtounited::modes::Cascade::kOmega;
      this->addSelection(config.bachelorTpcKaon.value, kBachelorTpcKaon, limits::kAbsUpperLimit, false, false);
    }

    this->addSelection(config.cascadeCpaMin.value, kCascadeCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.cascadeTransRadMin.value, kCascadeTransRadMin, limits::kUpperLimit, true, true);
    this->addSelection(config.cascadeDcaDauMax.value, kCascadeDcaDaughMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaCpaMin.value, kLambdaCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaTransRadMin.value, kLambdaTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaDcaDauMax.value, kLambdaDcaDauMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaDcaToPvMin.value, kLambdaDcaToPvMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaMin.value, kDauDcaMin, limits::kAbsLowerLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.posDauTpc.value, kPosDauTpc, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.negDauTpc.value, kNegDauTpc, limits::kAbsUpperLimit, false, false);
  };

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

  template <typename T>
  bool checkHypothesis(T const& cascade)
  {

    if (mType == o2::analysis::femtounited::modes::Cascade::kXi) {
      return (this->passesOptionalCut(kBachelorTpcPion) && this->passesOptionalCut(kPosDauTpc) && this->passesOptionalCut(kNegDauTpc)) && // check PID of bachelor and lambda daughters
             (mXiMassLowerLimit < cascade.mXi() && mXiMassUpperLimit > cascade.mXi()) &&                                                  // inside xi mass window
             (cascade.mOmega() < mOmegaMassLowerLimit || cascade.mOmega() > mOmegaMassUpperLimit);                                        // outside omega mass window
    }

    if (mType == o2::analysis::femtounited::modes::Cascade::kOmega) {
      return (this->passesOptionalCut(kBachelorTpcKaon) && this->passesOptionalCut(kPosDauTpc) && this->passesOptionalCut(kNegDauTpc)) && // check PID of bachelor and lambda daughters
             (mXiMassLowerLimit < cascade.mOmega() && mXiMassUpperLimit > cascade.mOmega()) &&                                            // inside omega mass window
             (cascade.mXi() < mXiMassLowerLimit || cascade.mXi() > mXiMassUpperLimit);                                                    // outside xi mass window
    }
    return false;
  }

 protected:
  float mXiMassLowerLimit = 0.f;
  float mXiMassUpperLimit = 999.f;

  float mOmegaMassLowerLimit = 0.f;
  float mOmegaMassUpperLimit = 999.f;

  o2::analysis::femtounited::modes::Cascade mType;
};
} // namespace cascadeselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
