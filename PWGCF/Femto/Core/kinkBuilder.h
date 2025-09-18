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

/// \file kinkBuilder.h
/// \brief kink builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTO_CORE_KINKBUILDER_H_
#define PWGCF_FEMTO_CORE_KINKBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace kinkbuilder
{

// filters applied in the producer task
struct ConfKinkFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMinSigma{"massMinSigma", 1.1f, "Minimum mass for Sigma hypothesis"};
  o2::framework::Configurable<float> massMaxSigma{"massMaxSigma", 1.3f, "Maximum mass for Sigma hypothesis"};
};

// selections bits for all kinks
#define KINK_DEFAULT_BITS                                                                                                                     \
  o2::framework::Configurable<std::vector<float>> kinkTopoDcaMax{"kinkTopoDcaMax", {2.0f}, "Maximum kink topological DCA"};                   \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                       \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"DauAbsEtaMax", {0.8f}, "Maximum absolute pseudorapidity for daughter track"}; \
  o2::framework::Configurable<std::vector<float>> dauDcaPvMin{"dauDcaPvMin", {0.0f}, "Minimum DCA of daughter from primary vertex (cm)"};     \
  o2::framework::Configurable<std::vector<float>> mothDcaPvMax{"mothDcaPvMax", {1.0f}, "Maximum DCA of mother from primary vertex (cm)"};     \
  o2::framework::Configurable<std::vector<float>> alphaAPMax{"alphaAPMax", {0.0f}, "Maximum Alpha_AP for Sigma candidates"};                  \
  o2::framework::Configurable<std::vector<float>> qtAPMin{"qtAPMin", {0.15f}, "Minimum qT_AP for Sigma candidates"};

// derived selection bits for sigma
struct ConfSigmaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> chaDauTpcPion{"chaDauTpcPion", {5.f}, "Maximum |nsigma_Pion| TPC for charged daughter tracks"};
};

#undef KINK_DEFAULT_BITS

// base selection for analysis task for kinks
#define KINK_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                         \
  o2::framework::Configurable<int> pdgCode{"pdgCode", defaultPdgCode, "Kink PDG code"};                 \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                 \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                               \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                            \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                             \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                              \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"}; \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Sigma"};      \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Sigma"};      \
  o2::framework::Configurable<o2::aod::femtodatatypes::KinkMaskType> mask{"mask", 0, "Bitmask for kink selection"};

// base selection for analysis task for sigmas
template <const char* Prefix>
struct ConfSigmaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3112)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Sigma mother track (e.g. -1 for Sigma- or +1 for AntiSigma-)"};
};

#undef KINK_DEFAULT_SELECTIONS

constexpr const char PrefixSigmaSelection1[] = "SigmaSelection1";
using ConfSigmaSelection1 = ConfSigmaSelection<PrefixSigmaSelection1>;

/// The different selections for kinks
enum KinkSeles {
  kKinkTopoDcaMax,
  kTransRadMin,
  kTransRadMax,

  kDauAbsEtaMax,
  kDauDcaPvMin,
  kMothDcaPvMax,

  kChaDaughTpcPion,

  kAlphaAPMax,
  kQtAPMin,

  kKinkSelsMax
};

const char kinkSelsName[] = "Kink selection object";
const std::unordered_map<KinkSeles, std::string> kinkSelsToStrings = {
  {kKinkTopoDcaMax, "kinkTopoDcaMax"},
  {kTransRadMin, "transRadMin"},
  {kTransRadMax, "transRadMax"},
  {kDauAbsEtaMax, "dauAbsEtaMax"},
  {kDauDcaPvMin, "dauDcaPvMin"},
  {kMothDcaPvMax, "mothDcaPvMax"},
  {kChaDaughTpcPion, "chaDauTpcPion"},
  {kAlphaAPMax, "alphaAPMax"},
  {kQtAPMin, "qtAPMin"}
};

/// \class KinkCuts
/// \brief Cut class to contain and execute all cuts applied to kinks
template <modes::Kink kinkType>
class KinkSelection : public BaseSelection<float, o2::aod::femtodatatypes::KinkMaskType, kKinkSelsMax>
{
 public:
  KinkSelection() {}
  virtual ~KinkSelection() = default;

  template <typename T1, typename T2>
  void configure(T1& config, T2& filter)
  {
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value; 
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;
    
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      mMassSigmaLowerLimit = filter.massMinSigma.value;
      mMassSigmaUpperLimit = filter.massMaxSigma.value;
      // Only add PID selection if we need it - will be checked at runtime
      this->addSelection(config.chaDauTpcPion.value, kChaDaughTpcPion, limits::kAbsUpperLimit, true, true);
    }

    this->addSelection(config.kinkTopoDcaMax.value, kKinkTopoDcaMax, limits::kUpperLimit, true, true);
    this->addSelection(config.transRadMin.value, kTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMax.value, kTransRadMax, limits::kUpperLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaPvMin.value, kDauDcaPvMin, limits::kLowerLimit, true, true);
    this->addSelection(config.mothDcaPvMax.value, kMothDcaPvMax, limits::kUpperLimit, true, true);
    this->addSelection(config.alphaAPMax.value, kAlphaAPMax, limits::kUpperLimit, true, true);
    this->addSelection(config.qtAPMin.value, kQtAPMin, limits::kLowerLimit, true, true);
  };

  template <typename T1, typename T2>
  void applySelections(T1 const& kinkCand, T2 const& /*tracks*/)
  {
    this->reset();
    // kink selections
    std::array<float, 3> momMother = {kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
    std::array<float, 3> momDaughter = {kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

    // Alpha_AP
    std::array<float, 3> momMissing = {momMother[0] - momDaughter[0], momMother[1] - momDaughter[1], momMother[2] - momDaughter[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momDaughter.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    float alphaAP = (lQlP + lQlN != 0.f) ? (lQlP - lQlN) / (lQlP + lQlN) : 0.f;
    this->evaluateObservable(kAlphaAPMax, alphaAP);

    // qT_AP
    float dp = lQlP;
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momDaughter.begin(), momDaughter.end(), momDaughter.begin(), 0.f);
    float qtAP = std::sqrt(std::max(0.f, p2A - dp * dp / p2V0));
    this->evaluateObservable(kQtAPMin, qtAP);

    this->evaluateObservable(kKinkTopoDcaMax, kinkCand.dcaKinkTopo());
    
    // float transRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());
    this->evaluateObservable(kTransRadMin, kinkCand.transRadius());
    this->evaluateObservable(kTransRadMax, kinkCand.transRadius());

    this->evaluateObservable(kDauAbsEtaMax, std::fabs(kinkCand.eta()));
    this->evaluateObservable(kDauDcaPvMin, std::abs(kinkCand.dcaDaugPv()));
    this->evaluateObservable(kMothDcaPvMax, std::abs(kinkCand.dcaMothPv()));

    // Conditionally evaluate PID information if available
    if constexpr (requires { kinkCand.nSigmaTPCPi(); }) {
      this->evaluateObservable(kChaDaughTpcPion, std::abs(kinkCand.nSigmaTPCPi()));
    }

    this->assembleBitmask();
  };

  template <typename T>
  bool checkFilters(const T& kink) const
  {
    return ((kink.ptMoth() > mPtMin && kink.ptMoth() < mPtMax) &&
            (kink.eta() > mEtaMin && kink.eta() < mEtaMax) &&
            (kink.phi() > mPhiMin && kink.phi() < mPhiMax));
  }

  template <typename T>
  bool checkHypothesis(T const& kinkCand) const
  {
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      return (kinkCand.mSigmaMinus() > mMassSigmaLowerLimit && kinkCand.mSigmaMinus() < mMassSigmaUpperLimit);
    }
    return false;
  }
    
 protected:
  float mMassSigmaLowerLimit = 1.15f;
  float mMassSigmaUpperLimit = 1.25f;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
};

struct KinkBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FSigmas> producedSigmas;
  o2::framework::Produces<o2::aod::FSigmaMasks> producedSigmaMasks;
  o2::framework::Produces<o2::aod::FSigmaExtras> producedSigmaExtras;
};

struct ConfKinkTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkTables");
  o2::framework::Configurable<int> produceSigmas{"produceSigmas", -1, "Produce Sigmas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaMasks{"produceSigmaMasks", -1, "Produce SigmaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaExtras{"produceSigmaExtras", -1, "Produce SigmaExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::Kink kinkType>
class KinkBuilder
{
 public:
  KinkBuilder() {}
  virtual ~KinkBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(T1& config, T2& filter, T3& table, T4& initContext)
  {
    kinkSelection.configure(config, filter);
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      produceSigmas = utils::enableTable("FSigmas_001", table.produceSigmas.value, initContext);
      produceSigmaMasks = utils::enableTable("FSigmaMasks_001", table.produceSigmaMasks.value, initContext);
      produceSigmaExtras = utils::enableTable("FSigmaExtras_001", table.produceSigmaExtras.value, initContext);
    }
    
    if (produceSigmas || produceSigmaMasks || produceSigmaExtras) {
      mFillAnyTable = true;
      kinkSelection.printSelections(kinkSelsName, kinkSelsToStrings);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void fillKinks(T1& collisionProducts, T2& trackProducts, T3& kinkProducts, T4 const& kinks, T5 const& tracks, T6& trackBuilder, T7& indexMap)
  {
    if (!mFillAnyTable) {
      return;
    }
    
    int64_t daughterIndex = 0;
    
    for (const auto& kink : kinks) {
      if (!kinkSelection.checkFilters(kink)) {
        continue;
      }
      
      kinkSelection.applySelections(kink, tracks);
      
      bool passesAllSelections = kinkSelection.passesAllRequiredSelections();
      if (!passesAllSelections) {
        continue;
      }
      
      if (!kinkSelection.checkHypothesis(kink)) {
        continue;
      }
      
      auto daughter = kink.template trackDaug_as<T5>();
      daughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kKinkDaughter>(daughter, trackProducts, collisionProducts, indexMap);
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
        fillSigma(collisionProducts, kinkProducts, kink, daughterIndex);
      }
    }
  }
  template <typename T1, typename T2, typename T3>
  void fillSigma(T1& collisionProducts, T2& kinkProducts, T3 const& kink, int daughterIndex)
  {
    float mass;
    mass = kink.mSigmaMinus();
    if (produceSigmas) {
      kinkProducts.producedSigmas(collisionProducts.producedCollision.lastIndex(),
                                 kink.mothSign() * kink.ptMoth(),
                                 kink.eta(),
                                 kink.phi(),
                                 mass,
                                 daughterIndex);
    }
    if (produceSigmaMasks) {
      kinkProducts.producedSigmaMasks(kinkSelection.getBitmask());
    }
    if (produceSigmaExtras) {
      kinkProducts.producedSigmaExtras(
        kink.kinkAngle(),
        kink.dcaDaugPv(),
        kink.dcaMothPv(),
        kink.xDecVtx(),
        kink.yDecVtx(),
        kink.zDecVtx(),
        kink.transRadius());
    }
  }


  bool fillAnyTable() { return mFillAnyTable; }

 private:
  KinkSelection<kinkType> kinkSelection;
  bool mFillAnyTable = false;
  bool produceSigmas = false;
  bool produceSigmaMasks = false;
  bool produceSigmaExtras = false;
};
} // namespace kinkbuilder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKBUILDER_H_