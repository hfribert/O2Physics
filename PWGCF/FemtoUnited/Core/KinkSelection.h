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

/// \file KinkSelecion.h
/// \brief kink selections
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_KINKSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_KINKSELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/Configurable.h"
#include "Common/SystemOfUnits.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace kinkselection
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

#define KINK_DEFAULT_BITS                                                                                   \
  o2::framework::Configurable<std::vector<float>> kinkCosPAMin{"kinkCosPAMin", {0.99f}, "Minimum cosine of pointing angle"};         \
  o2::framework::Configurable<std::vector<float>> kinkTopoDcaMax{"kinkTopoDcaMax", {0.1f}, "Maximum kink topological DCA"};           \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};               \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"DauAbsEtaMax", {0.8f}, "Maximum absolute pseudorapidity for daughter track"}; \
  o2::framework::Configurable<std::vector<float>> dauDcaPvMin{"dauDcaPvMin", {0.05f}, "Minimum DCA of daughter from primary vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> mothDcaPvMax{"mothDcaPvMax", {0.1f}, "Maximum DCA of mother from primary vertex (cm)"};  \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"DauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter track"};

// derived selection bits for sigma
struct ConfSigmaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> chaDauTpcPion{"chaDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for charged daughter tracks"};
};

#undef KINK_DEFAULT_BITS

// base selection for analysis task for kinks
#define KINK_DEFAULT_SELECTIONS                                                                           \
  o2::framework::Configurable<int> pdgCode{"pdgCode", 3112, "Kink PDG code"};                             \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                 \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                               \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                            \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                             \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                              \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"}; \
  o2::framework::Configurable<float> massMin{"massMin", 1.1f, "Minimum invariant mass for Sigma"};      \
  o2::framework::Configurable<float> massMax{"massMax", 1.3f, "Maximum invariant mass for Sigma"};     \
  o2::framework::Configurable<o2::aod::femtodatatypes::KinkMaskType> mask{"mask", 6, "Bitmask for kink selection"};

// base selection for analysis task for sigmas
template <const char* Prefix>
struct ConfSigmaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Sigma (-1 for Sigma minus and +1 for Anti Sigma minus"};
};

#undef KINK_DEFAULT_SELECTIONS

// Define unique prefixes as constexpr string literals
constexpr const char PrefixSigmaSelection1[] = "SigmaSelection1";
using ConfSigmaSelection1 = ConfSigmaSelection<PrefixSigmaSelection1>;

/// The different selections for kinks
enum KinkSeles {
  kKinkCosPAMin,
  kKinkTopoDcaMax,
  kTransRadMin,
  kTransRadMax,

  kDauAbsEtaMax,
  kDauDcaPvMin,
  kMothDcaPvMax,
  kDauTpcClsMin,

  kChaDaughTpcPion,

  kKinkSelsMax
};

/// \class KinkSelection
/// \brief Cut class to contain and execute all cuts applied to tracks
class KinkSelection : public BaseSelection<float, o2::aod::femtodatatypes::KinkMaskType, kKinkSelsMax>
{
 public:
  KinkSelection() {}
  virtual ~KinkSelection() = default;

  template <o2::analysis::femtounited::modes::KinkCands KinkCands, typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (o2::analysis::femtounited::modes::isFlagSet(KinkCands, o2::analysis::femtounited::modes::KinkCands::kSigma)) {
      mMassSigmaLowerLimit = filter.massMinSigma.value;
      mMassSigmaUpperLimit = filter.massMaxSigma.value;
      this->addSelection(config.chaDauTpcPion.value, kChaDaughTpcPion, limits::kAbsUpperLimit, false, false);
    }

    this->addSelection(config.kinkCosPAMin.value, kKinkCosPAMin, limits::kLowerLimit, true, true);
    this->addSelection(config.kinkTopoDcaMax.value, kKinkTopoDcaMax, limits::kUpperLimit, true, true);
    this->addSelection(config.transRadMin.value, kTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMax.value, kTransRadMax, limits::kUpperLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaPvMin.value, kDauDcaPvMin, limits::kLowerLimit, true, true);
    this->addSelection(config.mothDcaPvMax.value, kMothDcaPvMax, limits::kUpperLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClsMin, limits::kLowerLimit, true, true);
  };

  template <typename KinkCand, typename Tracks>
  void applySelections(KinkCand const& kinkCand, Tracks const& tracks, const o2::aod::Collisions* collision)
  {
    this->reset();

    auto dauTrack = kinkCand.template trackDaug_as<Tracks>();
    auto mothTrack = kinkCand.template trackMoth_as<Tracks>();

    std::array<float, 3> decayCoordinates = {kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()};
    std::array<float, 3> pvCoordinates = {collision->posX(), collision->posY(), collision->posZ()};

    std::array<float, 3> decayVector = {decayCoordinates[0] - pvCoordinates[0], decayCoordinates[1] - pvCoordinates[1], decayCoordinates[2] - pvCoordinates[2]};
    std::array<float, 3> motherMomentum = {mothTrack.px(), mothTrack.py(), mothTrack.pz()};

    float decayVectorMag = std::sqrt(std::inner_product(decayVector.begin(), decayVector.end(), decayVector.begin(), 0.f));
    float motherMomentumMag = std::sqrt(std::inner_product(motherMomentum.begin(), motherMomentum.end(), motherMomentum.begin(), 0.f));
    float dotProduct = std::inner_product(decayVector.begin(), decayVector.end(), motherMomentum.begin(), 0.f);

    float cosPA = 0.f;
    if (decayVectorMag > 0.f && motherMomentumMag > 0.f) {
      cosPA = dotProduct / (decayVectorMag * motherMomentumMag);
    }
    this->evaluateObservable(kKinkCosPAMin, cosPA);

    this->evaluateObservable(kKinkTopoDcaMax, kinkCand.dcaKinkTopo());
    this->evaluateObservable(kTransRadMin, std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx()));
    this->evaluateObservable(kTransRadMax, std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx()));

    this->evaluateObservable(kDauAbsEtaMax, std::fabs(dauTrack.eta()));
    this->evaluateObservable(kDauDcaPvMin, kinkCand.dcaDaugPv());
    this->evaluateObservable(kMothDcaPvMax, kinkCand.dcaMothPv());
    this->evaluateObservable(kDauTpcClsMin, 1.f * dauTrack.tpcNClsFound());

    this->evaluateObservable(kChaDaughTpcPion, dauTrack.tpcNSigmaPi());

    this->assembleBitmask();
  };


  template <typename KinkCand, typename Tracks>
  bool checkSigmaHypothesis(KinkCand const& kinkCand, const o2::aod::Tracks* tracks)
  {
    if (!this->passesOptionalCut(kChaDaughTpcPion)) return false;
    return (kinkCand.mSigmaMinus() > mMassSigmaLowerLimit && kinkCand.mSigmaMinus() < mMassSigmaUpperLimit);
  }


 protected:
  float mMassSigmaLowerLimit = 0.f;
  float mMassSigmaUpperLimit = 99.f;
};
} // namespace kinkselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_KINKSELECTION_H_