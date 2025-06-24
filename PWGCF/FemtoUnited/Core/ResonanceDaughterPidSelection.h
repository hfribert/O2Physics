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

/// \file ResonanceDaughterPidSelection.h
/// \brief track pid selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_RESONANCEDAUGHTERPIDSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_RESONANCEDAUGHTERPIDSELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include "Framework/Configurable.h"

#include <cmath>

namespace o2::analysis::femtounited
{
namespace resonancedaughterpidselection
{

template <const char* Prefix>
struct ConfResonanceDaughterPidBits : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<std::vector<float>> posDaughItsPion{"posDaughItsPion", {}, "Maximum |nsimga_Pion| ITS for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcPion{"posDaughTpcPion", {}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTofPion{"posDaughTofPion", {}, "Maximum |nsimga_Pion| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpctofPion{"posDaughTpctofPion", {}, "Maximum |nsimga_Pion| TPCTOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughItsKaon{"posDaughItsKaon", {}, "Maximum |nsimga_Kaon| ITS for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpcKaon{"posDaughTpcKaon", {}, "Maximum |nsimga_Kaon| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTofKaon{"posDaughTofKaon", {}, "Maximum |nsimga_Kaon| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDaughTpctofKaon{"posDaughTpctofKaon", {}, "Maximum |nsimga_Kaon| TPCTOF for positive daughter tracks"};

  o2::framework::Configurable<std::vector<float>> negDaughItsPion{"negDaughItsPion", {}, "Maximum |nsimga_Pion| ITS for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcPion{"negDaughTpcPion", {}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTofPion{"negDaughTofPion", {}, "Maximum |nsimga_Pion| TOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpctofPion{"negDaughTpctofPion", {}, "Maximum |nsimga_Pion| TPCTOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughItsKaon{"negDaughItsKaon", {}, "Maximum |nsimga_Kaon| ITS for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpcKaon{"negDaughTpcKaon", {}, "Maximum |nsimga_Kaon| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTofKaon{"negDaughTofKaon", {}, "Maximum |nsimga_Kaon| TOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDaughTpctofKaon{"negDaughTpctofKaon", {}, "Maximum |nsimga_Kaon| TPCTOF for negative daughter tracks"};
};
constexpr const char PrefixRhoDaughterPidBits[] = "RhoDaughterPidBits";
constexpr const char PrefixPhiDaughterPidBits[] = "PhiDaughterPidBits";
constexpr const char PrefixKstarDaughterPidBits[] = "KstarDaughterPidBits";
using ConfRhoDaughterPidBits = ConfResonanceDaughterPidBits<PrefixRhoDaughterPidBits>;
using ConfPhiDaughterPidBits = ConfResonanceDaughterPidBits<PrefixPhiDaughterPidBits>;
using ConfKstarDaughterPidBits = ConfResonanceDaughterPidBits<PrefixKstarDaughterPidBits>;

enum ResonanceDaughterPidSels {
  kPosDaughItsPion,
  kPosDaughTpcPion,
  kPosDaughTofPion,
  kPosDaughTpctofPion,
  kPosDaughItsKaon,
  kPosDaughTpcKaon,
  kPosDaughTofKaon,
  kPosDaughTpctofKaon,

  kNegDaughItsPion,
  kNegDaughTpcPion,
  kNegDaughTofPion,
  kNegDaughTpctofPion,
  kNegDaughItsKaon,
  kNegDaughTpcKaon,
  kNegDaughTofKaon,
  kNegDaughTpctofKaon,

  kResonanceDaughterPidSelsMax
};

/// \class ResonanceDaughterPidSelectionPos
/// \brief Cut class to contain and execute all cuts applied to pid of vzerodaughters
class ResonanceDaughterPidSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackPidMaskType, kResonanceDaughterPidSelsMax>
{
 public:
  ResonanceDaughterPidSelection() {}
  virtual ~ResonanceDaughterPidSelection() = default;
  template <typename track>
  void applySelections(track const& posDaughter, track const& negDaughter)
  {
    this->reset();

    // apply bits for positive daughters
    this->setBitmaskForObservable(kPosDaughItsPion, posDaughter.itsNsigmaPi());
