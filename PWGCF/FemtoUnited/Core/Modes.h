// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Modes.h
/// \brief common modes
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_MODES_H_
#define PWGCF_FEMTOUNITED_CORE_MODES_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include <Rtypes.h>

#include <cstdint>
#include <type_traits>

namespace o2::analysis::femtounited
{
namespace modes
{

// check if flag is set
template <typename T>
constexpr bool isFlagSet(T value, T flag)
{
  using U = std::underlying_type_t<T>;
  return (static_cast<U>(value) & static_cast<U>(flag)) != 0;
}

// set a flag
template <typename T>
constexpr T setFlag(T value, T flag)
{
  using U = std::underlying_type_t<T>;
  return static_cast<T>(static_cast<U>(value) | static_cast<U>(flag));
}

// clear a flag
template <typename T>
constexpr T clearFlag(T value, T flag)
{
  using U = std::underlying_type_t<T>;
  return static_cast<T>(static_cast<U>(value) & ~static_cast<U>(flag));
}

// #define BIT(n) (1 << n)

enum class Mode : uint32_t {
  kANALYSIS = BIT(0),
  kQA = BIT(1),
  kMC = BIT(2),
  kANALYSIS_QA = kANALYSIS | kQA,
  kANALYSIS_MC = kANALYSIS | kMC,
  kANALYSIS_QA_MC = kANALYSIS | kQA | kMC,
};

enum class System : uint32_t {
  kPP = BIT(0),
  kPbPb = BIT(1),
  kMC = BIT(2),
  kRun3 = BIT(3),
  kRun2 = BIT(4),
  kNoCentCal = BIT(5),
  kPP_Run3 = kPP | kRun3,
  kPP_Run2 = kPP | kRun2,
  kPP_NoCentCal_Run3 = kPP | kRun3 | kNoCentCal,
  kPbPb_Run3 = kPbPb | kRun3,
  kPbPb_Run2 = kPbPb | kRun2,
};

enum class Track : o2::aod::femtodatatypes::TrackType {
  kPrimaryTrack = BIT(0),
  kLambdaDaugher = BIT(1),
  kCascadeBachelor = BIT(2),
  kResonanceDaughter = BIT(3),
};

enum class V0 : o2::aod::femtodatatypes::V0Type {
  kLambda = BIT(0),
  kAntiLambda = BIT(1),
  kK0short = BIT(2)
};

// enum of supported resonances
enum class TwoTrackResonace : o2::aod::femtodatatypes::TwoTrackResonaceType {
  kRho = BIT(0),
  kPhi = BIT(1),
  kKstar = BIT(2),
  kAntiKstar = BIT(3),
};

}; // namespace modes
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_MODES_H_
