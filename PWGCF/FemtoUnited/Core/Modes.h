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

#include <cstdint>
#include <type_traits>

namespace o2::analysis::femtounited
{
namespace modes
{

// Generic flag check
template <typename E>
constexpr bool isFlagSet(E value, E flag)
{
  using U = std::underlying_type_t<E>;
  return static_cast<U>(value) & static_cast<U>(flag);
}

enum class Mode : uint32_t {
  kANALYSIS = 0x1,
  kQA = 0x2,
  kMC = 0x4,
  kANALYSIS_QA = kANALYSIS | kQA,
  kANALYSIS_MC = kANALYSIS | kMC,
  kANALYSIS_QA_MC = kANALYSIS | kQA | kMC,
};

enum class System : uint32_t {
  kPP = 0x1,
  kPbPb = 0x2,
  kMC = 0x4,
  kRun3 = 0x8,
  kRun2 = 0x10,
  kNoCentCal = 0x20,
  kPP_Run3 = kPP | kRun3,
  kPP_Run2 = kPP | kRun2,
  kPP_NoCentCal_Run3 = kPP | kRun3 | kNoCentCal,
  kPbPb_Run3 = kPbPb | kRun3,
  kPbPb_Run2 = kPbPb | kRun2,
};

enum class TrackType : uint32_t {
  kPrimaryTrack,
  kLambdaDaugher,
  kCascadeBachelor,
  kTrackTypeMax
};

}; // namespace modes
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_MODES_H_
