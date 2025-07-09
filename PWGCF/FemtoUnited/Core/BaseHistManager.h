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

/// \file HistManager.h
/// \brief common structs for histogram managers
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_BASEHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_BASEHISTMANAGER_H_

#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace o2::analysis::femtounited
{
namespace basehistmanager
{

template <typename Enum>
struct HistInfo {
  Enum hist;
  o2::framework::HistType histtype;
  const char* histname;
  const char* histdesc;
};

template <const char* prefix>
class BaseHistManager
{
 protected:
  o2::framework::HistogramRegistry* mRegistry = nullptr;

  // Constructor
  BaseHistManager(o2::framework::HistogramRegistry* registry) : mRegistry(registry) {}

  template <typename EnumType, typename TableType>
  static constexpr const char* getName(EnumType hist, const TableType& table)
  {
    for (const auto& entry : table) {
      if (entry.hist == hist)
        return entry.histname;
    }
    return "";
  }

  template <typename EnumType, typename TableType>
  static constexpr const char* getDesc(EnumType hist, const TableType& table)
  {
    for (const auto& entry : table) {
      if (entry.hist == hist)
        return entry.histdesc.data();
    }
    return "";
  }

  template <typename EnumType, typename TableType>
  static constexpr o2::framework::HistType getType(EnumType hist, const TableType& table)
  {
    for (const auto& entry : table) {
      if (entry.hist == hist)
        return entry.histtype;
    }
    return o2::framework::kUndefinedHist;
  }

  template <typename EnumType, typename TableType>
  void add(const char* dir, EnumType hist, const TableType& table, const std::vector<o2::framework::AxisSpec>& axes)
  {
    mRegistry->add(std::string(prefix) + dir + std::string(getName(hist, table)), getDesc(hist, table), getType(hist, table), {axes});
  }

  template <typename EnumType, typename TableType, typename... Args>
  void fill(const char* dir, EnumType hist, const TableType& table, Args&&... args)
  {
    mRegistry->fill(HIST(prefix) + HIST(dir) + HIST(getName(hist, table)), std::forward<Args>(args)...);
  }
};
} // namespace basehistmanager
} // namespace o2::analysis::femtounited
