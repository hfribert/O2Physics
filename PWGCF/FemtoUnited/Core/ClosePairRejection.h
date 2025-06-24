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

/// \file ClosePairRejection.h
/// \brief Definition of ClosePairRejection class
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_
#define PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/FemtoUtils.h"
#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace o2::analysis::femtounited
{
namespace closepairrejection
{
// enum for track histograms
enum CprHist {
  // kinemtics
  kAverage,
  kRadius0,
  kRadius1,
  kCprHistogramLast
};

// masks to extract charge
constexpr o2::aod::femtodatatypes::TrackMaskType kSignPlusMask = 2;
constexpr o2::aod::femtodatatypes::TrackMaskType kSignMinusMask = 1;
// tpc radii for computing phistar
constexpr int kNradii = 9;
constexpr std::array<float, kNradii> kTpcRadius = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

// directory names
constexpr char PrefixTrackTrackSe[] = "TrackTrack/CPR/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/CPR/ME/";
constexpr char PrefixTrackVzeroSe[] = "TrackVzero/CPR/SE/";
constexpr char PrefixTrackVzeroMe[] = "TrackVzero/CPR/ME/";
constexpr std::string_view AnalysisDir = "Analysis/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CprHist>, kCprHistogramLast> HistTable = {
  {{kAverage, o2::framework::kTH2F, "hAverage", "Average: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius0, o2::framework::kTH2F, "hRadius1", "Radius 1: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius1, o2::framework::kTH2F, "hRadius2", "Radius 2: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"}}};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class ClosePairRejection
{
 public:
  /// Destructor
  virtual ~ClosePairRejection() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {Specs[kAverage]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {Specs[kRadius0]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {Specs[kRadius1]});
    }

    // if constexpr (isModeSet(mode, modes::Mode::kQA)) {
    // std::string qaDir = std::string(prefix) + std::string(QaDir);
    // }
  }

  void setMagField(float magField) { mMagField = magField; }
  void setLimits(float detaMax, float dphistarMax)
  {
    mDetaMax = detaMax;
    mDphistarMax = dphistarMax;
  };

  template <typename T1, typename T2>
  void setPair(T1 const& track1, T2 const& track2)
  {
    // set deta
    mDeta = track1.eta() - track2.eta();
    // set dphi at primary vertex
    mDphi = track1.phi() - track2.phi();

    // get charge of the tracks
    int chargeTrack1 = 0.;
    int chargeTrack2 = 0.;

    if ((track1.trackMask() & kSignPlusMask) == kSignPlusMask) {
      chargeTrack1 = 1;
    } else if ((track1.trackMask() & kSignMinusMask) == kSignMinusMask) {
      chargeTrack1 = -1;
    } else {
      LOG(warn) << "ClosePairRejection was not able to determine the charge of the track!";
    }
    if ((track2.trackMask() & kSignPlusMask) == kSignPlusMask) {
      chargeTrack2 = 1;
    } else if ((track2.trackMask() & kSignMinusMask) == kSignMinusMask) {
      chargeTrack2 = -1;
    } else {
      LOG(warn) << "ClosePairRejection was not able to determine the charge of the track!";
    }

    // compute dphistar at different TPC radii
    for (int i = 0; i < kNradii; i++) {
      mDphistar.at(i) = utils::dphistar(mMagField, kTpcRadius.at(i), chargeTrack1, track1.pt(), track1.phi()) -
                        utils::dphistar(mMagField, kTpcRadius.at(i), chargeTrack2, track2.pt(), track2.phi());
      ;
    }
    // get average
    mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / mDphistar.size();
  }

  bool isClosePair()
  {
    if (std::pow(mAverageDphistar / mDphistarMax, 2) + std::pow(mDeta / mDetaMax, 2) < 1.) {
      return true;
    } else {
      return false;
    }
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kAverage, HistTable)), mDeta, mAverageDphistar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius0, HistTable)), mDeta, mDphistar.at(0));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius1, HistTable)), mDeta, mDphistar.at(1));
    }

    // if constexpr (isModeSet(mode, modes::Mode::kQA)) {
    //  mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsDcaz, HistTable)), track.pt(), track.dcaZ());
    // }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
  float mMagField = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;
  float mDeta = 0;
  float mDphi = 0;
  float mAverageDphistar = 0;
  std::array<float, kNradii> mDphistar{0};
};
}; // namespace closepairrejection
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_
