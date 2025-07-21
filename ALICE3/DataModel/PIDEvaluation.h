#ifndef PIDEVALUATION_H
#define PIDEVALUATION_H

#include <string>
#include <vector>
#include <array>
#include <map>
#include <numeric>
#include <cmath>

#include <TProfile.h>
#include <TH1F.h>
#include <TVector3.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/McTrackLabels.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"

using namespace o2;
using namespace o2::framework;

class PidEvaluation : public AnalysisTask {

public:
  Configurable<float> mNsigmaCut{"nSigmaCut", 3.0f, "N-sigma cut for PID identification"};
  Configurable<bool> mDoTOF{"doTOF", true, "Include TOF PID in baseline/combined PID"};
  Configurable<bool> mDoRICH{"doRICH", true, "Include RICH PID in baseline/combined PID"};
  Configurable<bool> mDoTracker{"doTracker", true, "Include Tracker ToT PID in combined PID (always requires a baseline)"};

  Configurable<int> mNumLogBins{"numLogBins", 200, "Number of logarithmic momentum bins"};
  Configurable<float> mPMin{"pMin", 0.05f, "Minimum momentum for histograms"};
  Configurable<float> mPMax{"pMax", 10.0f, "Maximum momentum for histograms"};

  HistogramRegistry histos{"PidEvalHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;

  static constexpr int kNumHypothesisParticles = 9;
  std::array<int, kNumHypothesisParticles> mHypothesisPdgCodes = {
      11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040
  };

  std::array<std::string, kNumHypothesisParticles> mParticleNames = {
      "Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "Helium3", "Alpha"
  };

  std::array<std::string, kNumHypothesisParticles> mParticleNamesPretty = {
      "#it{e}", "#it{#mu}", "#it{#pi}", "#it{K}", "#it{p}", "#it{d}", "#it{t}", "#it{^{3}He}", "#it{^{4}He}"
  };

  std::vector<double> mLogBins;
  int mNumBins;

  std::vector<std::vector<std::vector<double>>> mPurityNumerator;
  std::vector<std::vector<std::vector<double>>> mPurityDenominator;
  std::vector<std::vector<double>> mEfficiencyDenominator;

  std::map<int, int> mPdgToIndexMap;

  void init(o2::framework::InitContext& ctx) override;
  void process(soa::Join<
                   aod::Tracks,
                   aod::McTrackLabels,
                   aod::UpgradeTofs,
                   aod::UpgradeRichs,
                   aod::UpgradeTrkPids
                   >::iterator const& trackJoin) override;
  void postRun(o2::framework::PostRunContext& ctx) override;

private:
  float getNsigma(int pdgCode, const aod::UpgradeTofs& tofPid,
                  const aod::UpgradeRichs& richPid, const aod::UpgradeTrkPids& trackerPid,
                  int detectorType) const;

  int getPdgIndex(int pdgCode) const;

};

#endif // PID_EVALUATION_TASK_H