#include "ALICE3/DataModel/PIDEvaluation.h"

namespace o2 {
namespace aod {
  struct UpgradeRichs {
    float getNsigmaElectron() const { return 999.f; }
    float getNsigmaMuon() const { return 999.f; }
    float getNsigmaPion() const { return 999.f; }
    float getNsigmaKaon() const { return 999.f; }
    float getNsigmaProton() const { return 999.f; }
    float getNsigmaDeuteron() const { return 999.f; }
    float getNsigmaTriton() const { return 999.f; }  
    float getNsigmaHelium3() const { return 999.f; } 
    float getNsigmaAlpha() const { return 999.f; }   
  };
}
}

enum DetectorType {
    kTOF,
    kRICH,
    kTracker
};

void PidEvaluation::init(o2::framework::InitContext& ctx) {
  for (int i = 0; i < kNumHypothesisParticles; ++i) {
    mPdgToIndexMap[mHypothesisPdgCodes[i]] = i;
  }

  mLogBins.clear();
  mNumBins = mNumLogBins.value;
  double logMin = std::log10(mPMin.value);
  double logMax = std::log10(mPMax.value);
  double dLog = (logMax - logMin) / mNumBins;
  for (int i = 0; i <= mNumBins; ++i) {
    mLogBins.push_back(std::pow(10, logMin + i * dLog));
  }

  mPurityNumerator.resize(2, std::vector<std::vector<double>>(kNumHypothesisParticles, std::vector<double>(mNumBins, 0.0)));
  mPurityDenominator.resize(2, std::vector<std::vector<double>>(kNumHypothesisParticles, std::vector<double>(mNumBins, 0.0)));
  mEfficiencyDenominator.resize(kNumHypothesisParticles, std::vector<double>(mNumBins, 0.0));

  const AxisSpec axisMomentum{mLogBins, "#it{p} (GeV/#it{c})"};

  for (int i = 0; i < kNumHypothesisParticles; ++i) {
    std::string particleName = mParticleNames[i];
    std::string particleNamePretty = mParticleNamesPretty[i];

    histos.add(Form("hPurity_%s_Baseline", particleName.c_str()),
               Form("Purity of identified %s (Baseline); #it{p} (GeV/#it{c}); Purity", particleNamePretty.c_str()),
               kTH1F, {axisMomentum});

    histos.add(Form("hPurity_%s_WithTracker", particleName.c_str()),
               Form("Purity of identified %s (Baseline + Tracker); #it{p} (GeV/#it{c}); Purity", particleNamePretty.c_str()),
               kTH1F, {axisMomentum});

    histos.add(Form("hEfficiency_%s_Baseline", particleName.c_str()),
               Form("Efficiency of %s (Baseline); #it{p} (GeV/#it{c}); Efficiency", particleNamePretty.c_str()),
               kTH1F, {axisMomentum});

    histos.add(Form("hEfficiency_%s_WithTracker", particleName.c_str()),
               Form("Efficiency of %s (Baseline + Tracker); #it{p} (GeV/#it{c}); Efficiency", particleNamePretty.c_str()),
               kTH1F, {axisMomentum});

    histos.add(Form("hNSigma_TOF_%s", particleName.c_str()),
               Form("N#sigma_{TOF} for true %s; #it{p} (GeV/#it{c}); N#sigma_{TOF}", particleNamePretty.c_str()),
               kTH2F, {axisMomentum, AxisSpec{200, -10., 10., "N#sigma"}});
    histos.add(Form("hNSigma_RICH_%s", particleName.c_str()),
               Form("N#sigma_{RICH} for true %s; #it{p} (GeV/#it{c}); N#sigma_{RICH}", particleNamePretty.c_str()),
               kTH2F, {axisMomentum, AxisSpec{200, -10., 10., "N#sigma"}});
  }

  // Specific 2D plot for Alpha-Deuteron separation using TOF vs Tracker
  int alphaIdx = getPdgIndex(1000020040); // Alpha
  int deuteronIdx = getPdgIndex(1000010020); // Deuteron
  if (alphaIdx != -1 && deuteronIdx != -1) {
      histos.add(Form("hNSigma_TOF_vs_Tracker_True%s_Hyp%s", mParticleNames[alphaIdx].c_str(), mParticleNames[deuteronIdx].c_str()),
                 Form("N#sigma_{TOF} vs N#sigma_{Tracker} for true %s (Hypothesis %s); N#sigma_{TOF}; N#sigma_{Tracker}",
                      mParticleNamesPretty[alphaIdx].c_str(), mParticleNamesPretty[deuteronIdx].c_str()),
                 kTH2F, {AxisSpec{100, -5., 5., "N#sigma_{TOF}"}, AxisSpec{100, -5., 5., "N#sigma_{Tracker}"}});
      histos.add(Form("hNSigma_TOF_vs_Tracker_True%s_Hyp%s", mParticleNames[deuteronIdx].c_str(), mParticleNames[alphaIdx].c_str()),
                 Form("N#sigma_{TOF} vs N#sigma_{Tracker} for true %s (Hypothesis %s); N#sigma_{TOF}; N#sigma_{Tracker}",
                      mParticleNamesPretty[deuteronIdx].c_str(), mParticleNamesPretty[alphaIdx].c_str()),
                 kTH2F, {AxisSpec{100, -5., 5., "N#sigma_{TOF}"}, AxisSpec{100, -5., 5., "N#sigma_{Tracker}"}});
  }
}

int PidEvaluation::getPdgIndex(int pdgCode) const {
    auto it = mPdgToIndexMap.find(std::abs(pdgCode));
    if (it != mPdgToIndexMap.end()) {
        return it->second;
    }
    return -1;
}

float PidEvaluation::getNsigma(int pdgCode, const aod::UpgradeTofs& tofPid,
                                   const aod::UpgradeRichs& richPid, const aod::UpgradeTrkPids& trackerPid,
                                   int detectorType) const {
    float nsigma = 999.f; // Default "bad" value

    switch (std::abs(pdgCode)) {
        case 11: // Electron
            if (detectorType == kTOF) nsigma = tofPid.getNsigmaElectron();
            else if (detectorType == kRICH) nsigma = richPid.getNsigmaElectron();
            else if (detectorType == kTracker) nsigma = trackerPid.getNsigmaElectron();
            break;
        case 13: // Muon
            if (detectorType == kTOF) nsigma = tofPid.getNsigmaMuon();
            else if (detectorType == kRICH) nsigma = richPid.getNsigmaMuon();
            else if (detectorType == kTracker) nsigma = trackerPid.getNsigmaMuon();
            break;
        case 211: // Pion
            if (detectorType == kTOF) nsigma = tofPid.getNsigmaPion();
            else if (detectorType == kRICH) nsigma = richPid.getNsigmaPion();
            else if (detectorType == kTracker) nsigma = trackerPid.getNsigmaPion();
            break;
        case 321: // Kaon
            if (detectorType == kTOF) nsigma = tofPid.getNsigmaKaon();
            else if (detectorType == kRICH) nsigma = richPid.getNsigmaKaon();
            else if (detectorType == kTracker) nsigma = trackerPid.getNsigmaKaon();
            break;
        case 2212: // Proton
            if (detectorType == kTOF) nsigma = tofPid.getNsigmaProton();
            else if (detectorType == kRICH) nsigma = richPid.getNsigmaProton();
            else if (detectorType == kTracker) nsigma = trackerPid.getNsigmaProton();
            break;
        case 1000010020: // Deuteron
            if (detectorType == kTracker) nsigma = trackerPid.getNsigmaDeuteron();
            // TOF/RICH may not have getNsigmaDeuteron() - no else if
            break;
        case 1000010030: // Triton
            if (detectorType == kTracker) nsigma = trackerPid.getNsigmaTriton();
            break;
        case 1000020030: // Helium-3
            if (detectorType == kTracker) nsigma = trackerPid.getNsigmaHelium3();
            break;
        case 1000020040: // Alpha
            if (detectorType == kTracker) nsigma = trackerPid.getNsigmaAlpha();
            break;
        default:
            break;
    }
    return nsigma;
}

void PidEvaluation::process(soa::Join<
                                   aod::Tracks,
                                   aod::McTrackLabels,
                                   aod::UpgradeTofs,
                                   aod::UpgradeRichs,
                                   aod::UpgradeTrkPids
                                   >::iterator const& trackJoin) {
  if (!trackJoin.has_mcParticle() || !trackJoin.has_UpgradeTofs() ||
      !trackJoin.has_UpgradeRichs() || !trackJoin.has_upgradeTrkPids()) {
    return;
  }

  const auto& track = trackJoin.get<aod::Tracks>();
  const auto& mcParticle = trackJoin.get<aod::McTrackLabels>().mcParticle();
  const auto& tofPid = trackJoin.get<aod::UpgradeTofs>();
  const auto& richPid = trackJoin.get<aod::UpgradeRichs>();
  const auto& trackerPid = trackJoin.get<aod::UpgradeTrkPids>();

  float p = track.p();
  int truePdg = std::abs(mcParticle.pdgCode());
  int truePdgIdx = getPdgIndex(truePdg);

  if (truePdgIdx == -1) {
      return;
  }

  auto it_p = std::upper_bound(mLogBins.begin(), mLogBins.end(), p);
  int pBin = std::distance(mLogBins.begin(), it_p) - 1;
  pBin = std::max(0, std::min(pBin, mNumBins - 1));

  mEfficiencyDenominator[truePdgIdx][pBin]++;

  histos.fill(HIST(Form("hNSigma_TOF_%s", mParticleNames[truePdgIdx].c_str())), p, getNsigma(truePdg, tofPid, richPid, trackerPid, kTOF));
  histos.fill(HIST(Form("hNSigma_RICH_%s", mParticleNames[truePdgIdx].c_str())), p, getNsigma(truePdg, tofPid, richPid, trackerPid, kRICH));

  int alphaIdx = getPdgIndex(1000020040);
  int deuteronIdx = getPdgIndex(1000010020);

  if (truePdg == 1000020040 && deuteronIdx != -1) { // True Alpha, comparing with Deuteron hypothesis
      histos.fill(HIST(Form("hNSigma_TOF_vs_Tracker_True%s_Hyp%s", mParticleNames[alphaIdx].c_str(), mParticleNames[deuteronIdx].c_str())),
                  getNsigma(1000020040, tofPid, richPid, trackerPid, kTOF),
                  getNsigma(1000020040, tofPid, richPid, trackerPid, kTracker));
  }
  if (truePdg == 1000010020 && alphaIdx != -1) { // True Deuteron, comparing with Alpha hypothesis
      histos.fill(HIST(Form("hNSigma_TOF_vs_Tracker_True%s_Hyp%s", mParticleNames[deuteronIdx].c_str(), mParticleNames[alphaIdx].c_str())),
                  getNsigma(1000010020, tofPid, richPid, trackerPid, kTOF),
                  getNsigma(1000010020, tofPid, richPid, trackerPid, kTracker));
  }

  for (int hypPdgIdx = 0; hypPdgIdx < kNumHypothesisParticles; ++hypPdgIdx) {
    int hypPdgCode = mHypothesisPdgCodes[hypPdgIdx];

    bool passBaselinePid = true;
    if (mDoTOF.value) {
        if (std::abs(getNsigma(hypPdgCode, tofPid, richPid, trackerPid, kTOF)) > mNsigmaCut.value) {
            passBaselinePid = false;
        }
    }
    if (mDoRICH.value) {
        if (std::abs(getNsigma(hypPdgCode, tofPid, richPid, trackerPid, kRICH)) > mNsigmaCut.value) {
            passBaselinePid = false;
        }
    }

    bool passWithTrackerPid = passBaselinePid;
    if (mDoTracker.value && passWithTrackerPid) {
        if (std::abs(getNsigma(hypPdgCode, tofPid, richPid, trackerPid, kTracker)) > mNsigmaCut.value) {
            passWithTrackerPid = false;
        }
    }

    if (passBaselinePid) {
      mPurityDenominator[0][hypPdgIdx][pBin]++;
      if (truePdgIdx == hypPdgIdx) {
        mPurityNumerator[0][hypPdgIdx][pBin]++;
      }
    }

    if (passWithTrackerPid) {
      mPurityDenominator[1][hypPdgIdx][pBin]++;
      if (truePdgIdx == hypPdgIdx) {
        mPurityNumerator[1][hypPdgIdx][pBin]++;
      }
    }
  }
}

void PidEvaluation::postRun(o2::framework::PostRunContext& ctx) {
  for (int i = 0; i < kNumHypothesisParticles; ++i) {
    std::string particleName = mParticleNames[i];

    TH1F* hPurityBaseline = histos.get<TH1F>(Form("hPurity_%s_Baseline", particleName.c_str()));
    TH1F* hEfficiencyBaseline = histos.get<TH1F>(Form("hEfficiency_%s_Baseline", particleName.c_str()));

    for (int bin = 0; bin < mNumBins; ++bin) {
      double purityNum = mPurityNumerator[0][i][bin];
      double purityDen = mPurityDenominator[0][i][bin];
      if (purityDen > 0) {
        hPurityBaseline->SetBinContent(bin + 1, purityNum / purityDen);
        hPurityBaseline->SetBinError(bin + 1, std::sqrt(purityNum / purityDen * (1 - purityNum / purityDen) / purityDen));
      } else {
        hPurityBaseline->SetBinContent(bin + 1, 0.0);
        hPurityBaseline->SetBinError(bin + 1, 0.0);
      }

      double effNum = mPurityNumerator[0][i][bin];
      double effDen = mEfficiencyDenominator[i][bin];
      if (effDen > 0) {
        hEfficiencyBaseline->SetBinContent(bin + 1, effNum / effDen);
        hEfficiencyBaseline->SetBinError(bin + 1, std::sqrt(effNum / effDen * (1 - effNum / effDen) / effDen));
      } else {
        hEfficiencyBaseline->SetBinContent(bin + 1, 0.0);
        hEfficiencyBaseline->SetBinError(bin + 1, 0.0);
      }
    }

    TH1F* hPurityWithTracker = histos.get<TH1F>(Form("hPurity_%s_WithTracker", particleName.c_str()));
    TH1F* hEfficiencyWithTracker = histos.get<TH1F>(Form("hEfficiency_%s_WithTracker", particleName.c_str()));

    for (int bin = 0; bin < mNumBins; ++bin) {
      double purityNum = mPurityNumerator[1][i][bin];
      double purityDen = mPurityDenominator[1][i][bin];
      if (purityDen > 0) {
        hPurityWithTracker->SetBinContent(bin + 1, purityNum / purityDen);
        hPurityWithTracker->SetBinError(bin + 1, std::sqrt(purityNum / purityDen * (1 - purityNum / purityDen) / purityDen));
      } else {
        hPurityWithTracker->SetBinContent(bin + 1, 0.0);
        hPurityWithTracker->SetBinError(bin + 1, 0.0);
      }

      double effNum = mPurityNumerator[1][i][bin];
      double effDen = mEfficiencyDenominator[i][bin];
      if (effDen > 0) {
        hEfficiencyWithTracker->SetBinContent(bin + 1, effNum / effDen);
        hEfficiencyWithTracker->SetBinError(bin + 1, std::sqrt(effNum / effDen * (1 - effNum / effDen) / effDen));
      } else {
        hEfficiencyWithTracker->SetBinContent(bin + 1, 0.0);
        hEfficiencyWithTracker->SetBinError(bin + 1, 0.0);
      }
    }
  }
}

// Dummy PID tasks for local compilation if not using real O2 tasks
// REMOVE these in a full O2 environment where actual TOF/RICH PID tasks exist
struct OnTheFlyTofPidDummy : public o2::framework::AnalysisTask {
    o2::framework::Produces<o2::aod::UpgradeTofs> mOutput;
    void init(o2::framework::InitContext&) override {}
    void process(soa::Join<aod::Tracks>::iterator const& track) override {
        o2::aod::UpgradeTofs pid;
        pid.setNsigmaElectron(999.f); pid.setNsigmaMuon(999.f); pid.setNsigmaPion(999.f);
        pid.setNsigmaKaon(999.f); pid.setNsigmaProton(999.f); pid.setNsigmaDeuteron(999.f);
        pid.setNsigmaTriton(999.f); pid.setNsigmaHelium3(999.f); pid.setNsigmaAlpha(999.f);
        mOutput(pid);
    }
};

struct OnTheFlyRichPidDummy : public o2::framework::AnalysisTask {
    o2::framework::Produces<o2::aod::UpgradeRichs> mOutput;
    void init(o2::framework::InitContext&) override {}
    void process(soa::Join<aod::Tracks>::iterator const& track) override {
        o2::aod::UpgradeRichs pid;
        pid.setNsigmaElectron(999.f); pid.setNsigmaMuon(999.f); pid.setNsigmaPion(999.f);
        pid.setNsigmaKaon(999.f); pid.setNsigmaProton(999.f); pid.setNsigmaDeuteron(999.f);
        pid.setNsigmaTriton(999.f); pid.setNsigmaHelium3(999.f); pid.setNsigmaAlpha(999.f);
        mOutput(pid);
    }
};

// Assuming OnTheFlyTrackerPid definition is available from its own header
// #include "ALICE3/TableProducer/OTF/onTheFlyTrackerPid.h"

// This defineDataProcessing should be in your main analysis executable,
// not necessarily in PidEvaluation.cxx itself if it's compiled as a library.
// For testing, it can be here.
// You need to replace OnTheFlyTofPidDummy and OnTheFlyRichPidDummy
// with the actual O2 PID tasks in a real O2 environment.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
    return WorkflowSpec{
        adaptAnalysisTask<OnTheFlyTofPidDummy>(cfgc, TaskName{"OnTheFlyTofPid"}),
        adaptAnalysisTask<OnTheFlyRichPidDummy>(cfgc, TaskName{"OnTheFlyRichPid"}),
        adaptAnalysisTask<OnTheFlyTrackerPid>(cfgc, TaskName{"OnTheFlyTrackerPid"}), // Your Tracker PID task
        adaptAnalysisTask<PidEvaluation>(cfgc, TaskName{"PidEvaluation"}),
    };
}