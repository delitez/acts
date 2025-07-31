// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootJetWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <algorithm>
#include <cmath>
#include <ios>
#include <limits>
#include <numbers>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <HepMC3/GenParticle.h>
#include <TFile.h>
#include <TTree.h>

const Acts::MagneticFieldContext magFieldContext;

namespace {
double calcSumPt2(const Acts::Vertex& vtx) {
  double sumPt2 = 0;
  const double minTrkWeight = 0.0;
  for (const auto& trk : vtx.tracks()) {
    if (trk.trackWeight > minTrkWeight) {
      double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                      ->transverseMomentum();
      sumPt2 += pt * pt;
    }
  }
  return sumPt2;
}
}  // namespace

Acts::MagneticFieldProvider::Cache magFieldCache() {
  return Acts::NullBField{}.makeCache(magFieldContext);
}

namespace ActsExamples {

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

RootJetWriter::RootJetWriter(const RootJetWriter::Config& config,
                             Acts::Logging::Level level)
    : WriterT(config.inputTrackJets, "RootJetWriter", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
  }
  if (m_cfg.inputTrackJets.empty()) {
    throw std::invalid_argument("Missing track jets input collection");
  }
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Missing vertices input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument(
        "Missing track to particle map input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Initialize the input handles
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackJets.initialize(m_cfg.inputTrackJets);
  m_inputVertices.initialize(m_cfg.inputVertices);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  m_inputParticles.initialize(m_cfg.inputParticles);

  // Tools
  Acts::EigenStepper<> stepper(m_cfg.field);
  m_propagator = std::make_shared<Propagator>(stepper);
  Acts::ImpactPointEstimator::Config ipEstCfg(m_cfg.field, m_propagator);
  m_ipEst = std::make_shared<Acts::ImpactPointEstimator>(ipEstCfg);

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  m_outputTree->Branch("EventNumber", &m_eventNr);

  m_outputTree->Branch("recovertex_x", &m_recovtx_x);
  m_outputTree->Branch("recovertex_y", &m_recovtx_y);
  m_outputTree->Branch("recovertex_z", &m_recovtx_z);
  m_outputTree->Branch("recovertex_t", &m_recovtx_t);
  m_outputTree->Branch("recovertex_sumPt2", &m_recovtx_sumPt2);
  m_outputTree->Branch("recovertex_isHS", &m_recovtx_isHS);
  m_outputTree->Branch("recovertex_isPU", &m_recovtx_isPU);
  m_outputTree->Branch("recovertex_isSec", &m_recovtx_isSec);
  m_outputTree->Branch("recovertex_eta", &m_recovtx_eta);
  m_outputTree->Branch("recovertex_theta", &m_recovtx_theta);
  m_outputTree->Branch("recovertex_phi", &m_recovtx_phi);

  m_outputTree->Branch("secvertex_x", &m_secvtx_x);
  m_outputTree->Branch("secvertex_y", &m_secvtx_y);
  m_outputTree->Branch("secvertex_z", &m_secvtx_z);
  m_outputTree->Branch("secvertex_t", &m_secvtx_t);
  m_outputTree->Branch("secvtx_Lxy", &m_secvtx_Lxy);
  m_outputTree->Branch("secvertex_eta", &m_secvtx_eta);
  m_outputTree->Branch("secvertex_theta", &m_secvtx_theta);
  m_outputTree->Branch("secvertex_phi", &m_secvtx_phi);
  m_outputTree->Branch("jet_secvtx_deltaR", &m_jet_secvtx_deltaR);
  m_outputTree->Branch("secvtx_pt", &m_secvtx_pt);

  // The jets
  m_outputTree->Branch("jet_pt", &m_jet_pt);
  m_outputTree->Branch("jet_eta", &m_jet_eta);
  m_outputTree->Branch("jet_phi", &m_jet_phi);
  m_outputTree->Branch("jet_ncomponents", &m_jet_ncomponents);
  m_outputTree->Branch("jet_components", &m_jet_components);
  // for each jet, the indices of the  tracks it contains
  m_outputTree->Branch("jet_tracks_idx", &m_jet_tracks_idx);
  m_outputTree->Branch("jet_ntracks", &m_jet_ntracks);
  m_outputTree->Branch("jet_isPU", &m_jet_isPU);
  m_outputTree->Branch("jet_isHS", &m_jet_isHS);
  m_outputTree->Branch("jet_label", &m_jet_label);
  m_outputTree->Branch("jet_label_hadron_pt", &m_jet_label_hadron_pt);
  m_outputTree->Branch("jet_num_sec_vtx", &m_jet_num_sec_vtx);

  // Tracks in jets
  // for each track (that is matched to a jet), the index of the vertex it
  // belongs to
  m_outputTree->Branch("track_matched_secvtx_idx", &m_matched_secvtx_idx);
  m_outputTree->Branch("track_matched_jet_idx",
                       &m_matched_jet_idx);  // for each track, the index of the
                                             // jet it belongs to
  m_outputTree->Branch("track_prob", &m_tracks_prob);
  m_outputTree->Branch("track_d0", &m_trk_d0);
  m_outputTree->Branch("track_z0", &m_trk_z0);
  m_outputTree->Branch("track_signedd0", &m_trk_signed_d0);
  m_outputTree->Branch("track_signedd0sig", &m_trk_signed_d0sig);
  m_outputTree->Branch("track_signedz0sinTheta", &m_trk_signed_z0sinTheta);
  m_outputTree->Branch("track_signedz0sinThetasig",
                       &m_trk_signed_z0sinThetasig);
  m_outputTree->Branch("track_eta", &m_trk_eta);
  m_outputTree->Branch("track_theta", &m_trk_theta);
  m_outputTree->Branch("track_phi", &m_trk_phi);
  m_outputTree->Branch("track_pt", &m_trk_pt);
  m_outputTree->Branch("track_qOverP", &m_trk_qOverP);
  m_outputTree->Branch("track_t", &m_trk_t);

  m_outputTree->Branch("track_var_d0", &m_trk_var_d0);
  m_outputTree->Branch("track_var_z0", &m_trk_var_z0);
  m_outputTree->Branch("track_var_phi", &m_trk_var_phi);
  m_outputTree->Branch("track_var_theta", &m_trk_var_theta);
  m_outputTree->Branch("track_var_qOverP", &m_trk_var_qOverP);

  m_outputTree->Branch("track_cov_d0z0", &m_trk_cov_d0z0);
  m_outputTree->Branch("track_cov_d0phi", &m_trk_cov_d0phi);
  m_outputTree->Branch("track_cov_d0theta", &m_trk_cov_d0theta);
  m_outputTree->Branch("track_cov_d0qOverP", &m_trk_cov_d0qOverP);
  m_outputTree->Branch("track_cov_z0phi", &m_trk_cov_z0phi);
  m_outputTree->Branch("track_cov_z0theta", &m_trk_cov_z0theta);
  m_outputTree->Branch("track_cov_z0qOverP", &m_trk_cov_z0qOverP);
  m_outputTree->Branch("track_cov_phitheta", &m_trk_cov_phitheta);
  m_outputTree->Branch("track_cov_phiqOverP", &m_trk_cov_phiqOverP);
  m_outputTree->Branch("track_cov_tehtaqOverP", &m_trk_cov_thetaqOverP);

  m_outputTree->Branch("jet_track_deltaR_all",
                       &m_jet_track_deltaR_all);  // deltaR for all tracks in
                                                  // the jet, not only the ones
                                                  // matched to a jet

  m_outputTree->Branch("jet_track_deltaR_matched",
                       &m_jet_track_deltaR_matched);  // deltaR for tracks that
                                                      // are matched to a jet
}

ProcessCode RootJetWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote states of trajectories to tree '"
            << m_cfg.treeName << "' in '" << m_cfg.treeName << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootJetWriter::writeT(const AlgorithmContext& ctx,
                                  const TrackJetContainer& trackJets) {
  // Exclusive access to all the writing procedure
  std::lock_guard<std::mutex> lock(m_writeMutex);

  clear();

  if (m_outputFile == nullptr) {
    return ProcessCode::ABORT;  // make it an error condition
  }

  // Read input collections

  const auto& tracks = m_inputTracks(ctx);
  const auto& vertices = m_inputVertices(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& particles = m_inputParticles(ctx);

  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputTracks << " "
                                        << tracks.size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputTrackJets << " "
                                        << trackJets.size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputVertices << " "
                                        << m_inputVertices(ctx).size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputTrackParticleMatching
                                        << " " << trackParticleMatching.size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputParticles << " "
                                        << particles.size());

  m_eventNr = ctx.eventNumber;

  // Read the vertices and find the HS vertex
  // Find the vertex with the highest sumPt2
  double maxSumPt2 = -1.;

  const Acts::Vertex* hsVtx = nullptr;

  for (const auto& vtx : vertices) {
    double sumPt2 = calcSumPt2(vtx);
    if (sumPt2 > maxSumPt2) {
      maxSumPt2 = sumPt2;
      hsVtx = &vtx;
    }
  }

  if (hsVtx == nullptr) {
    ACTS_ERROR("No hard-scatter vertex found");
    return ProcessCode::SUCCESS;
  }

  ACTS_VERBOSE("Hard-scatter vertex found at position: ("
               << hsVtx->position().x() << ", " << hsVtx->position().y() << ", "
               << hsVtx->position().z() << ")");

  Acts::ImpactPointEstimator::State state{magFieldCache()};
  std::vector<Acts::Vector4> secondaryVertices;

  std::map<std::size_t, std::vector<int>> secVerticesByJet;

  double secVtx_hs_Lxy = -9999;  // Lxy for each secondary vertex

  // Loop over the tracks
  for (std::size_t itrk = 0; itrk < tracks.size(); itrk++) {
    double signed_d0 = -9999;
    double signed_z0SinTheta = -9999;
    double signed_d0_err = 1.;
    double signed_z0SinTheta_err = 1.;

    int matched_jet_idx = -999;
    int isecvtx = -999;  // flag to check if we found a secondary vertex

    double covd0 = -999;
    double covz0 = -999;
    double covphi = -999;
    double covtheta = -999;
    double covqOverP = -999;
    double covd0z0 = -999;
    double covd0phi = -999;
    double covd0theta = -999;
    double covd0qOverP = -999;
    double covz0phi = -999;
    double covz0theta = -999;
    double covz0qOverP = -999;
    double covphitheta = -999;
    double covphiqOverP = -999;
    double covthetaqOverP = -999;
    double deltaR_all = -999;

    const auto trk = tracks.at(itrk);
    const auto params = trk.parameters();

    auto boundParams = trk.createParametersAtReference();

    for (std::size_t ijet = 0; ijet < trackJets.size(); ++ijet) {
      auto& jet = trackJets.at(ijet);
      Acts::Vector3 jetMom = jet.getFourMomentum().template head<3>();
      Acts::Vector3 trkMom = trk.momentum();
      // Before jet-track matching, calculate delta R for all tracks and jets
      deltaR_all = Acts::VectorHelpers::deltaR(jetMom, trkMom);

      // Check if this track belongs to a jet and compute the IPs
      auto jtrks = jet.getTracks();

      if (std::ranges::find(jtrks, itrk) == jtrks.end()) {
        continue;
      }

      ACTS_DEBUG("Track " << itrk << " is in jet " << ijet);
      matched_jet_idx = ijet;

      auto ipAndSigma =
          m_ipEst->getImpactParameters(boundParams, *hsVtx, gctx_, mctx_, true);
      auto vszs = m_ipEst->getLifetimeSignOfTrack(
          boundParams, *hsVtx, jet.getDirection(), gctx_, mctx_);

      if (!ipAndSigma.ok() || !vszs.ok()) {
        continue;
      }

      auto deltaR_matched = Acts::VectorHelpers::deltaR(jetMom, trkMom);
      m_jet_track_deltaR_matched.push_back(deltaR_matched);

      double absD0 = std::abs((*ipAndSigma).d0);
      double absZ0 = std::abs((*ipAndSigma).z0);

      signed_d0 = absD0 * (*vszs).first;
      double theta = boundParams.theta();
      signed_z0SinTheta = absZ0 * std::sin(theta) * (*vszs).second;
      signed_d0_err = (*ipAndSigma).sigmaD0;
      signed_z0SinTheta_err =
          (*ipAndSigma).sigmaZ0;  // not correct yet - have to take into
                                  // account the sin(theta) factor

      // Get the associated truth particle for this track
      // if matched, take this truth particle's initial ? position
      // make a vector for 2ndvtxpos (out of the loop) am I close enough to
      // the position of the truth particle?

      auto match = trackParticleMatching.find(itrk);

      if (match == trackParticleMatching.end() ||
          !match->second.particle.has_value()) {
        continue;
      }

      auto barcode = match->second.particle.value();

      auto findParticle = particles.find(barcode);
      if (findParticle == particles.end()) {
        continue;
      }

      auto initParticle = *findParticle;  // This is the initial particle
                                          // matched to the track with type

      Acts::Vector4 initParticlePosition = initParticle.fourPosition();
      Acts::Vector4 hsPosition = hsVtx->fullPosition();

      if (jet.getLabel() == JetLabel::BJet) {
        auto hs_it = std::ranges::find_if(
            secondaryVertices.begin(), secondaryVertices.end(),
            [&hsPosition](const Acts::Vector4& vtx) {
              return (vtx - hsPosition).norm() < 1e-3;  // tolerance
            });

        if (hs_it == secondaryVertices.end()) {
          auto it = std::ranges::find_if(
              secondaryVertices.begin(), secondaryVertices.end(),
              [&initParticlePosition](const Acts::Vector4& vtx) {
                return (vtx - initParticlePosition).norm() <
                       1e-3;  // tolerance // can also make it configurable
              });

          if (it == secondaryVertices.end()) {
            double initParticleR =
                std::sqrt(initParticlePosition[Acts::ePos0] *
                              initParticlePosition[Acts::ePos0] +
                          initParticlePosition[Acts::ePos1] *
                              initParticlePosition[Acts::ePos1]);
            if (m_cfg.useOnlyBeamPipe &&
                (initParticleR < 1.0 || initParticleR > 30.0)) {
              ACTS_DEBUG("Skipping secondary vertex matching for track "
                         << itrk << " in jet " << ijet
                         << " due to beam pipe radius condition");
              continue;
            }
            secondaryVertices.push_back(initParticlePosition);

            isecvtx = secondaryVertices.size() - 1;

          } else {
            isecvtx = std::distance(secondaryVertices.begin(), it);
          }
          secVerticesByJet[ijet].push_back(isecvtx);
        } else {
          ACTS_VERBOSE("Hard scatter vertex already exists at position: ("
                       << hsPosition[Acts::ePos0] << ", "
                       << hsPosition[Acts::ePos1] << ", "
                       << hsPosition[Acts::ePos2] << ")");
        }
      } else {
        ACTS_VERBOSE("Jet is not a B-jet, skipping secondary vertex matching");
      }
      secVtx_hs_Lxy =
          std::sqrt(std::pow(std::abs(initParticlePosition[Acts::ePos0] -
                                      hsPosition[Acts::ePos0]),
                             2) +
                    std::pow(std::abs(initParticlePosition[Acts::ePos1] -
                                      hsPosition[Acts::ePos1]),
                             2));

    }  // loop on jets

    m_jet_track_deltaR_all.push_back(deltaR_all);
    float trk_theta = params[Acts::eBoundTheta];
    float trk_eta = std::atanh(std::cos(trk_theta));
    float trk_qop = params[Acts::eBoundQOverP];
    float trk_p = std::abs(1.0 / trk_qop);
    float trk_pt = trk_p * std::sin(trk_theta);

    m_tracks_prob.push_back(1.);  // todo
    m_trk_d0.push_back(params[Acts::eBoundLoc0]);
    m_trk_z0.push_back(params[Acts::eBoundLoc1]);

    m_trk_signed_d0.push_back(signed_d0);
    m_trk_signed_d0sig.push_back(signed_d0 / signed_d0_err);
    m_trk_signed_z0sinTheta.push_back(signed_z0SinTheta);
    m_trk_signed_z0sinThetasig.push_back(signed_z0SinTheta /
                                         signed_z0SinTheta_err);

    // for each track, the index of the jet it belongs to
    m_matched_jet_idx.push_back(matched_jet_idx);
    // fill number of secondary vertices for each track (PROBABLY WRONG:
    // MATCHING MUST BE IN THE JET LOOP-that are matched to a jet) HINT: this
    // might also include the hard scatter vertex
    m_matched_secvtx_idx.push_back(isecvtx);

    m_trk_eta.push_back(trk_eta);
    m_trk_theta.push_back(trk_theta);
    m_trk_phi.push_back(params[Acts::eBoundPhi]);
    m_trk_pt.push_back(trk_pt);
    m_trk_qOverP.push_back(trk_qop);
    m_trk_t.push_back(params[Acts::eBoundTime]);

    auto cov = trk.covariance();

    covd0 = cov(Acts::eBoundLoc0, Acts::eBoundLoc0);
    covz0 = cov(Acts::eBoundLoc1, Acts::eBoundLoc1);
    covphi = cov(Acts::eBoundPhi, Acts::eBoundPhi);
    covtheta = cov(Acts::eBoundTheta, Acts::eBoundTheta);
    covqOverP = cov(Acts::eBoundQOverP, Acts::eBoundQOverP);
    covd0z0 = cov(Acts::eBoundLoc0, Acts::eBoundLoc1);
    covd0phi = cov(Acts::eBoundLoc0, Acts::eBoundPhi);
    covd0theta = cov(Acts::eBoundLoc0, Acts::eBoundTheta);
    covd0qOverP = cov(Acts::eBoundLoc0, Acts::eBoundQOverP);
    covz0phi = cov(Acts::eBoundLoc1, Acts::eBoundPhi);
    covz0theta = cov(Acts::eBoundLoc1, Acts::eBoundTheta);
    covz0qOverP = cov(Acts::eBoundLoc1, Acts::eBoundQOverP);
    covphitheta = cov(Acts::eBoundPhi, Acts::eBoundTheta);
    covphiqOverP = cov(Acts::eBoundPhi, Acts::eBoundQOverP);
    covthetaqOverP = cov(Acts::eBoundTheta, Acts::eBoundQOverP);

    m_trk_var_d0.push_back(covd0);
    m_trk_var_z0.push_back(covz0);
    m_trk_var_phi.push_back(covphi);
    m_trk_var_theta.push_back(covtheta);
    m_trk_var_qOverP.push_back(covqOverP);
    m_trk_cov_d0z0.push_back(covd0z0);
    m_trk_cov_d0phi.push_back(covd0phi);
    m_trk_cov_d0theta.push_back(covd0theta);
    m_trk_cov_d0qOverP.push_back(covd0qOverP);
    m_trk_cov_z0phi.push_back(covz0phi);
    m_trk_cov_z0theta.push_back(covz0theta);
    m_trk_cov_z0qOverP.push_back(covz0qOverP);
    m_trk_cov_phitheta.push_back(covphitheta);
    m_trk_cov_phiqOverP.push_back(covphiqOverP);
    m_trk_cov_thetaqOverP.push_back(covthetaqOverP);

  }  // loop over tracks

  for (const auto& secVtx : secondaryVertices) {
    double R = std::sqrt(secVtx[Acts::ePos0] * secVtx[Acts::ePos0] +
                         secVtx[Acts::ePos1] * secVtx[Acts::ePos1]);

    m_secvtx_x.push_back(secVtx[Acts::ePos0]);
    m_secvtx_y.push_back(secVtx[Acts::ePos1]);
    m_secvtx_z.push_back(secVtx[Acts::ePos2]);
    m_secvtx_t.push_back(secVtx[Acts::eTime]);

    double sec_vtx_theta = std::atan2(R, secVtx[Acts::ePos2]);
    double eta = -std::log(std::tan(sec_vtx_theta / 2.0));
    m_secvtx_eta.push_back(eta);
    m_secvtx_theta.push_back(sec_vtx_theta);
    m_secvtx_phi.push_back(phi(secVtx));
    m_secvtx_Lxy.push_back(secVtx_hs_Lxy);

    double deltaR_jet_secvtx = -9999;
    for (std::size_t j = 0; j < trackJets.size(); ++j) {
      // Calculate deltaR between jet and secondary vertex
      auto& jet = trackJets.at(j);
      Acts::Vector4 jet_4mom = jet.getFourMomentum();
      Acts::Vector3 jet_3mom{jet_4mom[0], jet_4mom[1], jet_4mom[2]};
      double jetTheta = theta(jet_3mom);
      double jetEta = std::atanh(std::cos(jetTheta));
      double jetPhi = phi(jet_4mom);

      deltaR_jet_secvtx =
          std::sqrt((jetEta - eta) * (jetEta - eta) +
                    (jetPhi - phi(secVtx)) * (jetPhi - phi(secVtx)));
    }
    m_jet_secvtx_deltaR.push_back(deltaR_jet_secvtx);
  }

  for (auto& [ijet, secVtxIndices] : secVerticesByJet) {
    std::ranges::sort(secVtxIndices);
    const auto ret = std::ranges::unique(secVtxIndices);
    secVtxIndices.erase(ret.begin(), ret.end());
    std::size_t nSecVtx = secVtxIndices.size();
    m_jet_num_sec_vtx.push_back(nSecVtx);
  }

  for (std::size_t ijets = 0; ijets < trackJets.size(); ++ijets) {
    Acts::Vector4 jet_4mom = trackJets[ijets].getFourMomentum();
    Acts::Vector3 jet_3mom{jet_4mom[0], jet_4mom[1], jet_4mom[2]};
    float jet_theta = theta(jet_3mom);
    m_jet_pt.push_back(perp(jet_4mom));
    m_jet_eta.push_back(std::atanh(std::cos(jet_theta)));
    m_jet_phi.push_back(phi(jet_4mom));

    const auto& jet = trackJets[ijets];

    m_jet_ncomponents.push_back(jet.getConstituents().size());
    m_jet_components.push_back(jet.getConstituents());

    auto jtrks = jet.getTracks();
    m_jet_tracks_idx.emplace_back(jtrks.begin(), jtrks.end());

    m_jet_ntracks.push_back(jet.getTracks().size());
    m_jet_label.push_back(static_cast<int>(jet.getLabel()));
    m_jet_isHS.push_back(1);  // this is not correct

    const auto* labelHadron = jet.getLabelHadron();

    if (labelHadron != nullptr) {
      m_jet_label_hadron_pt.push_back(labelHadron->momentum().pt());
    }

  }  // jets

  // Fill the reconstructed vertex information
  for (const auto& vtx : vertices) {
    const auto pos4 = vtx.fullPosition();
    const auto pos3 = vtx.position();
    m_recovtx_x.push_back(pos4[Acts::ePos0]);
    m_recovtx_y.push_back(pos4[Acts::ePos1]);
    m_recovtx_z.push_back(pos4[Acts::ePos2]);
    m_recovtx_t.push_back(pos4[Acts::eTime]);
    m_recovtx_eta.push_back(eta(pos3));
    m_recovtx_theta.push_back(theta(pos3));
    m_recovtx_phi.push_back(phi(pos3));
    m_recovtx_sumPt2.push_back(calcSumPt2(vtx));

    m_recovtx_isHS.push_back(static_cast<int>(&vtx == hsVtx));
  }  // vertices

  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

void ActsExamples::RootJetWriter::clear() {
  // Vertices
  m_recovtx_x.clear();
  m_recovtx_y.clear();
  m_recovtx_z.clear();
  m_recovtx_t.clear();
  m_recovtx_sumPt2.clear();
  m_recovtx_isHS.clear();
  m_recovtx_isPU.clear();
  m_recovtx_isSec.clear();
  m_recovtx_eta.clear();
  m_recovtx_theta.clear();
  m_recovtx_phi.clear();
  m_matched_secvtx_idx.clear();    // for each track (that is matched to a jet),
                                   // the index of the vertex it belongs to
  m_jet_track_deltaR_all.clear();  // deltaR for all tracks in the jet, not only
                                   // the ones matched to a jet
  m_jet_track_deltaR_matched
      .clear();  // deltaR for tracks that are matched to a jet

  // clear secondary vertex information
  m_secvtx_x.clear();
  m_secvtx_y.clear();
  m_secvtx_z.clear();
  m_secvtx_t.clear();
  m_secvtx_Lxy.clear();  // Lxy for each secondary vertex
  m_secvtx_eta.clear();
  m_secvtx_theta.clear();
  m_secvtx_phi.clear();
  m_jet_secvtx_deltaR.clear();  // deltaR for each jet and secondary vertex
  m_secvtx_pt.clear();          // pt for each secondary vertex

  // Jets
  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_ncomponents.clear();
  m_jet_components.clear();
  m_jet_tracks_idx.clear();
  m_matched_jet_idx
      .clear();  // for each track, the index of the jet it belongs to
  m_jet_ntracks.clear();
  m_jet_isPU.clear();
  m_jet_isHS.clear();
  m_jet_label.clear();
  m_jet_label_hadron_pt.clear();
  m_jet_num_sec_vtx.clear();

  // Tracks
  m_tracks_prob.clear();
  m_trk_d0.clear();
  m_trk_z0.clear();
  m_trk_signed_d0.clear();
  m_trk_signed_d0sig.clear();
  m_trk_signed_z0sinTheta.clear();
  m_trk_signed_z0sinThetasig.clear();
  m_trk_eta.clear();
  m_trk_theta.clear();
  m_trk_phi.clear();
  m_trk_pt.clear();
  m_trk_qOverP.clear();
  m_trk_t.clear();
  m_trk_var_d0.clear();
  m_trk_var_z0.clear();
  m_trk_var_phi.clear();
  m_trk_var_theta.clear();
  m_trk_var_qOverP.clear();
  m_trk_cov_d0z0.clear();
  m_trk_cov_d0phi.clear();
  m_trk_cov_d0theta.clear();
  m_trk_cov_d0qOverP.clear();
  m_trk_cov_z0phi.clear();
  m_trk_cov_z0theta.clear();
  m_trk_cov_z0qOverP.clear();
  m_trk_cov_phitheta.clear();
  m_trk_cov_phiqOverP.clear();
  m_trk_cov_thetaqOverP.clear();
}

double ActsExamples::RootJetWriter::deltaR(
    TrackJet jet,
    Acts::TrackProxy<Acts::ConstVectorTrackContainer,
                     Acts::ConstVectorMultiTrajectory, std::shared_ptr, true>
        trk) {
  double jetPx = jet.getFourMomentum().x();
  double jetPy = jet.getFourMomentum().y();
  double jetPz = jet.getFourMomentum().z();
  double trackPx = trk.momentum().x();
  double trackPy = trk.momentum().y();
  double trackPz = trk.momentum().z();
  double pjet = std::sqrt(jetPx * jetPx + jetPy * jetPy + jetPz * jetPz);
  double ptrack =
      std::sqrt(trackPx * trackPx + trackPy * trackPy + trackPz * trackPz);

  // Calculate eta and phi for the jet and track
  // Note: eta = arctanh(pz/|p|), phi = atan2(py, px)
  // where theta is the polar angle of the momentum vector
  // and phi is the azimuthal angle in the xy-plane.
  // Here we use the four-momentum to calculate eta and phi.

  double jetEta = std::atanh(jetPz / pjet);
  double jetPhi = std::atan2(jetPy, jetPx);
  double trackEta = std::atanh(trackPz / ptrack);
  double trackPhi = std::atan2(trackPy, trackPx);

  // Calculate delta R
  double deltaEta = jetEta - trackEta;
  double deltaPhi = jetPhi - trackPhi;
  double deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
  return deltaR;
}
