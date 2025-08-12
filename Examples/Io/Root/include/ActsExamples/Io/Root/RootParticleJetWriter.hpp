// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
// This I need to compute the ip2d and ip3d - to do: move this into an algorithm
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include <array>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include <HepMC3/GenParticle.h>

class TFile;
class TTree;

namespace ActsExamples {

using TrackJetWriter = WriterT<TrackJetContainer>;
using VertexContainer = std::vector<Acts::Vertex>;

using Propagator = Acts::Propagator<Acts::EigenStepper<>>;

class RootParticleJetWriter final : public TrackJetWriter {
 public:
  struct Config {
    /// Input tracks
    std::string inputTracks;
    /// Input track jets
    std::string inputTrackJets;
    /// Input vertices
    std::string inputVertices;
    /// Input measurement particles map
    std::string inputTrackParticleMatching;
    /// Input particles
    std::string inputParticles;
    /// output filename.
    std::string filePath = "events.root";
    /// name of the output tree.
    std::string treeName = "events";
    /// file access mode.
    std::string fileMode = "RECREATE";
    /// Magnetic field provider
    std::shared_ptr<Acts::MagneticFieldProvider> field;
    /// Use only beam pipe radius for secondary vertex matching
    bool useOnlyBeamPipe = false;
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootParticleJetWriter(const Config& config, Acts::Logging::Level level);

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] tracks are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackJetContainer& trackJets) override;

 private:
  /// The config class
  Config m_cfg;

  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number of
  void clear();

  // Handles

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "inputTracks"};
  ReadDataHandle<TrackJetContainer> m_inputTrackJets{this, "inputTrackJets"};
  ReadDataHandle<VertexContainer> m_inputVertices{this, "inputVertices"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "inputTrackParticleMatching"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "inputParticles"};

  // Vertices//
  std::vector<float> m_recovtx_x;
  std::vector<float> m_recovtx_y;
  std::vector<float> m_recovtx_z;
  std::vector<float> m_recovtx_t;
  std::vector<float> m_recovtx_sumPt2;
  std::vector<int> m_recovtx_isHS;
  std::vector<int> m_recovtx_isPU;
  std::vector<int> m_recovtx_isSec;
  std::vector<int>
      m_matched_secvtx_idx;  // for each track (that is matched to a jet), the
                             // index of the vertex it belongs to

  std::vector<float> m_recovtx_eta;
  std::vector<float> m_recovtx_theta;
  std::vector<float> m_recovtx_phi;

  std::vector<float> m_secvtx_x;
  std::vector<float> m_secvtx_y;
  std::vector<float> m_secvtx_z;
  std::vector<float> m_secvtx_t;
  std::vector<float> m_secvtx_hs_Lxy;  // Lxy for each secondary vertex
  std::vector<float> m_secvtx_eta;
  std::vector<float> m_secvtx_theta;
  std::vector<float> m_secvtx_phi;
  std::vector<float> m_jet_secvtx_deltaR;  // deltaR for each jet and secondary
                                           // vertex
  std::vector<float> m_secvtx_pt;          // pt for each secondary vertex
  // Jets//

  std::vector<int>
      m_lightjet_nsecvtx;  // number of secondary vertices for light jets
  std::vector<int> m_bjet_nsecvtx;  // number of secondary vertices for b-jets
  std::vector<int> m_cjet_nsecvtx;  // number of secondary vertices for c-jets
  std::vector<int>
      m_ntracks_per_secvtx;  // number of tracks for each secondary vertex

  // skipping jet_m, jet_q

  std::vector<float> m_jet_pt, m_jet_eta, m_jet_phi;
  std::vector<int> m_jet_ncomponents;
  std::vector<std::vector<int>> m_jet_components;
  std::vector<std::vector<int>> m_jet_tracks_idx;
  std::vector<int> m_jet_ntracks;
  std::vector<int> m_jet_label;
  std::vector<int> m_jet_isPU;
  std::vector<int> m_jet_isHS;
  std::vector<float> m_jet_label_hadron_pt;
  std::vector<int> m_jet_num_sec_vtx;
  std::vector<float> m_jet_secvtx_dPhi;

  // Tracks in jets
  std::vector<float> m_tracks_prob;
  std::vector<float> m_trk_d0;
  std::vector<float> m_trk_z0;
  std::vector<float> m_trk_signed_d0;
  std::vector<float> m_trk_signed_z0sinTheta;
  std::vector<float> m_trk_signed_d0sig;
  std::vector<float> m_trk_signed_z0sinThetasig;
  std::vector<float> m_trk_eta;
  std::vector<float> m_trk_theta;
  std::vector<float> m_trk_phi;
  std::vector<float> m_trk_pt;
  std::vector<float> m_trk_qOverP;
  std::vector<float> m_trk_t;

  std::vector<float> m_trk_var_d0;
  std::vector<float> m_trk_var_z0;
  std::vector<float> m_trk_var_phi;
  std::vector<float> m_trk_var_theta;
  std::vector<float> m_trk_var_qOverP;
  std::vector<float> m_trk_cov_d0z0;
  std::vector<float> m_trk_cov_d0phi;
  std::vector<float> m_trk_cov_d0theta;
  std::vector<float> m_trk_cov_d0qOverP;
  std::vector<float> m_trk_cov_z0phi;
  std::vector<float> m_trk_cov_z0theta;
  std::vector<float> m_trk_cov_z0qOverP;
  std::vector<float> m_trk_cov_phitheta;
  std::vector<float> m_trk_cov_phiqOverP;
  std::vector<float> m_trk_cov_thetaqOverP;

  // the index of the jet the tracks belong to
  std::vector<int> m_matched_jet_idx;

  // deltaR for all tracks in the jet, not only the ones matched to a jet
  std::vector<float> m_jet_track_deltaR_all;
  // deltaR for tracks that are matched to a jet
  std::vector<float> m_jet_track_deltaR_matched;

  // Tools

  std::shared_ptr<Propagator> m_propagator;
  std::shared_ptr<Acts::ImpactPointEstimator> m_ipEst;

  Acts::GeometryContext gctx_;
  Acts::MagneticFieldContext mctx_;
  double deltaR(
      TrackJet jet,
      Acts::TrackProxy<Acts::ConstVectorTrackContainer,
                       Acts::ConstVectorMultiTrajectory, std::shared_ptr, true>
          trk);

  double correctedPhi(double phi1, double phi2);
  double correctedSinglePhi(double phi1);
};

}  // namespace ActsExamples
