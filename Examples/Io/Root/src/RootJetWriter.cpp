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
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsFatras/EventData/Barcode.hpp"


#include <cmath>
#include <ios>
#include <limits>
#include <numbers>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>


const Acts::MagneticFieldContext magFieldContext;

Acts::MagneticFieldProvider::Cache magFieldCache() {
  return Acts::NullBField{}.makeCache(magFieldContext );
}

namespace ActsExamples {

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

RootJetWriter::RootJetWriter(
    const RootJetWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputJets, "RootJetWriter", level),
      m_cfg(config) {

    if (m_cfg.inputJets.empty()) {
    throw std::invalid_argument("Missing jets input collection");
    }
    // if (m_cfg.inputProtoTracks.empty()) {
    // throw std::invalid_argument("Missing proto tracks input collection");
    // }
    // if (m_cfg.inputParticles.empty()) {
    // throw std::invalid_argument("Missing particles input collection");
    // }
    // if (m_cfg.inputSimHits.empty()) {
    // throw std::invalid_argument("Missing simulated hits input collection");
    // }
    if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
    }
    if (m_cfg.inputJets.empty()) {
    throw std::invalid_argument("Missing jets input collection");
    }
    if (m_cfg.inputTrackJets.empty()) {
    throw std::invalid_argument("Missing track jets input collection");
    }
    if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
    }
    if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
    }

    m_inputJets.initialize(m_cfg.inputJets);
    m_inputTracks.initialize(m_cfg.inputTracks);
    m_inputTrackJets.initialize(m_cfg.inputTrackJets);
    

  Acts::EigenStepper<> stepper(m_cfg.field);
  m_propagator = std::make_shared<Propagator>(stepper);

      //Setting up the ImpactPointEstimator
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


    // The Truth HS Vertex
    
    // The reco vertices
    
    // m_outputTree->Branch("recovertex_x",&m_vtx_x);
    // m_outputTree->Branch("recovertex_y",&m_vtx_y);
    // m_outputTree->Branch("recovertex_z",&m_vtx_z);
    // m_outputTree->Branch("recovertex_t",&m_vtx_t);
    // m_outputTree->Branch("recovertex_sumPt2",&m_vtx_sumPt2);
    // m_outputTree->Branch("recovertex_isHS",&m_vtx_isHS);
    // m_outputTree->Branch("recovertex_isPU",&m_vtx_isPU);

    /*
    m_outputTree->Branch("recovertex_cov_xx",&m_vtx_cov_xx);
    m_outputTree->Branch("recovertex_cov_xy",&m_vtx_cov_xy);
    m_outputTree->Branch("recovertex_cov_xz",&m_vtx_cov_xz);
    m_outputTree->Branch("recovertex_cov_xt",&m_vtx_cov_xt);
    m_outputTree->Branch("recovertex_cov_yy",&m_vtx_cov_yy);
    m_outputTree->Branch("recovertex_cov_yz",&m_vtx_cov_yz);
    m_outputTree->Branch("recovertex_cov_yt",&m_vtx_cov_yt);
    m_outputTree->Branch("recovertex_cov_zz",&m_vtx_cov_zz);
    m_outputTree->Branch("recovertex_cov_zt",&m_vtx_cov_zt);
    m_outputTree->Branch("recovertex_cov_tt",&m_vtx_cov_tt);
    */
    
    
    //m_outputTree->Branch("recovertex_tracks_idx",m_vtx_tracks_idx);
    //m_outputTree->Branch("recovertex_tracks_weight",m_vtx_tracks_weight);
    
    
    
    // The jets
    m_outputTree->Branch("jet_pt",&m_jet_pt);
    m_outputTree->Branch("jet_eta",&m_jet_eta);
    m_outputTree->Branch("jet_phi",&m_jet_phi);
    m_outputTree->Branch("jet_ncomponents",&m_jet_ncomponents);
    m_outputTree->Branch("jet_components",&m_jet_components);
    m_outputTree->Branch("jet_tracks_idx",&m_jet_tracks_idx);
    m_outputTree->Branch("jet_isPU",&m_jet_isPU);
    m_outputTree->Branch("jet_isHS",&m_jet_isHS);
    m_outputTree->Branch("jet_label",&m_jet_label);
    
    //Tracks in jets
    m_outputTree->Branch("track_prob",       &m_tracks_prob);
    m_outputTree->Branch("track_d0",         &m_trk_d0);
    m_outputTree->Branch("track_z0",         &m_trk_z0);
     m_outputTree->Branch("track_signedd0",             &m_trk_signed_d0);
    // m_outputTree->Branch("track_signedd0sig",          &m_trk_signed_d0sig);
     m_outputTree->Branch("track_signedz0sinTheta",     &m_trk_signed_z0sinTheta);
    //m_outputTree->Branch("track_signedz0sinThetasig",  &m_trk_signed_z0sinThetasig);
    m_outputTree->Branch("track_eta",        &m_trk_eta);
    m_outputTree->Branch("track_theta",      &m_trk_theta);
    m_outputTree->Branch("track_phi",        &m_trk_phi);
    m_outputTree->Branch("track_pt",         &m_trk_pt);
    m_outputTree->Branch("track_qOverP",     &m_trk_qOverP);
    m_outputTree->Branch("track_t",          &m_trk_t);
    m_outputTree->Branch("track_t30",        &m_trk_t30);
    m_outputTree->Branch("track_t60",        &m_trk_t60);
    m_outputTree->Branch("track_t90",        &m_trk_t90);
    m_outputTree->Branch("track_t120",       &m_trk_t120);
    m_outputTree->Branch("track_t180",       &m_trk_t180);
    m_outputTree->Branch("track_z",          &m_trk_z);

    // m_outputTree->Branch("track_var_d0",     &m_trk_var_d0);
    // m_outputTree->Branch("track_var_z0",     &m_trk_var_z0);
    // m_outputTree->Branch("track_var_phi",    &m_trk_var_phi);
    // m_outputTree->Branch("track_var_theta",  &m_trk_var_theta);
    // m_outputTree->Branch("track_var_qOverP", &m_trk_var_qOverP);

    // m_outputTree->Branch("track_cov_d0z0"       ,&m_trk_cov_d0z0);
    // m_outputTree->Branch("track_cov_d0phi"      ,&m_trk_cov_d0phi);
    // m_outputTree->Branch("track_cov_d0theta"    ,&m_trk_cov_d0theta);
    // m_outputTree->Branch("track_cov_d0qOverP"   ,&m_trk_cov_d0qOverP);
    // m_outputTree->Branch("track_cov_z0phi"      ,&m_trk_cov_z0phi);
    // m_outputTree->Branch("track_cov_z0theta"    ,&m_trk_cov_z0theta);
    // m_outputTree->Branch("track_cov_z0qOverP"   ,&m_trk_cov_z0qOverP);
    // m_outputTree->Branch("track_cov_phitheta"   ,&m_trk_cov_phitheta);
    // m_outputTree->Branch("track_cov_phiqOverP"  ,&m_trk_cov_phiqOverP);
    // m_outputTree->Branch("track_cov_tehtaqOverP",&m_trk_cov_thetaqOverP);

    // // Number of Innermost Pixel Layer Hits
    // m_outputTree->Branch("track_numPix1L" ,&m_trk_numPix1L);
    // // Number of Next to Innermost Pixel Layer Hits
    // m_outputTree->Branch("track_numPix2L" ,&m_trk_numPix2L);   
    // m_outputTree->Branch("track_numPix"   ,&m_trk_numPix);   
    // m_outputTree->Branch("track_numSCT"   ,&m_trk_numSCT);  
    // m_outputTree->Branch("track_numLSCT"  ,&m_trk_numLSCT);

    // The truth track parameters
    /*
    m_outputTree->Branch("eventNr", &m_eventNr);
    m_outputTree->Branch("t_loc0", &m_t_loc0);
    m_outputTree->Branch("t_loc1", &m_t_loc1);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_qop", &m_t_qop);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("truthMatched", &m_truthMatched);
    */
  
}

RootJetWriter::~RootJetWriter() {
  m_outputFile->Close();
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
                                          const TrackJetContainer& tp) {

// Exclusive access to all the writing procedure
std::lock_guard<std::mutex> lock(m_writeMutex);
  
Clear();

if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }
  
// Read input collections
// FIX: remove truth jets and keep only track jets!!
const auto& jets = m_inputJets(ctx);
const auto& tracks = m_inputTracks(ctx);
const auto& trackJets = m_inputTrackJets(ctx);

ACTS_DEBUG("RootWriter::Number of " << m_cfg.inputTracks << " " << tracks.size());
ACTS_DEBUG("RootWriter::Number of " << m_cfg.inputJets << " " << jets.size());
ACTS_DEBUG("RootWriter::Number of " << m_cfg.inputTrackJets << " " << trackJets.size());


  for (size_t ijets = 0; ijets < jets.size(); ++ijets) {

    //Only store jets with tracks associated to them and with label > -99
    //For now we don't have a label for the jets, so we skip this check

    // int jLabel = static_cast<int>(jets[ijets].getLabel());
    // int nTracksOnJet =  jets[ijets].getTracks().size();

    // if (jLabel < 0 || nTracksOnJet < 1 ) {
    //   continue;
    // }
    
    Acts::Vector4 jet_4mom = jets[ijets].getFourMomentum();
    Acts::Vector3 jet_3mom{jet_4mom[0],jet_4mom[1],jet_4mom[2]};
    float jet_theta = theta(jet_3mom);
    m_jet_pt.push_back(perp(jet_4mom));
    //m_jet_eta.push_back(jets[ijets].eta());
    m_jet_eta.push_back(std::atanh(std::cos(jet_theta)));
    m_jet_phi.push_back(phi(jet_4mom));

    m_jet_ncomponents.push_back(jets[ijets].getConstituents().size());
    m_jet_components.push_back(jets[ijets].getConstituents());
    m_jet_tracks_idx.push_back(jets[ijets].getTracks());
    
    m_jet_label.push_back(static_cast<int>(jets[ijets].getLabel()));
    
    m_jet_isPU.push_back(0);
    m_jet_isHS.push_back(1);

    std::vector<int> jetTracks = jets[ijets].getTracks();
    // ACTS_DEBUG("Jet " << ijets << " with pt: " << m_jet_pt[ijets] << " GeV"
    //            << " with label: " << m_jet_label[ijets]
    //            << " has " << jetTracks.size() << " tracks associated");

    
  } // jets


    //Do not fill the event if there is a jet without tracks associated to it
    bool keepEvent = true;
    for (auto links : m_jet_tracks_idx) {
      if (links.size() == 0)
        keepEvent = false;
      break;
    }

Acts::Vector3 vertexPosition{0., 0., 0.};
Acts::Vector4 stddev;
stddev[Acts::ePos0] = 10 * Acts::UnitConstants::um;
stddev[Acts::ePos1] = 10 * Acts::UnitConstants::um;
stddev[Acts::ePos2] = 75 * Acts::UnitConstants::um;
stddev[Acts::eTime] = 1 * Acts::UnitConstants::ns;
Acts::SquareMatrix4 vertexCov = stddev.cwiseProduct(stddev).asDiagonal();

Acts::Vertex ip_vtx(vertexPosition);


  Acts::ImpactPointEstimator::State state{magFieldCache()};

// Loop over the tracks
for (size_t itrk = 0; itrk < tracks.size(); itrk++) {

    double signed_d0             = -9999;
    double signed_z0SinTheta     = -9999;
    double signed_d0_err         = 1.;
    double signed_z0SinTheta_err = 1.;

    double covd0         = -999;
    double covz0         = -999;
    double covphi        = -999;
    double covtheta      = -999;
    double covqOverP     = -999;
    double covd0z0       = -999;
    double covd0phi      = -999;
    double covd0theta    = -999;
    double covd0qOverP   = -999;
    double covz0phi      = -999;
    double covz0theta    = -999;
    double covz0qOverP   = -999;
    double covphitheta   = -999;
    double covphiqOverP  = -999;
    double covthetaqOverP= -999;
        

    const auto trk = tracks.at(itrk);
    // const auto hit_infos  = hitInfos[itrk];
    const auto params     = trk.parameters();

    auto boundParams = trk.createParametersAtReference();



    //Check if this track belongs to a jet and compute the IPs
    for (size_t ijet = 0; ijet<trackJets.size(); ++ijet) {
      std::vector<int> jtrks = trackJets[ijet].getTracks();
      
      
      if (std::find(jtrks.begin(), jtrks.end(),itrk) != jtrks.end()) {

        Acts::Vector3 jetDir{jets[ijet].getFourMomentum()[0],
          trackJets[ijet].getFourMomentum()[1],
          trackJets[ijet].getFourMomentum()[2]};
      
      auto ipAndSigma = m_ipEst->estimate3DImpactParameters(gctx_, mctx_, boundParams, vertexPosition, state);
      auto vszs = m_ipEst->getLifetimeSignOfTrack(boundParams, ip_vtx, jetDir, gctx_, mctx_);
      
    
        if (!ipAndSigma.ok() || !vszs.ok()) {
          continue;
        }

          signed_d0 = std::fabs((*ipAndSigma).parameters()(0))  * (*vszs).first;
          auto z0 = std::fabs((*ipAndSigma).parameters()(1));
          auto theta = (*ipAndSigma).theta();
          signed_z0SinTheta = z0 * std::sin(theta) * (*vszs).second;
          
          
        //This is not unbiased!
        // signed_d0             = std::fabs((*ipAndSigma).IPd0)  * (*vszs).first;
        // signed_z0SinTheta     = std::fabs((*ipAndSigma).IPz0SinTheta) * (*vszs).second;
        // signed_d0_err         = (*ipAndSigma).sigmad0;
        // signed_z0SinTheta_err = (*ipAndSigma).sigmaz0SinTheta;

      }//track in jet
    }//loop on jets

// FIX: Still need to add vertex information to check covariances

//const auto& cov  = *trk.covariance();
      
      float trk_theta = params[Acts::eBoundTheta];
      float trk_eta   = std::atanh(std::cos(trk_theta));
      float trk_qop   = params[Acts::eBoundQOverP];
      float trk_p     = std::abs(1.0 / trk_qop);
      float trk_pt    = trk_p * std::sin(trk_theta); 
      
      m_tracks_prob.push_back(1.);   // todo
      m_trk_d0.push_back(params[Acts::eBoundLoc0]);             
      m_trk_z0.push_back(params[Acts::eBoundLoc1]);
      
      m_trk_signed_d0.push_back(signed_d0);
      //m_trk_signed_d0sig.push_back(signed_d0 / signed_d0_err);
      m_trk_signed_z0sinTheta.push_back(signed_z0SinTheta);
      //m_trk_signed_z0sinThetasig.push_back(signed_z0SinTheta / signed_z0SinTheta_err);

      // m_trk_numPix1L.push_back(hit_infos.nPixInnermost);
      // m_trk_numPix2L.push_back(hit_infos.nPixNextToInnermost);
      // m_trk_numPix  .push_back(hit_infos.nPix);
      // m_trk_numSCT  .push_back(hit_infos.nSStrip);
      // m_trk_numLSCT .push_back(hit_infos.nLStrip);
      
      
      m_trk_eta.push_back(trk_eta);            
      m_trk_theta.push_back(trk_theta);          
      m_trk_phi.push_back(params[Acts::eBoundPhi]);            
      m_trk_pt.push_back(trk_pt);             
      m_trk_qOverP.push_back(trk_p);         
      m_trk_t.push_back(params[Acts::eBoundTime]);              
      m_trk_t30.push_back(1.);            //todo
      m_trk_t60.push_back(1.);            //todo
      m_trk_t90.push_back(1.);            //todo
      m_trk_t120.push_back(1.);           //todo
      m_trk_t180.push_back(1.);           //todo
      m_trk_z.push_back(1.);              //todo
      
      // //This give un-initialized warnings
      
      // covd0          = cov(Acts::eBoundLoc0,Acts::eBoundLoc0);
      // covz0          = cov(Acts::eBoundLoc1,Acts::eBoundLoc1);
      // covphi         = cov(Acts::eBoundPhi,Acts::eBoundPhi);
      // covtheta       = cov(Acts::eBoundTheta,Acts::eBoundTheta);
      // covqOverP      = cov(Acts::eBoundQOverP,Acts::eBoundQOverP);
      // covd0z0        = cov(Acts::eBoundLoc0,Acts::eBoundLoc1);
      // covd0phi       = cov(Acts::eBoundLoc0,Acts::eBoundPhi);
      // covd0theta     = cov(Acts::eBoundLoc0,Acts::eBoundTheta);
      // covd0qOverP    = cov(Acts::eBoundLoc0,Acts::eBoundQOverP);
      // covz0phi       = cov(Acts::eBoundLoc1,Acts::eBoundPhi);
      // covz0theta     = cov(Acts::eBoundLoc1,Acts::eBoundTheta);
      // covz0qOverP    = cov(Acts::eBoundLoc1,Acts::eBoundQOverP);
      // covphitheta    = cov(Acts::eBoundPhi,Acts::eBoundTheta);
      // covphiqOverP   = cov(Acts::eBoundPhi,Acts::eBoundQOverP);
      // covthetaqOverP = cov(Acts::eBoundTheta,Acts::eBoundQOverP);
            
      
      // m_trk_var_d0.push_back(covd0);
      // m_trk_var_z0.push_back(covz0);         
      // m_trk_var_phi.push_back(covphi);        
      // m_trk_var_theta.push_back(covtheta);      
      // m_trk_var_qOverP.push_back(covqOverP);     
      // m_trk_cov_d0z0.push_back(covd0z0);       
      // m_trk_cov_d0phi.push_back(covd0phi);      
      // m_trk_cov_d0theta.push_back(covd0theta);    
      // m_trk_cov_d0qOverP.push_back(covd0qOverP);   
      // m_trk_cov_z0phi.push_back(covz0phi);      
      // m_trk_cov_z0theta.push_back(covz0theta);    
      // m_trk_cov_z0qOverP.push_back(covz0qOverP);   
      // m_trk_cov_phitheta.push_back(covphitheta);   
      // m_trk_cov_phiqOverP.push_back(covphiqOverP);  
      // m_trk_cov_thetaqOverP.push_back(covthetaqOverP);
      

}// loop over tracks



// const auto& inputTrajectories = m_inputTrajectories(ctx);
// const auto& reco_vertices = m_recoVertices(ctx);


    m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

void ActsExamples::RootJetWriter::Clear() {
    //Vertices
//   m_vtx_x.clear();
//   m_vtx_y.clear();
//   m_vtx_z.clear();
//   m_vtx_t.clear();
//   m_vtx_sumPt2.clear();
//   m_vtx_isHS.clear();
//   m_vtx_isPU.clear();
  
  //Jets
  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_ncomponents.clear();
  m_jet_components.clear();
  m_jet_tracks_idx.clear();
  m_jet_isPU.clear();
  m_jet_isHS.clear();
  m_jet_label.clear();

//   //Tracks
//   m_tracks_prob.clear();        
//   m_trk_d0.clear();             
//   m_trk_z0.clear();
  m_trk_signed_d0.clear();
//   m_trk_signed_d0sig.clear();
  m_trk_signed_z0sinTheta.clear();
//  m_trk_signed_z0sinThetasig.clear();
//   m_trk_eta.clear();            
//   m_trk_theta.clear();          
//   m_trk_phi.clear();            
//   m_trk_pt.clear();
//   m_trk_qOverP.clear();             
//   m_trk_t.clear();              
//   m_trk_t30.clear();            
//   m_trk_t60.clear();            
//   m_trk_t90.clear();            
//   m_trk_t120.clear();           
//   m_trk_t180.clear();           
//   m_trk_z.clear();              
//   m_trk_var_d0.clear();         
//   m_trk_var_z0.clear();         
//   m_trk_var_phi.clear();        
//   m_trk_var_theta.clear();      
//   m_trk_var_qOverP.clear();     
//   m_trk_cov_d0z0.clear();       
//   m_trk_cov_d0phi.clear();      
//   m_trk_cov_d0theta.clear();    
//   m_trk_cov_d0qOverP.clear();   
//   m_trk_cov_z0phi.clear();      
//   m_trk_cov_z0theta.clear();    
//   m_trk_cov_z0qOverP.clear();   
//   m_trk_cov_phitheta.clear();   
//   m_trk_cov_phiqOverP.clear();  
//   m_trk_cov_thetaqOverP.clear();

//   m_trk_numPix1L.clear();
//   m_trk_numPix2L.clear();
//   m_trk_numPix  .clear();
//   m_trk_numSCT  .clear();
//   m_trk_numLSCT .clear();
}

