// -*- C++ -*-
//
// Package:    Analysis/TriggerStudies
// Class:      TriggerStudies
// 
/**\class TriggerStudies TriggerStudies.cc Analysis/TriggerStudies/plugins/TriggerStudies.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Tue, 02 Jun 2015 11:49:34 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/Utilities.h"

#include <sstream>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

//
// class declaration
//

class TriggerStudies : public edm::EDAnalyzer {
  public:
    explicit TriggerStudies(const edm::ParameterSet&);
    ~TriggerStudies();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual bool passJetIDLoose (const double& chf, const double& nhf, const double& cef, const double& nef, 
        const int& nch, const int& nconstituents) ; 

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    std::string hltProcName_ ; 
    edm::InputTag ak8jetLabel_; 
    edm::InputTag ak4jetLabel_; 
    std::vector<std::string> hltpaths_ ; 

    HLTConfigProvider hltConfig_;
    int triggerBit_;
    std::vector<std::string>hltPaths_ ; 

    edm::Service<TFileService> fs ; 
    std::map<std::string, TH1D*> h1_ ; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerStudies::TriggerStudies(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))), 
  hltProcName_(iConfig.getParameter<std::string>("hltProcName")), 
  ak8jetLabel_ (iConfig.getParameter<edm::InputTag>("ak8jetLabel")),
  ak4jetLabel_ (iConfig.getParameter<edm::InputTag>("ak4jetLabel")), 
  hltPaths_ (iConfig.getParameter<std::vector<std::string>>    ("hltPaths"))  
{
  //now do what ever initialization is needed

}


TriggerStudies::~TriggerStudies()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TriggerStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  //// Event collections
  vlq::JetCollection goodAK8Jets, goodAK4Jets, btaggedlooseAK4, btaggedmediumAK4 ;

  //// Event selection variables
  double ptak8leading(0) ; 
  double ptak82nd(0) ;
  double etaforwardmostak4(0) ; 
  double htak4(0) ; 
  double ptak4bjetleading(0) ; 

  bool changedConfig = false;

  InputTag jLabel_; 
  if (!hltConfig_.init(iEvent.getRun(), iSetup, hltProcName_, changedConfig)) {
    std::cout << ">>>>>>ERROR: Initialization of HLTConfigProvider failed!!" << std::endl;
    return ;
  }

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  //std::cout << " Trigger bit size = " << triggerBits->size() << std::endl ; 

  //for (unsigned int ii = 0, nn = triggerBits->size(); ii < nn; ++ii) {
  //  std::cout << ">>>>>>HLT bit " << ii << " path: " <<  (hltConfig_.triggerNames()[ii]) << " pass? " << triggerBits->accept(ii)  
  //    << std::endl ; 
  //  //<< " has prescale index " << triggerPrescales->getPrescaleForIndex(ii) << std::endl ; 
  //}

  //// Do event selection
  edm::Handle<std::vector<pat::Jet> > ak8jetHandle, ak4jetHandle;
  iEvent.getByLabel(ak8jetLabel_, ak8jetHandle);
  iEvent.getByLabel(ak4jetLabel_, ak4jetHandle);

  std::auto_ptr<vector<pat::Jet> > ak8jetColl( new vector<pat::Jet> (*ak8jetHandle) );
  for (size_t i = 0; i< ak8jetColl->size(); i++){
    pat::Jet & jet = (*ak8jetColl)[i];

    int idx = i;
    double pt = jet.pt() ; 
    double eta = jet.eta() ; 
    double rapidity = jet.rapidity() ; 
    double phi = jet.phi() ; 
    double energy = jet.energy() ; 
    double mass = jet.mass() ; 

    if (pt < 300 || abs(eta) > 2.5 || mass < 50) continue ; 

    double bdiscCSV =  jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") ; 
    double softdropmass = jet.userFloat("softDropMass") ; 
    double tau1 = jet.userFloat("userFloat('NjettinessAK8:tau1')") ; 
    double tau2 = jet.userFloat("userFloat('NjettinessAK8:tau2')") ; 

    TLorentzVector jetp4 ; 
    jetp4.SetPtEtaPhiM(pt, eta, phi, mass) ; 
    vlq::Jet myjet;  
    myjet.setP4(jetp4) ;
    myjet.setCSV(bdiscCSV) ;
    goodAK8Jets.push_back(myjet) ;

  }

  std:: auto_ptr<vector<pat::Jet> > ak4jetColl( new vector<pat::Jet> (*ak4jetHandle) );
  for (size_t i = 0; i< ak4jetColl->size(); i++){
    pat::Jet & jet = (*ak4jetColl)[i];

    int idx = i;
    double pt = jet.pt() ; 
    double eta = jet.eta() ; 

    if (pt < 50 || abs(eta) > 5) continue ; 

    double rapidity = jet.rapidity() ; 
    double phi = jet.phi() ; 
    double energy = jet.energy() ; 
    double mass = jet.mass() ; 
    double bdiscCSV = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") ; 

    TLorentzVector jetp4 ; 
    jetp4.SetPtEtaPhiM(pt, eta, phi, mass) ; 
    vlq::Jet myjet;  
    myjet.setP4(jetp4) ;
    myjet.setCSV(bdiscCSV) ;
    goodAK4Jets.push_back(myjet) ;
    if (bdiscCSV > 0.605 && abs(eta) < 2.4) btaggedlooseAK4.push_back(myjet) ; 
    if (bdiscCSV > 0.890 && abs(eta) < 2.4) btaggedmediumAK4.push_back(myjet) ; 

  }

  //// Event selection NAK8jets >= 2
  if ( goodAK8Jets.size() < 2 ) return ; 

  //// Event selection NAK4jets >= 5
  if ( goodAK4Jets.size() < 5 ) return ; 

  //// Event selection Nbjets >= 3
  if ( btaggedmediumAK4.size() < 1 || btaggedlooseAK4.size() < 3 ) return ; 

  //// Forwardmost jet eta > 2.5
  std::sort(goodAK4Jets.begin(), goodAK4Jets.end(), Utilities::sortByEta<vlq::Jet>) ; 
  etaforwardmostak4 = abs((goodAK4Jets.at(0)).getEta()) ; 
  if ( etaforwardmostak4 < 2.5 ) return ; 

  std::sort(goodAK4Jets.begin(), goodAK4Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 
  std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 
  std::sort(btaggedmediumAK4.begin(), btaggedmediumAK4.end(), Utilities::sortByPt<vlq::Jet>) ; 
  
  /// Event selection HT  
  HT HTAK4(goodAK4Jets) ; 
  htak4 = HTAK4.getHT() ;
  if ( htak4 < 600 ) return ; 

  ptak8leading = (goodAK8Jets.at(0)).getPt() ; 
  h1_["ptak8leading"] -> Fill(ptak8leading) ; 
  ptak82nd = (goodAK8Jets.at(1)).getPt() ; 
  h1_["ptak82nd"] -> Fill(ptak82nd) ; 
  h1_["ht"] -> Fill(htak4) ; 
  ptak4bjetleading = (btaggedmediumAK4.at(0)).getPt() ;
  h1_["ptak4bjetleading"] -> Fill(ptak4bjetleading) ; 

  ////  Get HLT decisions
  for ( const std::string& myhltpath : hltPaths_ ) {
    for (unsigned int ii = 0, nn = triggerBits->size(); ii < nn; ++ii) { 
      std::string trigname = (hltConfig_.triggerNames()[ii]) ; 
      bool hltdecision = triggerBits->accept(ii) ; 
      if ( myhltpath == trigname && hltdecision ) {
        std::stringstream ss ;
        ss << "ptak8leading_" << myhltpath ; 
        h1_[ss.str()] -> Fill(ptak8leading) ; 

        ss.clear() ; 
        ss.str("") ; 
        ss << "ptak82nd_" << myhltpath ; 
        h1_[ss.str()] -> Fill (ptak82nd) ;  

        ss.clear() ; 
        ss.str("") ; 
        ss << "ht_" << myhltpath ; 
        h1_[ss.str()] -> Fill (htak4) ;  

        ss.clear() ; 
        ss.str("") ; 
        ss << "ptak4bjetleading_" << myhltpath ; 
        h1_[ss.str()] -> Fill (ptak4bjetleading) ;  

      }
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void TriggerStudies::beginJob() {

  h1_["ptak8leading"]  = fs->make<TH1D>("ptak8leading"  ,";p_T(leading AK8 jet) [GeV];;" , 40, 0., 2000.) ; 
  h1_["ptak82nd"]  = fs->make<TH1D>("ptak82nd"  ,";p_T(2nd AK8 jet) [GeV];;" , 40, 0., 2000.) ; 
  h1_["ht"] = fs->make<TH1D>("ht" ,";H_T (AK4 jets) [GeV]", 50, 0., 2000.) ; 
  h1_["ptak4bjetleading"] = fs->make<TH1D>("ptak4bjetleading"  ,";p_T(leading AK4 b jet medium OP) [GeV];;" , 100, 0., 1000.) ; 

  for ( const std::string& myhltpath : hltPaths_ ) {
    std::stringstream ss ;
    ss << "ptak8leading_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_T (leading AK8 jet) [GeV];;" ,40, 0., 2000.) ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak82nd_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_T (2nd AK8 jet) [GeV];;" ,40, 0., 2000.) ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ht_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,"H_T (AK4 jets) [GeV]" ,50, 0., 2000.) ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak4bjetleading_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_T (leading AK4 b jet medium OP) [GeV];;" ,100, 0., 1000.) ; 
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerStudies::endJob() {

  for ( const std::string& myhltpath : hltPaths_ ) {
    std::stringstream ss ;
    ss << "ptak8leading_" << myhltpath ; 
    TGraphAsymmErrors* greffak8leading = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["ptak8leading"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_ptak8leading_" << myhltpath ;
    greffak8leading->SetName((ss.str()).c_str()) ;  
    greffak8leading->Write() ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak82nd_" << myhltpath ; 
    TGraphAsymmErrors* greffak82nd = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["ptak82nd"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_ptak82nd_" << myhltpath ;
    greffak82nd->SetName((ss.str()).c_str()) ; 
    greffak82nd->Write() ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ht_" << myhltpath ; 
    TGraphAsymmErrors* greffht = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["ht"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_ht_" << myhltpath ;
    greffht->SetName((ss.str()).c_str()) ; 
    greffht->Write() ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak4bjetleading_" << myhltpath ; 
    TGraphAsymmErrors* greffak4bjetleading = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["ptak4bjetleading"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_ptak4bjetleading_" << myhltpath ;
    greffak4bjetleading->SetName((ss.str()).c_str()) ; 
    greffak4bjetleading->Write() ; 
  }

}

bool TriggerStudies::passJetIDLoose (const double& chf, const double& nhf, const double& cef, const double& nef, const int& nch, const int& nconstituents) {
  bool decision(false) ; 

  double _chf = 0.0;
  double _nhf = 0.99;
  double _cef = 0.99;
  double _nef = 0.99;
  int    _nch = 0;
  int    _nconstituents = 1;

  if (chf > _chf && nhf < _nhf && cef < _cef && nef < _nef && nch > _nch && nconstituents > _nconstituents) decision = true ;

  return decision ;
}

// ------------ method called when starting to processes a run  ------------
/*
   void TriggerStudies::beginRun(edm::Run const&, edm::EventSetup const&) {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void  TriggerStudies::endRun(edm::Run const&, edm::EventSetup const&) {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void  TriggerStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void  TriggerStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerStudies);
