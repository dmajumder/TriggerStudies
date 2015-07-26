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
    virtual bool filter(const edm::Event&, const std::string, const double, const int) ; 

    virtual bool passJetIDLoose (const double& chf, const double& nhf, const double& cef, const double& nef, 
        const int& nch, const int& nconstituents) ; 

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

    std::string origPath_;
    std::vector<double> newThreshs_;
    std::vector<int> triggerTypes_;

    std::string hltProcName_ ; 
    edm::InputTag ak8jetLabel_; 
    edm::InputTag ak4jetLabel_; 
    std::vector<std::string> hltpaths_ ; 

    HLTConfigProvider hltConfig_;
    int triggerBit_;
    std::vector<std::string>hltPaths_ ; 

    edm::Service<TFileService> fs ; 
    std::map<std::string, TH1D*> h1_ ; 

    static const float ptbins_[13] ; 
    static const float ptbbins_[19] ; 
    static const float htbins_[19] ; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
const float TriggerStudies::ptbins_[13] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 900., 1200.} ; 
const float TriggerStudies::ptbbins_[19] = {0., 25., 50., 75., 100., 125., 150., 175., 200., 225., 250., 275., 300., 325., 350., 375., 400., 500., 800.} ;
const float TriggerStudies::htbins_[19] = {200., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1100., 1200., 1400., 1600., 2000.} ; 

template <typename T>
struct iterator_extractor { typedef typename T::iterator type; };

template <typename T>
struct iterator_extractor<T const> { typedef typename T::const_iterator type; };


template <typename T>
class Indexer {
  public:
    class iterator {
      typedef typename iterator_extractor<T>::type inner_iterator;

      typedef typename std::iterator_traits<inner_iterator>::reference inner_reference;
      public:
      typedef std::pair<size_t, inner_reference> reference;

      iterator(inner_iterator it): _pos(0), _it(it) {}

      reference operator*() const { return reference(_pos, *_it); }

      iterator& operator++() { ++_pos; ++_it; return *this; }
      iterator operator++(int) { iterator tmp(*this); ++*this; return tmp; }

      bool operator==(iterator const& it) const { return _it == it._it; }
      bool operator!=(iterator const& it) const { return !(*this == it); }

      private:
      size_t _pos;
      inner_iterator _it;
    };

    Indexer(T& t): _container(t) {}

    iterator begin() const { return iterator(_container.begin()); }
    iterator end() const { return iterator(_container.end()); }

  private:
    T& _container;
}; // class Indexer

template <typename T>
Indexer<T> index(T& t) { return Indexer<T>(t); }

//
// constructors and destructor
//
TriggerStudies::TriggerStudies(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))), 
  origPath_(iConfig.getParameter<std::string>("origpath")),//original path to use as a base
  newThreshs_(iConfig.getParameter<std::vector<double>>("newthreshs")),//new threshold to use on top of original path
  triggerTypes_(iConfig.getParameter<std::vector<int>>("triggertypes")),//trigger type used
  hltProcName_(iConfig.getParameter<std::string>("hltProcName")), 
  ak8jetLabel_ (iConfig.getParameter<edm::InputTag>("ak8jetLabel")),
  ak4jetLabel_ (iConfig.getParameter<edm::InputTag>("ak4jetLabel")), 
  hltPaths_ (iConfig.getParameter<std::vector<std::string>>    ("hltPaths"))  
{
  //now do what ever initialization is needed
  if (hltPaths_.size() != triggerTypes_.size() && triggerTypes_.size() != newThreshs_.size() ) 
    edm::LogError("SizeMismatch") << ">>>>> Sizes of hltPaths_ " << hltPaths_.size() 
      << " triggerTypes_ " << triggerTypes_.size() << " newThreshs_" << newThreshs_.size() << " do not match\n" ; 

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
  double ptak4bjetleading(0) ; 
  double htak4(0) ; 
  double htak8(0) ; 

  //// Event selection flags
  bool ispassedLeadingAK8Pt(false) ; 
  bool ispassed2ndAK8Pt(false) ; 
  bool ispassedBTag(false) ; 
  bool ispassedHT(false) ; 

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
    double softdropmass = jet.userFloat("ak8PFJetsCHSSoftDropMass") ; 

    if (abs(eta) > 2.5 || softdropmass < 50) continue ; 

    double tau1 = jet.userFloat("userFloat('NjettinessAK8:tau1')") ; 
    double tau2 = jet.userFloat("userFloat('NjettinessAK8:tau2')") ; 

    if (tau2/tau1 > 0.5) continue ; 

    double bdiscCSV =  jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") ; 

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

    if (pt < 30 || abs(eta) > 5) continue ; 

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

  //// Event pre-selection 
  if ( goodAK8Jets.size() < 1 ) return ; 
  if ( goodAK4Jets.size() < 4 ) return ; 

  if (goodAK4Jets.size() > 0) std::sort(goodAK4Jets.begin(), goodAK4Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 
  if (goodAK8Jets.size() > 0) std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 
  if (btaggedmediumAK4.size() > 0) std::sort(btaggedmediumAK4.begin(), btaggedmediumAK4.end(), Utilities::sortByPt<vlq::Jet>) ; 

  if (goodAK8Jets.size() > 0)  ptak8leading = (goodAK8Jets.at(0)).getPt() ; 
  if (goodAK8Jets.size() > 1) ptak82nd = (goodAK8Jets.at(1)).getPt() ;  
  if (btaggedmediumAK4.size() > 0)  ptak4bjetleading = (btaggedmediumAK4.at(0)).getPt() ;
  HT HTAK4(goodAK4Jets) ; 
  htak4 = HTAK4.getHT() ;
  HT HTAK8(goodAK8Jets) ; 
  if  (goodAK8Jets.size() > 1) htak8 = HTAK8.getHT() ;

  //// Event selection flags set 
  if (ptak8leading > 300. ) ispassedLeadingAK8Pt = true ;
  if (ptak82nd > 250.) ispassed2ndAK8Pt = true ;  
  if ( btaggedmediumAK4.size() > 0 || btaggedlooseAK4.size() > 2 ) ispassedBTag = true ; 
  if (htak4 > 700) ispassedHT = true ; 

  //// Fill "N-1" plots 
  if (ispassedBTag && ispassedHT ) h1_["ptak8leading"] -> Fill(ptak8leading) ; 
  if (ispassedLeadingAK8Pt && ispassedBTag && ispassedHT ) h1_["ptak82nd"] -> Fill(ptak82nd) ; 
  if (ispassedBTag) h1_["htak4"] -> Fill(htak4) ; 
  if (goodAK8Jets.size() > 1 && ispassedBTag) h1_["htak8"] -> Fill(htak8) ; 
  if (ispassedLeadingAK8Pt && ispassed2ndAK8Pt && ispassedBTag && ispassedHT) h1_["ptak4bjetleading"] -> Fill(ptak4bjetleading) ; 

  ////  Get HLT decisions

  for ( auto hltpath : index(hltPaths_) ) {

    const std::string myhltpath = hltpath.second ;  
    const int triggerType = triggerTypes_.at(hltpath.first) ; 
    const double newThresh = newThreshs_.at(hltpath.first) ; 

    //std::cout << " hltpath = " << myhltpath << " triggertype = " << triggerType << std::endl ;

    if ( !filter(iEvent, myhltpath, newThresh, triggerType) ) continue ; 

    std::stringstream ss ;

    if (ispassedBTag && ispassedHT )  {
      ss.clear() ; 
      ss.str("") ; 
      ss << "ptak8leading_" << myhltpath ; 
      h1_[ss.str()] -> Fill(ptak8leading) ; 
    }

    if (ispassedLeadingAK8Pt && ispassedBTag && ispassedHT ) { 
      ss.clear() ; 
      ss.str("") ; 
      ss << "ptak82nd_" << myhltpath ; 
      h1_[ss.str()] -> Fill (ptak82nd) ;  
    }

    if (ispassedBTag) { 
      ss.clear() ; 
      ss.str("") ; 
      ss << "htak4_" << myhltpath ; 
      h1_[ss.str()] -> Fill (htak4) ;  
    }

    if ( goodAK8Jets.size() > 1 && ispassedBTag ) {
      ss.clear() ; 
      ss.str("") ; 
      ss << "htak8_" << myhltpath ; 
      h1_[ss.str()] -> Fill (htak8) ;  
    }

    if (ispassedLeadingAK8Pt && ispassed2ndAK8Pt && ispassedBTag && ispassedHT) { 
      ss.clear() ; 
      ss.str("") ; 
      ss << "ptak4bjetleading_" << myhltpath ; 
      h1_[ss.str()] -> Fill (ptak4bjetleading) ;  
    }

  }

  ////  Get HLT decisions
  //DMfor ( const std::string& origPath_ : hltPaths_ ) {
  //DM  for (unsigned int ii = 0, nn = triggerBits->size(); ii < nn; ++ii) { 
  //DM    std::string trigname = (hltConfig_.triggerNames()[ii]) ; 
  //DM    bool hltdecision = triggerBits->accept(ii) ; 
  //DM    if ( origPath_ == trigname && hltdecision ) {

  //DM    }
  //DM  }
  //DM}

}


// ------------ method called once each job just before starting event loop  ------------
void TriggerStudies::beginJob() {

  //h1_["ptak8leading"]  = fs->make<TH1D>("ptak8leading"  ,";p_{T}(leading AK8 jet) [GeV];;" , 40, 0., 2000.) ; 
  //h1_["ptak82nd"]  = fs->make<TH1D>("ptak82nd"  ,";p_{T}(2nd AK8 jet) [GeV];;" , 40, 0., 2000.) ; 
  //h1_["htak4"] = fs->make<TH1D>("htak4" ,";H_{T} (AK4 jets) [GeV]", 50, 0., 2000.) ; 
  //h1_["htak8"] = fs->make<TH1D>("htak8" ,";H_{T} (AK8 jets) [GeV]", 50, 0., 2000.) ; 
  //h1_["ptak4bjetleading"] = fs->make<TH1D>("ptak4bjetleading"  ,";p_{T}(leading AK4 b jet medium OP) [GeV];;" , 100, 0., 1000.) ; 

  h1_["ptak8leading"]  = fs->make<TH1D>("ptak8leading"  ,";p_{T}(leading AK8 jet) [GeV];;" , 12, ptbins_ ) ;
  h1_["ptak82nd"]  = fs->make<TH1D>("ptak82nd"  ,";p_{T}(2nd AK8 jet) [GeV];;" , 12, ptbins_ ) ;
  h1_["htak4"] = fs->make<TH1D>("htak4" ,";H_{T} (AK4 jets) [GeV]", 18, htbins_ ) ;  
  h1_["htak8"] = fs->make<TH1D>("htak8" ,";H_{T} (AK8 jets) [GeV]", 18, htbins_ ) ;  
  h1_["ptak4bjetleading"] = fs->make<TH1D>("ptak4bjetleading"  ,";p_{T}(leading AK4 b jet medium OP) [GeV];;" , 18, ptbbins_ ) ; 

  for ( const std::string& myhltpath : hltPaths_ ) {
    std::stringstream ss ;
    ss << "ptak8leading_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_{T} (leading AK8 jet) [GeV];;" , 12, ptbins_ ) ;

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak82nd_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_{T} (2nd AK8 jet) [GeV];;" ,12, ptbins_ ) ;

    ss.clear() ; 
    ss.str("") ; 
    ss << "htak4_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,"H_{T} (AK4 jets) [GeV]" ,18, htbins_ ) ;

    ss.clear() ; 
    ss.str("") ; 
    ss << "htak8_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,"H_{T} (AK8 jets) [GeV]" ,18, htbins_ ) ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "ptak4bjetleading_" << myhltpath ; 
    h1_[ss.str()] = fs->make<TH1D>((ss.str()).c_str() ,";p_{T} (leading AK4 b jet medium OP) [GeV];;" ,18, ptbbins_ ) ; 
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
    ss << "htak4_" << myhltpath ; 
    TGraphAsymmErrors* greffhtak4 = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["htak4"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_htak4_" << myhltpath ;
    greffhtak4->SetName((ss.str()).c_str()) ; 
    greffhtak4->Write() ; 

    ss.clear() ; 
    ss.str("") ; 
    ss << "htak8_" << myhltpath ; 
    TGraphAsymmErrors* greffhtak8 = fs->make<TGraphAsymmErrors>(h1_[ss.str()], h1_["htak8"], "cp") ;
    ss.clear() ; 
    ss.str("") ; 
    ss << "eff_htak8_" << myhltpath ;
    greffhtak8->SetName((ss.str()).c_str()) ; 
    greffhtak8->Write() ; 

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

bool TriggerStudies::filter(const edm::Event& iEvent, const std::string myhltpath, const double newThresh, const int triggerType) {
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    //std::cout << ">>>> HLT path " << names.triggerName(i) << " decision = " << triggerBits->accept(i) << std::endl ;
    if (names.triggerName(i)==myhltpath && triggerBits->accept(i)) {
      //std::cout <<" Found path pass " << names.triggerName(i) << std::endl ;
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
        obj.unpackPathNames(names);
        for (unsigned h = 0; h < obj.filterIds().size(); ++h) {
          if (obj.filterIds()[h]==triggerType && obj.hasPathName( myhltpath, true, true )) { //look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h for an explanation of trigger types
            //std::cout << "Found correct type and path object... ";
            if (obj.pt()>newThresh) {
              //std::cout << "Found passing object!" << std::endl;
              return true;
            }
          }
        }
      }
    }
  }
  //std::cout << std::endl;

  return false;
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
