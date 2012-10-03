// -*- C++ -*-
//
// Package:    EDMElectronSelector
// Class:      EDMElectronSelector
// 
/**\class EDMElectronSelector EDMElectronSelector.cc EDMUtilities/EDMElectronSelector/src/EDMElectronSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kaufman Nicolas
//         Created:  Sun Sep 30 12:46:48 CDT 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/Polynomial.h"
#include "Math/GenVector/Boost.h"
#include "TMath.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

//
// class declaration
//

class EDMElectronSelector : public edm::EDProducer {
   public:
      explicit EDMElectronSelector(const edm::ParameterSet&);
      ~EDMElectronSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      virtual double electronIsolation(const pat::Electron& electron);

      edm::InputTag electronSrc_          ;
      std::string electronIsolationType_  ;
      double electronIsolationMax_        ;

      // ----------member data ---------------------------
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
EDMElectronSelector::EDMElectronSelector(const edm::ParameterSet& iConfig)
{
  electronSrc_            = iConfig.getParameter<edm::InputTag>( "electronSrc" );
  electronIsolationType_  = iConfig.getParameter<std::string>( "electronIsolationType" );
  electronIsolationMax_   = iConfig.getParameter<double>( "maximumElectronIsolation" );

  produces<std::vector<pat::Electron> >();

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


EDMElectronSelector::~EDMElectronSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EDMElectronSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<View<pat::Electron> > inputElectrons;
   iEvent.getByLabel( electronSrc_ , inputElectrons );
   if (inputElectrons.failedToGet()) cout <<"not found some electrons"<<endl;

   std::vector<pat::Electron> outputElectrons;

   //cout << "beginning loop over electrons" << endl;

   for(unsigned int iElectron = 0; iElectron < inputElectrons->size(); iElectron++)
     {
       pat::Electron electron (inputElectrons->at(iElectron));
       double iso = electronIsolation( electron );
       //cout << "electron isolation: " << iso << endl;
       if (iso < electronIsolationMax_ ) {
	 //cout << "keeping this electron " << endl;
	 outputElectrons.push_back( electron );
       }
     }

   //cout << "number of electrons passing reliso cut: " << outputElectrons.size() << endl;

   auto_ptr<std::vector<pat::Electron> > selectedElectrons( new std::vector<pat::Electron>(outputElectrons) );
   iEvent.put( selectedElectrons );

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

double
EDMElectronSelector::electronIsolation(const pat::Electron& electron)
{

  reco::isodeposit::AbsVetos vetos_ch;
  reco::isodeposit::AbsVetos vetos_nh;
  reco::isodeposit::AbsVetos vetos_ph;
  reco::isodeposit::Direction Dir = reco::isodeposit::Direction(electron.superCluster()->eta(), electron.superCluster()->phi());
  if( abs( electron.superCluster()->eta() ) > 1.479 ){
    vetos_ch.push_back(new reco::isodeposit::ConeVeto( Dir, 0.015 ));
    vetos_ph.push_back(new reco::isodeposit::ConeVeto( Dir, 0.08 ));
  }
  //cone size 0.3
  const double chIso03 = electron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos_ch).first;
  const double nhIso03 = electron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos_nh).first;
  const double phIso03 = electron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos_ph).first;

  const double puChIso03 = electron.isoDeposit(pat::PfPUChargedHadronIso)->depositAndCountWithin(0.3, vetos_ch).first;

  const double relIso = ( chIso03 + nhIso03 + phIso03 ) / electron.pt() ;
  double maxDepIso03 = (nhIso03 + phIso03 - 0.5*puChIso03 > 0.0) ? nhIso03 + phIso03 - 0.5*puChIso03 : 0.0;
  const double relIsodb = ( chIso03 + maxDepIso03 ) / electron.pt();

  //cout << "The relative isolation is " << relIso << endl;

  if(electronIsolationType_ == "relIso") return relIso ;
  if(electronIsolationType_ == "relIsodb") return relIsodb ;
  return 0.;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EDMElectronSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EDMElectronSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
EDMElectronSelector::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EDMElectronSelector::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EDMElectronSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EDMElectronSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EDMElectronSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDMElectronSelector);
