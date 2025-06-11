// All the user defined functions for the nTuplelizer are supposed to be stored here
// Include necessary header files

#ifndef ZEE_RecHit_NTuplizer_H
#define ZEE_RecHit_NTuplizer_H

#include <memory>
#include <fstream>
#include <iostream>
#include <random>
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TFile.h"
#include "correction.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"


#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <Math/Vector4D.h>

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "fastjet/ClusterSequence.hh"
#include "DataFormats/Math/interface/deltaR.h"
#include <algorithm>
#include "DataFormats/VertexReco/interface/Vertex.h"

//Namespaces

using namespace edm;
using namespace std;
using namespace reco;

using reco::TrackCollection;



// A customised implementation of the Cluster Lazy tools
class MyEcalClusterLazyTools : public noZS::EcalClusterLazyTools {
public:
    MyEcalClusterLazyTools(const edm::Event &ev,
                        ESData const &esData,
                        edm::EDGetTokenT<EcalRecHitCollection> token1,
                        edm::EDGetTokenT<EcalRecHitCollection> token2)
        : noZS::EcalClusterLazyTools(ev, esData, token1, token2) {}

    // Public method to access the protected method
    const EcalRecHitCollection* getRecHits(const BasicCluster& b_cluster) {
        return this->getEcalRecHitCollection(b_cluster);
    }
};


//Function to print jet information
void printJetInfo(const pat::Jet* jet) {
    // Print all available user floats
    std::cout << "Available user floats:" << std::endl;
    for (const auto& name : jet->userFloatNames()) {
        std::cout << "  - " << name << ": " << jet->userFloat(name) << std::endl;
	}

    // Print all available user integers
    std::cout << "Available user integers:" << std::endl;
    for (const auto& name : jet->userIntNames()) {
        std::cout << "  - " << name << ": " << jet->userInt(name) << std::endl;
    }

    // Accessing b-tagging info (example)
    std::cout << "Jet b-tagging discriminators:" << std::endl;
    for (const auto& discName : jet->getPairDiscri()) {
        std::cout << "  - " << discName.first << ": " << discName.second << std::endl;
    }
}


// Function to get dR values 
double dR(const reco::Candidate* a1, const reco::GenParticle* a2){
	return deltaR(*a1, *a2);
}

// Function to get dR between 4 vectors
double dR_4vec(ROOT::Math::PtEtaPhiEVector a1, ROOT::Math::PtEtaPhiEVector a2){
        return deltaR(a1.eta(), a1.phi(), a2.eta(), a2.phi());
}


// Function to reject jets which are too close to electrons
bool isJet(const ROOT::Math::PtEtaPhiEVector& a_1, 
                         const std::vector<ROOT::Math::PtEtaPhiEVector>& elec_vec, 
                         double tolerance = 0.1) {
    for (const auto& vec : elec_vec) {
        double dR_value = dR_4vec(a_1, vec);  
        if (dR_value <= tolerance) {
            return false;  
        }
    }
    return true;  
}

// Function to print all mothers of a GenParticle
void findMothers(const reco::GenParticle* particle, int depth =1){
        if(particle->numberOfMothers()!=0){
                //cout<<particle->pdgId();
                for(size_t i=0; i<particle->numberOfMothers(); i++){
                        const reco::Candidate* mother = particle->mother(i);
                        //cout<<std::string(depth * 2, '#')<<"<=="<<mother->pdgId();
                        cout<<"<=="<<mother->pdgId()<<" (Level: "<<depth<<")";
                        findMothers(dynamic_cast<const reco::GenParticle*>(mother), depth+1);
                        cout<<endl;
                }
        }
        else{
                cout<<"<==end";
        }
}


// Function to find all mothers of a Candidate class object
void findMothers_candidate(const reco::Candidate* particle, int depth =1){
        if(particle->numberOfMothers()!=0){
                //cout<<particle->pdgId();
                for(size_t i=0; i<particle->numberOfMothers(); i++){
                        const reco::Candidate* mother = particle->mother(i);
                        //cout<<std::string(depth * 2, '#')<<"<=="<<mother->pdgId();
                        cout<<"<=="<<mother->pdgId()<<" (Level: "<<depth<<")";
                        findMothers_candidate(dynamic_cast<const reco::Candidate*>(mother), depth+1);
                        cout<<endl;
                }
        }
        else{
                cout<<"<==end";
        }
}

// Function for dR matching between GenJet and RecoJet
double dR_gen_reco_jets(const reco::GenJet* gj1, const reco::PFJet* rj1){
  	return deltaR(*gj1, *rj1);	
}

//Function to get DressedLepton from Born and Naked
const reco::GenParticle* GetDressedLepton(
            const reco::GenParticle * const BORN_LEP,
            const reco::GenParticle * const NAKED_LEP,
            const double MAX_DELTA_R,
            const int PDGID
            ) {
        double MASS;
        if (fabs(PDGID) == 11){
        MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. GetDressedLepton will return a nullptr"<<endl;
                return nullptr;
        }
        // Make a 4 vector for the dressed electron
        math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return nullptr;
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == PDGID) {
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //cout<<"Status of photon: "<<test_particle->status()<<endl;
                    if (DELTA_R < MAX_DELTA_R) {
                        dressed_p4 += math::PtEtaPhiMLorentzVector(
                                test_particle->pt(),
                                test_particle->eta(),
                                test_particle->phi(),
                                MASS
                                );
                    }
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
            }
        }

        // Make a GenParticle from the vector and return it
        reco::GenParticle* dressed_lep = new reco::GenParticle(
                NAKED_LEP->charge(),
                dressed_p4,
                NAKED_LEP->vertex(),
                NAKED_LEP->pdgId(),
                NAKED_LEP->status(),
                1
            );

        return dressed_lep;
    }

//Function to return the total fsr energy of the leptons
double total_fsr(const reco::GenParticle* const BORN_LEP, const reco::GenParticle* const NAKED_LEP, const int PDGID){
        // Make a 4 vector for the dressed electron
        //double MASS;
        if (fabs(PDGID) == 11){
        //MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        //MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        //MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. Will return -999.0"<<endl;
                return -999.0;
        }
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        double energy=0;
        //assert(energy.size() == 0);
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return -999.0;
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == PDGID) {
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    //const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //cout<<"Status of photon: "<<test_particle->status()<<endl;
                    
                        energy=energy+test_particle->energy();
                    
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
            }
        }

        return energy;
    }

//Function to return the FSR energy vector of the leptons
std::vector<double> fsr_energy(
            const reco::GenParticle * const BORN_LEP,
            const reco::GenParticle * const NAKED_LEP,
            const int PDGID
            ) {
        //double MASS;
        if (fabs(PDGID) == 11){
        //MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        //MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        //MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. GetDressedLepton will return a nullptr"<<endl;
                return std::vector<double>();
        }
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        std::vector<double> rad_photon_energy;
        assert(rad_photon_energy.size() == 0);
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return std::vector<double>();
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                //cout<<"Number of daughters: "<<tmp_lep->numberOfDaughters()<<endl;
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == PDGID) {
                //        cout<<"Lepton found putting it into swap"<<endl;
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    //const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    rad_photon_energy.push_back(test_particle->energy());
//                    cout<<"Size of rad_photon_energy: "<<rad_photon_energy.size()<<endl;
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
//                cout<<"Swapping the lepton"<<endl;
            }
        }

        //
        return rad_photon_energy;
    }  

//Function to find the FSR photons pt vector
std::vector<double> GetRadPhotonPt(
            const reco::GenParticle * const BORN_LEP,
            const reco::GenParticle * const NAKED_LEP,
            const int PDGID
            ) {
        //double MASS;
        if (fabs(PDGID) == 11){
        //MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        //MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        //MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. GetDressedLepton will return a nullptr"<<endl;
                return std::vector<double>();
        }
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        std::vector<double> rad_photon_pt;
        assert(rad_photon_pt.size() == 0);
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return std::vector<double>();
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == PDGID) {
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    //const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //cout<<"Status of photon: "<<test_particle->status()<<endl;
                    rad_photon_pt.push_back(test_particle->pt());
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
            }
        }

        return rad_photon_pt;
    }

////Function to store the dR values of all photons emitted in the born -> naked process: Using a modified version of the GeDressed Function////
vector<double> GetdR( const reco::GenParticle * const BORN_LEP, const reco::GenParticle * const NAKED_LEP, const int PDGID){
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_ELECTRON->pt(), NAKED_ELECTRON->eta(), NAKED_ELECTRON->phi(), ELECTRON_MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        
        //double MASS;
        if (fabs(PDGID) == 11){
        //MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        //MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        //MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. GetDressedLepton will return a nullptr"<<endl;
                return std::vector<double>();
        }
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        std::vector<double> dR;
        assert(dR.size() == 0);
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return std::vector<double>();
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == PDGID) {
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    //const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //cout<<"Status of photon: "<<test_particle->status()<<endl;
                    double delta = deltaR(tmp_lep->eta(), tmp_lep->phi(), test_particle->eta(), test_particle->phi());
                    dR.push_back(delta);
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
            }
        }

        return dR;       
}

//Function to return the ranking of elements in vector
vector<int> index_list_db(vector<double> v) {
    vector<double> l = v;
    vector<int> index;
    sort(l.begin(), l.end(), greater<double>());  // Sort in descending order
    set<int> used_indices;  // To keep track of used indices

    for (size_t j = 0; j < l.size(); j++) {
        for (size_t i = 0; i < v.size(); i++) {
            // Check if the value matches and the index hasn't been used yet
            if (l[j] == v[i] && used_indices.find(i) == used_indices.end()) {
                index.push_back(i);
                used_indices.insert(i);  // Mark this index as used
                break;
            }
        }
    }

    return index;
}


//Fucntion which states the preselection condition
bool preselection_db(vector<double> pt, vector<double> eta, vector<double> phi, vector<double> energy){
	bool total_condition;
	bool condition1 = energy.size()>=2;
	if (condition1){
		size_t e_lead=index_list_db(pt)[0];
		size_t e_sublead=index_list_db(pt)[1];
		ROOT::Math::PtEtaPhiEVector ele_1((pt)[e_lead], (eta)[e_lead],(phi)[e_lead], (energy)[e_lead]);
		ROOT::Math::PtEtaPhiEVector ele_2((pt)[e_sublead], (eta)[e_sublead],(phi)[e_sublead], (energy)[e_sublead]);
		assert(ele_1.Pt()>=ele_2.Pt());

		bool condition2= ele_1.Pt()>=32;
		bool condition3= ele_2.Pt()>=25;
		bool condition4= fabs(ele_1.Eta())<2.5;
		bool condition5= fabs(ele_2.Eta())<2.5;
		bool condition6= (ele_1+ele_2).M()<110.0;
		bool condition7= (ele_1+ele_2).M()>70.0;
		total_condition=condition1 && condition2 && condition3 && condition4 && condition5 && condition6 && condition7;
	}
	else{
		total_condition= false;
	}
	return total_condition;
}

////////////Function Definition: GetNakedLepton///////////////////////////////////
const reco::GenParticle* GetNakedLepton(const reco::GenParticle* const BORN_LEP, const int PDGID){
        //Initialising the pointer
        const reco::GenParticle* naked_lep=BORN_LEP;
        //Now iterate over the daughters untill find a stable electron
        while (naked_lep ->status()!=1){
                //if it does not have any daughters then return nullptr
                if(naked_lep->numberOfDaughters()==0){
                        return nullptr;
                }
                //if not then we continue our search
                for(size_t i=0; i< naked_lep->numberOfDaughters(); i++){
                        const reco::Candidate* test_part = naked_lep->daughter(i);
                        if (fabs(test_part->pdgId())==fabs(PDGID)){
                                naked_lep=dynamic_cast<const reco::GenParticle*>(test_part);
                                break;
                        }

                }
                
        }
        return naked_lep;

}

//index_list which is supposed to work on floats
bool descending(double i, double j) { return i > j; }

vector<int> index_list(vector<float> v) {
    vector<float> l = v;
    vector<int> index;
    sort(l.begin(), l.end(), greater<float>());  // Sort in descending order
    set<int> used_indices;  // To keep track of used indices

    for (size_t j = 0; j < l.size(); j++) {
        for (size_t i = 0; i < v.size(); i++) {
            // Check if the value matches and the index hasn't been used yet
            if (l[j] == v[i] && used_indices.find(i) == used_indices.end()) {
                index.push_back(i);
                used_indices.insert(i);  // Mark this index as used
                break;
            }
        }
    }

    return index;
}

bool preselection_vector(const std::vector<ROOT::Math::PtEtaPhiEVector>& elec_vec){
        std::vector<double> pt;
        std::vector<double> eta;
        std::vector<double> phi;
        std::vector<double> energy;
        assert(pt.size()==0);
        assert(eta.size()==0);
        assert(phi.size()==0);
        assert(energy.size()==0);
        for(const auto& vec : elec_vec){
                pt.push_back(vec.pt());
                eta.push_back(vec.eta());
                phi.push_back(vec.phi());
                energy.push_back(vec.energy());
        }
	bool total_condition;
	bool condition1 = energy.size()>=2;
	if (condition1){
		size_t e_lead=index_list_db(pt)[0];
		size_t e_sublead=index_list_db(pt)[1];
		ROOT::Math::PtEtaPhiEVector ele_1((pt)[e_lead], (eta)[e_lead],(phi)[e_lead], (energy)[e_lead]);
		ROOT::Math::PtEtaPhiEVector ele_2((pt)[e_sublead], (eta)[e_sublead],(phi)[e_sublead], (energy)[e_sublead]);
		assert(ele_1.Pt()>=ele_2.Pt());

		bool condition2= ele_1.Pt()>=32;
		bool condition3= ele_2.Pt()>=25;
		bool condition4= fabs(ele_1.Eta())<2.4;
		bool condition5= fabs(ele_2.Eta())<2.4;
		bool condition6= (ele_1+ele_2).M()<110.0;
		bool condition7= (ele_1+ele_2).M()>70.0;
		total_condition=condition1 && condition2 && condition3 && condition4 && condition5 && condition6 && condition7;
	}
	else{
		total_condition= false;
	}
	return total_condition;
}

#endif