/*
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
//#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
*/

#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>


using namespace reco;


inline edm::ParameterSet CombinedSVComputer::dropDeltaR(const edm::ParameterSet &pset) const
{
	edm::ParameterSet psetCopy(pset);
	psetCopy.addParameter<double>("jetDeltaRMax", 99999.0);
	return psetCopy;
}

CombinedSVComputer::CombinedSVComputer(const edm::ParameterSet &params) :
	trackFlip(params.getParameter<bool>("trackFlip")),
	vertexFlip(params.getParameter<bool>("vertexFlip")),
	charmCut(params.getParameter<double>("charmCut")),
	sortCriterium(TrackSorting::getCriterium(params.getParameter<std::string>("trackSort"))),
	trackSelector(params.getParameter<edm::ParameterSet>("trackSelection")),
	trackNoDeltaRSelector(dropDeltaR(params.getParameter<edm::ParameterSet>("trackSelection"))),
	trackPseudoSelector(params.getParameter<edm::ParameterSet>("trackPseudoSelection")),
	pseudoMultiplicityMin(params.getParameter<unsigned int>("pseudoMultiplicityMin")),
	trackMultiplicityMin(params.getParameter<unsigned int>("trackMultiplicityMin")),
	minTrackWeight(params.getParameter<double>("minimumTrackWeight")),
	useTrackWeights(params.getParameter<bool>("useTrackWeights")),
	vertexMassCorrection(params.getParameter<bool>("correctVertexMass")),
	pseudoVertexV0Filter(params.getParameter<edm::ParameterSet>("pseudoVertexV0Filter")),
	trackPairV0Filter(params.getParameter<edm::ParameterSet>("trackPairV0Filter")),
	
	//For scikit tagger
	//sl_computer_(params.getParameter<edm::ParameterSet>("slComputerCfg")),
        mva_name_( params.getParameter<std::string >("mvaName") ),
        use_condDB_(params.getParameter<bool>("useCondDB")),
        gbrForest_label_(params.getParameter<std::string>("gbrForestLabel")),
        weight_file_(params.getParameter<edm::FileInPath>("weightFile")),
        use_GBRForest_(params.getParameter<bool>("useGBRForest")),
        use_adaBoost_(params.getParameter<bool>("useAdaBoost"))
{
	vpset vars_definition = params.getParameter<vpset>("variables");
	for(auto &var : vars_definition) {
		MVAVar mva_var;
		mva_var.name = var.getParameter<std::string>("name");
		mva_var.id = reco::getTaggingVariableName(
			var.getParameter<std::string>("taggingVarName")
			);
		mva_var.has_index = var.existsAs<int>("idx") ;
		mva_var.index = mva_var.has_index ? var.getParameter<int>("idx") : 0;
		mva_var.default_value = var.getParameter<double>("default");

		variables_.push_back(mva_var);
	}

	uses(0, "pfImpactParameterTagInfos");
	uses(1, "pfInclusiveSecondaryVertexFinderCvsLTagInfos");
	uses(2, "softPFMuonsTagInfos");
	uses(3, "softPFElectronsTagInfos");

}

inline double CombinedSVComputer::flipValue(double value, bool vertex) const
{
	return (vertex ? vertexFlip : trackFlip) ? -value : value;
}

inline CombinedSVComputer::IterationRange CombinedSVComputer::flipIterate(
						int size, bool vertex) const
{
	IterationRange range;
	if (vertex ? vertexFlip : trackFlip) {
		range.begin = size - 1;
		range.end = -1;
		range.increment = -1;
	} else {
		range.begin = 0;
		range.end = size;
		range.increment = +1;
	}

	return range;
}

const btag::TrackIPData &
CombinedSVComputer::threshTrack(const CandIPTagInfo &trackIPTagInfo,
                                const btag::SortCriteria sort,
                                const reco::Jet &jet,
                                const GlobalPoint &pv) const
{
        const CandIPTagInfo::input_container &tracks =
                                        trackIPTagInfo.selectedTracks();
        const std::vector<btag::TrackIPData> &ipData =
                                        trackIPTagInfo.impactParameterData();
        std::vector<std::size_t> indices = trackIPTagInfo.sortedIndexes(sort);

        IterationRange range = flipIterate(indices.size(), false);
        TrackKinematics kin;
        range_for(i, range) {
                std::size_t idx = indices[i];
                const btag::TrackIPData &data = ipData[idx];
                const CandidatePtr &track = tracks[idx];

                if (!trackNoDeltaRSelector(track, data, jet, pv))
                        continue;

                kin.add(track);
                if (kin.vectorSum().M() > charmCut)
                        return data;
        }

        static const btag::TrackIPData dummy = {
                GlobalPoint(),
                GlobalPoint(),
                Measurement1D(-1.0, 1.0),
                Measurement1D(-1.0, 1.0),
                Measurement1D(-1.0, 1.0),
                Measurement1D(-1.0, 1.0),
                0.
        };
        return dummy;
}

const btag::TrackIPData &
CombinedSVComputer::threshTrack(const TrackIPTagInfo &trackIPTagInfo,
                                const btag::SortCriteria sort,
                                const reco::Jet &jet,
                                const GlobalPoint &pv) const
{
	const edm::RefVector<TrackCollection> &tracks =
					trackIPTagInfo.selectedTracks();
	const std::vector<btag::TrackIPData> &ipData =
					trackIPTagInfo.impactParameterData();
	std::vector<std::size_t> indices = trackIPTagInfo.sortedIndexes(sort);

	IterationRange range = flipIterate(indices.size(), false);
	TrackKinematics kin;
	range_for(i, range) {
		std::size_t idx = indices[i];
		const btag::TrackIPData &data = ipData[idx];
		const Track &track = *tracks[idx];

		if (!trackNoDeltaRSelector(track, data, jet, pv))
			continue;

		kin.add(track);
		if (kin.vectorSum().M() > charmCut) 
			return data;
	}

	static const btag::TrackIPData dummy = {
 		GlobalPoint(),
		GlobalPoint(),
		Measurement1D(-1.0, 1.0),
		Measurement1D(-1.0, 1.0),
		Measurement1D(-1.0, 1.0),
		Measurement1D(-1.0, 1.0),
		0.
	};
	return dummy;
}


TaggingVariableList
CombinedSVComputer::operator () (const TrackIPTagInfo &ipInfo,
                                 const SecondaryVertexTagInfo &svInfo) const
{
	using namespace ROOT::Math;

	edm::RefToBase<Jet> jet = ipInfo.jet();
	math::XYZVector jetDir = jet->momentum().Unit();
	TaggingVariableList vars;

        TrackKinematics vertexKinematics;
	
	double vtx_track_ptSum = 0.;
	double vtx_track_ESum  = 0.;
	
	// the following is specific depending on the type of vertex
        int vtx = -1;
        unsigned int numberofvertextracks = 0;

	IterationRange range = flipIterate(svInfo.nVertices(), true);
	range_for(i, range) {

		numberofvertextracks = numberofvertextracks + (svInfo.secondaryVertex(i)).nTracks();

		const Vertex &vertex = svInfo.secondaryVertex(i);
		bool hasRefittedTracks = vertex.hasRefittedTracks();
		for(reco::Vertex::trackRef_iterator track = vertex.tracks_begin(); track != vertex.tracks_end(); track++) {
			double w = vertex.trackWeight(*track);
			if (w < minTrackWeight)
				continue;
			if (hasRefittedTracks) {
				const Track actualTrack = vertex.refittedTrack(*track);
				vertexKinematics.add(actualTrack, w);
				vars.insert(btau::trackEtaRel, reco::btau::etaRel(jetDir,actualTrack.momentum()), true);
				if(vtx < 0) // calculate this only for the first vertex
				{
					vtx_track_ptSum += std::sqrt(actualTrack.momentum().Perp2());
					vtx_track_ESum  += std::sqrt(actualTrack.momentum().Mag2() + ROOT::Math::Square(ParticleMasses::piPlus));
				}
			} else {
				vertexKinematics.add(**track, w);
				vars.insert(btau::trackEtaRel, reco::btau::etaRel(jetDir,(*track)->momentum()), true);
				if(vtx < 0) // calculate this only for the first vertex
				{
					vtx_track_ptSum += std::sqrt((*track)->momentum().Perp2());
					vtx_track_ESum  += std::sqrt((*track)->momentum().Mag2() + ROOT::Math::Square(ParticleMasses::piPlus));
				}
			}
		}
		
		if (vtx < 0) vtx = i;
        }
	if(vtx>=0){
		vars.insert(btau::vertexNTracks, numberofvertextracks, true);
		vars.insert(btau::vertexFitProb,(svInfo.secondaryVertex(vtx)).normalizedChi2(), true);
	}

	// after we collected vertex information we let the common code complete the job
	fillCommonVariables(vars,vertexKinematics,ipInfo,svInfo,vtx_track_ptSum,vtx_track_ESum);

	vars.finalize();
	return vars;
}

TaggingVariableList
CombinedSVComputer::operator () (const CandIPTagInfo &ipInfo,
                                 const CandSecondaryVertexTagInfo &svInfo) const
{
        using namespace ROOT::Math;

        edm::RefToBase<Jet> jet = ipInfo.jet();
        math::XYZVector jetDir = jet->momentum().Unit();
        TaggingVariableList vars;

        TrackKinematics vertexKinematics;
	
	double vtx_track_ptSum = 0.;
	double vtx_track_ESum  = 0.;
	
	// the following is specific depending on the type of vertex
        int vtx = -1;
        unsigned int numberofvertextracks = 0;

	IterationRange range = flipIterate(svInfo.nVertices(), true);
	range_for(i, range) {

		numberofvertextracks = numberofvertextracks + (svInfo.secondaryVertex(i)).numberOfSourceCandidatePtrs();

		const reco::VertexCompositePtrCandidate &vertex = svInfo.secondaryVertex(i);
		const std::vector<CandidatePtr> & tracks = vertex.daughterPtrVector();
		for(std::vector<CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
			vertexKinematics.add(*track);
			vars.insert(btau::trackEtaRel, reco::btau::etaRel(jetDir,(*track)->momentum()), true);
			if(vtx < 0) // calculate this only for the first vertex
			{
				vtx_track_ptSum += std::sqrt((*track)->momentum().Perp2());
				vtx_track_ESum  += std::sqrt((*track)->momentum().Mag2() + ROOT::Math::Square(ParticleMasses::piPlus));
			}
		}
		
		if (vtx < 0) vtx = i;
	}
	if(vtx>=0){
		vars.insert(btau::vertexNTracks, numberofvertextracks, true);
		vars.insert(btau::vertexFitProb,(svInfo.secondaryVertex(vtx)).vertexNormalizedChi2(), true);
	}
	
	// after we collected vertex information we let the common code complete the job
	fillCommonVariables(vars,vertexKinematics,ipInfo,svInfo,vtx_track_ptSum,vtx_track_ESum);
	
	vars.finalize();
	return vars;
}

//For scikit tagger
void CombinedSVComputer::initialize(const JetTagComputerRecord & record)
{
	mvaID_.reset(new TMVAEvaluator());

	std::vector<std::string> variable_names;
	variable_names.reserve(variables_.size());

	for(auto &var : variables_) {
		variable_names.push_back(var.name);
	}
	std::vector<std::string> spectators;

  if(use_condDB_) {
		const GBRWrapperRcd & gbrWrapperRecord = record.getRecord<GBRWrapperRcd>();

		edm::ESHandle<GBRForest> gbrForestHandle;
		gbrWrapperRecord.get(gbrForest_label_.c_str(), gbrForestHandle);

		mvaID_->initializeGBRForest(
			gbrForestHandle.product(), variable_names, 
			spectators, use_adaBoost_
			);
  }
  else {
    mvaID_->initialize(
			"Color:Silent:Error", mva_name_.c_str(),
			weight_file_.fullPath(), variable_names, 
			spectators, use_GBRForest_, use_adaBoost_
			);
  }
}

/// b-tag a jet based on track-to-jet parameters in the extened info collection
float CombinedSVComputer::discriminator(const TagInfoHelper & tagInfo) const {
  // default value, used if there are no leptons associated to this jet
  const reco::CandIPTagInfo & ip_info = tagInfo.get<reco::CandIPTagInfo>(0);
	const reco::CandSecondaryVertexTagInfo & sv_info = tagInfo.get<reco::CandSecondaryVertexTagInfo>(1);
	//const reco::CandSoftLeptonTagInfo& softmu_info = tagInfo.get<reco::CandSoftLeptonTagInfo>(2);
	//const reco::CandSoftLeptonTagInfo& softel_info = tagInfo.get<reco::CandSoftLeptonTagInfo>(3);
	//reco::TaggingVariableList vars = sl_computer_(ip_info, sv_info, softmu_info, softel_info);
	reco::TaggingVariableList vars = CombinedSVComputer::operator()(ip_info, sv_info);

	// Loop over input variables
	std::map<std::string, float> inputs;
	for(auto &mva_var : variables_){
		//vectorial tagging variable
		if(mva_var.has_index){
			std::vector<float> vals = vars.getList(mva_var.id, false);
			inputs[mva_var.name] = (vals.size() > mva_var.index) ? vals[mva_var.index] : mva_var.default_value;
		}
		//single value tagging var
		else {
			inputs[mva_var.name] = vars.get(mva_var.id, mva_var.default_value);
		}
	}

	//get the MVA output
	float tag = mvaID_->evaluate(inputs);
	return tag;
}

