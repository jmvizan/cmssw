#include "RecoBTag/CSVscikit/interface/CSVscikitTagger.h"

#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"


//#include <ext/functional>
//#include "DataFormats/BTauReco/interface/TaggingVariable.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>


CSVscikitTagger::CSVscikitTagger(const edm::ParameterSet & configuration):
	sl_computer_(configuration.getParameter<edm::ParameterSet>("slComputerCfg")),
	sv_computer_(configuration.getParameter<edm::ParameterSet>("slComputerCfg")),
  mva_name_( configuration.getParameter<std::string >("mvaName") ),
  use_condDB_(configuration.getParameter<bool>("useCondDB")),
  gbrForest_label_(configuration.getParameter<std::string>("gbrForestLabel")),
  weight_file_(configuration.getParameter<edm::FileInPath>("weightFile")),
  use_GBRForest_(configuration.getParameter<bool>("useGBRForest")),
  use_adaBoost_(configuration.getParameter<bool>("useAdaBoost"))
{
	vpset vars_definition = configuration.getParameter<vpset>("variables");

	std::cout << "CSVscikit::CSVscikit" << std::endl;
	for(auto &var : vars_definition) {
		MVAVar mva_var;
		mva_var.name = var.getParameter<std::string>("name");
		mva_var.id = reco::getTaggingVariableName(
			var.getParameter<std::string>("taggingVarName")
			);
		//if (mva_var.id == btau::lastTaggingVariable) {
		//if (mva_var.id == 1) {
		
		std::cout << "CSVscikit:mva_var.id =" << mva_var.id << " var.getParameter(taggingVarName)=" <<  var.getParameter<std::string>("taggingVarName") << std::endl;
		//}
		
		
		
		mva_var.has_index = var.existsAs<int>("idx") ;
		mva_var.index = mva_var.has_index ? var.getParameter<int>("idx") : 0;
		mva_var.default_value = var.getParameter<double>("default");

		variables_.push_back(mva_var);
	}

	uses(0, "pfImpactParameterTagInfos");
	//uses(1, "pfInclusiveSecondaryVertexFinderCvsLTagInfos");
	uses(1, "pfInclusiveSecondaryVertexFinderTagInfos");
	//uses(2, "softPFMuonsTagInfos");
	//uses(3, "softPFElectronsTagInfos");
}

void CSVscikitTagger::initialize(const JetTagComputerRecord & record)
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

CSVscikitTagger::~CSVscikitTagger()
{
}

/// b-tag a jet based on track-to-jet parameters in the extened info collection
float CSVscikitTagger::discriminator(const TagInfoHelper & tagInfo) const {
  // default value, used if there are no leptons associated to this jet
  const reco::CandIPTagInfo & ip_info = tagInfo.get<reco::CandIPTagInfo>(0);
	const reco::CandSecondaryVertexTagInfo & sv_info = tagInfo.get<reco::CandSecondaryVertexTagInfo>(1);
	//const reco::CandSoftLeptonTagInfo& softmu_info = tagInfo.get<reco::CandSoftLeptonTagInfo>(2);
	//const reco::CandSoftLeptonTagInfo& softel_info = tagInfo.get<reco::CandSoftLeptonTagInfo>(3);
	//reco::TaggingVariableList vars = sl_computer_(ip_info, sv_info, softmu_info, softel_info);
	reco::TaggingVariableList vars = sv_computer_(ip_info, sv_info);

	// Loop over input variables
	std::map<std::string, float> inputs;
	//std::cout << "jetvars->";
	for(auto &mva_var : variables_){
		//vectorial tagging variable
		if(mva_var.has_index){
			std::vector<float> vals = vars.getList(mva_var.id, false);
			inputs[mva_var.name] = (vals.size() > mva_var.index) ? vals[mva_var.index] : mva_var.default_value;
			////std::cout << " var[" << mva_var.name << "]=" << inputs[mva_var.name];
			//std::cout << inputs[mva_var.name] << "\t";
		}
		//single value tagging var
		else {
			inputs[mva_var.name] = vars.get(mva_var.id, mva_var.default_value);
			////std::cout << " var[" << mva_var.name << "]=" << inputs[mva_var.name];
			//std::cout << inputs[mva_var.name] << "\t";
		}
	}

	//get the MVA output
	float tag = mvaID_->evaluate(inputs);
	////std::cout << "---> tag = " << tag << std::endl;
	//std::cout << tag <<"\n";
	return tag;
}
