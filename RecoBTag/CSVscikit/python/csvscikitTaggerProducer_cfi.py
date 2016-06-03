import FWCore.ParameterSet.Config as cms
#use import as to mask it to process.load() 
import RecoBTag.SecondaryVertex.candidateCombinedSecondaryVertexSoftLeptonComputer_cfi as sl_cfg 
from RecoBTag.CSVscikit.training_settings import csvscikit_vpset

#charmTagsComputerCvsL = cms.ESProducer(
CSVscikitTags = cms.ESProducer(
   #'CharmTaggerESProducer',
   'CSVscikitESProducer',
   #clone the cfg only
   slComputerCfg = cms.PSet(
      **sl_cfg.candidateCombinedSecondaryVertexSoftLeptonComputer.parameters_()
      ),
   weightFile = cms.FileInPath('RecoBTag/CSVscikit/data/TMVAClassification_BDTG.weights.xml'),
   variables = csvscikit_vpset,
   computer = cms.ESInputTag('combinedSecondaryVertexSoftLeptonComputer'),
   tagInfos = cms.VInputTag(
      cms.InputTag('pfImpactParameterTagInfos'),
      cms.InputTag('pfInclusiveSecondaryVertexFinderCvsLTagInfos'),
      cms.InputTag('softPFMuonsTagInfos'),
      cms.InputTag('softPFElectronsTagInfos'),
      ),
   mvaName = cms.string('BDT'),
   useCondDB = cms.bool(False),
   gbrForestLabel = cms.string(''),
   useGBRForest = cms.bool(True),
   useAdaBoost = cms.bool(False)
   )

#charmTagsComputerCvsL.slComputerCfg.correctVertexMass = False
CSVscikitTags.slComputerCfg.correctVertexMass = False

#charmTagsComputerCvsB = charmTagsComputerCvsL.clone(
#CSVscikitTagsdummy = CSVscikitTags.clone(
#   weightFile = cms.FileInPath('RecoBTag/CSVscikit/data/TMVAClassification_BDTG.weights.xml'),   
#   variables = c_vs_b_vars_vpset
#   )

