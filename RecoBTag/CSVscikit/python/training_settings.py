import FWCore.ParameterSet.Config as cms

## IMPORTANT!
## This file was automatically generated by RecoBTag/CTagging/test/dump_training_vars_cfg.py
## with input xml files:
##    - C vs L: ../data/c_vs_udsg.weight.xml sha1 checksum: 1b50773894bf3c64e41694bd48bda5f6f0e3795b
##    - C vs B: ../data/c_vs_b.weight.xml sha1 checksum: c342f54c6448d488e6e2b483a3a3956e34ad8ea1 
csvscikit_vpset = cms.VPSet(
    
		cms.PSet(
		    default = cms.double(-0.1),#this one I am not sure
		    idx = cms.int32(0),
		    name = cms.string('Jet_pt'),
		    taggingVarName = cms.string('jetPt')
		 ),
		cms.PSet(
		    default = cms.double(-5),#this one I am not sure
		    idx = cms.int32(0),
		    name = cms.string('Jet_eta'),
		    taggingVarName = cms.string('jetEta')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    #idx = cms.int32(0),
		    name = cms.string('TagVarCSV_jetNTracks'),
		    taggingVarName = cms.string('jetNTracks')
		 ),
		cms.PSet(
		    default = cms.double(-100),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackSip3dSig_0'),
		    taggingVarName = cms.string('trackSip3dSig')
		 ),
		cms.PSet(
		    default = cms.double(-100),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackSip3dSig_1'),
		    taggingVarName = cms.string('trackSip3dSig')
		 ),
		cms.PSet(
		    default = cms.double(-100),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackSip3dSig_2'),
		    taggingVarName = cms.string('trackSip3dSig')
		 ),
		cms.PSet(
		    default = cms.double(-100),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackSip3dSig_3'),
		    taggingVarName = cms.string('trackSip3dSig')
		 ),
		cms.PSet(
		    default = cms.double(-999),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackSip3dSigAboveCharm'),
		    taggingVarName = cms.string('trackSip3dSigAboveCharm')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackPtRel_0'),
		    taggingVarName = cms.string('trackPtRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackPtRel_1'),
		    taggingVarName = cms.string('trackPtRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackPtRel_2'),
		    taggingVarName = cms.string('trackPtRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackPtRel_3'),
		    taggingVarName = cms.string('trackPtRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackEtaRel_0'),
		    taggingVarName = cms.string('trackEtaRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackEtaRel_1'),
		    taggingVarName = cms.string('trackEtaRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackEtaRel_2'),
		    taggingVarName = cms.string('trackEtaRel')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackEtaRel_3'),
		    taggingVarName = cms.string('trackEtaRel')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackDeltaR_0'),
		    taggingVarName = cms.string('trackDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackDeltaR_1'),
		    taggingVarName = cms.string('trackDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackDeltaR_2'),
		    taggingVarName = cms.string('trackDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackDeltaR_3'),
		    taggingVarName = cms.string('trackDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackPtRatio_0'),
		    taggingVarName = cms.string('trackPtRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackPtRatio_1'),
		    taggingVarName = cms.string('trackPtRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackPtRatio_2'),
		    taggingVarName = cms.string('trackPtRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackPtRatio_3'),
		    taggingVarName = cms.string('trackPtRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackJetDist_0'),
		    taggingVarName = cms.string('trackJetDist')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackJetDist_1'),
		    taggingVarName = cms.string('trackJetDist')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackJetDist_2'),
		    taggingVarName = cms.string('trackJetDist')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackJetDist_3'),
		    taggingVarName = cms.string('trackJetDist')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackDecayLenVal_0'),
		    taggingVarName = cms.string('trackDecayLenVal')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(1),
		    name = cms.string('TagVarCSV_trackDecayLenVal_1'),
		    taggingVarName = cms.string('trackDecayLenVal')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(2),
		    name = cms.string('TagVarCSV_trackDecayLenVal_2'),
		    taggingVarName = cms.string('trackDecayLenVal')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(3),
		    name = cms.string('TagVarCSV_trackDecayLenVal_3'),
		    taggingVarName = cms.string('trackDecayLenVal')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    #idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackSumJetEtRatio'),
		    taggingVarName = cms.string('trackSumJetEtRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    #idx = cms.int32(0),
		    name = cms.string('TagVarCSV_trackSumJetDeltaR'),
		    taggingVarName = cms.string('trackSumJetDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_vertexMass'),
		    taggingVarName = cms.string('vertexMass')
		 ),
		cms.PSet(
		    default = cms.double(0),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_vertexNTracks'),
		    taggingVarName = cms.string('vertexNTracks')
		 ),
		cms.PSet(
		    default = cms.double(-10),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_vertexEnergyRatio'),
		    taggingVarName = cms.string('vertexEnergyRatio')
		 ),
		cms.PSet(
		    default = cms.double(-0.1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_vertexJetDeltaR'),
		    taggingVarName = cms.string('vertexJetDeltaR')
		 ),
		cms.PSet(
		    default = cms.double(-1),
		    idx = cms.int32(0),
		    name = cms.string('TagVarCSV_flightDistance2dSig'),
		    taggingVarName = cms.string('flightDistance2dSig')
		 ),
		cms.PSet(
		    default = cms.double(0),
		    #idx = cms.int32(0),
		    name = cms.string('TagVarCSV_jetNSecondaryVertices'),
		    taggingVarName = cms.string('jetNSecondaryVertices')
		 ),
    
    )
