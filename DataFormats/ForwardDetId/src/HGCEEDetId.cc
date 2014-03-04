#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <ostream>
#include <iostream>

const HGCEEDetId HGCEEDetId::Undefined(ForwardEmpty,0,0,0,0,0);

HGCEEDetId::HGCEEDetId() : DetId() {
}

HGCEEDetId::HGCEEDetId(uint32_t rawid) : DetId(rawid) {
}

HGCEEDetId::HGCEEDetId(ForwardSubdetector subdet, int zp, int lay, int mod, int subsec, int cell) : DetId(Forward,subdet) 
{  
  uint32_t rawid=0;
  rawid |= ((cell   & 0xffff) << 0 );
  rawid |= ((mod    & 0x1f)   << 16);
  rawid |= ((subsec & 0x1)    << 21);
  rawid |= ((lay    & 0x1f)   << 22);
  if(zp>0) rawid |= ((zp     & 0x1)    << 27);
  rawid |= (1                 << 28);
  rawid |= ((Forward & 0x7)   << 29);
  id_=rawid;
}

HGCEEDetId::HGCEEDetId(const DetId& gen) {
  if (!gen.null()) {
    ForwardSubdetector subdet=(ForwardSubdetector(gen.subdetId()));
    if (gen.det()!=Forward || (subdet!=HGCEE)) {
      throw cms::Exception("Invalid DetId") << "Cannot initialize HGCEEDetId from " << std::hex << gen.rawId() << std::dec; 
    }  
  }
  id_ = gen.rawId();
}

HGCEEDetId& HGCEEDetId::operator=(const DetId& gen) {
  if (!gen.null()) {
    ForwardSubdetector subdet=(ForwardSubdetector(gen.subdetId()));
    if (gen.det()!=Forward || (subdet!=HGCEE)) {
      throw cms::Exception("Invalid DetId") << "Cannot assign HGCEEDetId from " << std::hex << gen.rawId() << std::dec; 
    }  
  }
  id_ = gen.rawId();
  return (*this);
}

std::ostream& operator<<(std::ostream& s,const HGCEEDetId& id) {
  switch (id.subdet()) {
  case(HGCEE) : return s << "isEE=" << id.isEE() 
			 << " zpos=" << id.zside() 
			 << " layer=" << id.layer() 
			 << " phi sub-sector" << id.subsector()
			 << " module=" << id.module() 
			 << " cell=" << id.cell();
  default : return s << id.rawId();
  }
}


