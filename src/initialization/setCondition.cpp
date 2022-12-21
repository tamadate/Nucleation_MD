#include "../md.hpp"

void
MD::setCondition(char* condfile){
	readCondFile(condfile);
	pp->readIonProp(atomFile);
	pp->readVaporProp(vaporFile);
	pp->setPhysicalProp(gastype,T,p);
}
