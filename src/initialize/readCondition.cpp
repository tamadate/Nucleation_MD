#include "../md.hpp"


void
MD::readCondFile(char* condfile){
	ifstream stream(condfile);
	string str;
	int iflag=0;
	cout<<"**************************************************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"************  Calculation conditions  ************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"**************************************************"<<endl;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		std::vector<string> readings;
		istringstream stream(str);
		string reading;
		while(getline(stream,reading,' ')) {
			readings.push_back(reading);
		}
		if(readings[0]=="Input"){
			strcpy(atomFile,readings[1].c_str());
			cout<<"Atom file -->\t\t"<<readings[1]<<endl;
		}
		if(readings[0]=="Vapor"){
			strcpy(vaporFile,readings[1].c_str());
			con->Nof_around_vapor=stoi(readings[2]);
			cout<<"Vapor file -->\t\t"<<readings[1]<<endl;
			cout<<"Number of vapors\t"<<con->Nof_around_vapor<<endl;
		}
		if(readings[0]=="VaporStickPositions"){
			strcpy(vaporStickFile,readings[1].c_str());
			con->positionLogStep=stoi(readings[2]);
			ifstream stream2(vaporStickFile);
			string str2;
			while(getline(stream2,str2)) {
				if(str2.length()==0) continue;
				stickPositionList.push_back(stoi(str2)-1);
			}
			cout<<"Vapor stick position file -->\t\t"<<readings[1]<<endl;
		}
		if(readings[0]=="Gas"){
			Atom_type a;
			if(readings[1]=="He"){
				gastype=1;
				a.mass=4.027;
				a.name="He";
				a.coeff1=0.0203;
				a.coeff2=2.556;
			}
			if(readings[1]=="N2"){
				gastype=2;
				a.mass=14.01;
				a.name="N";
				a.coeff1=0.1098;
				a.coeff2=3.27351824993;
				IntraInter.push_back(new PotentialGasIntra());
			}
			if(readings[1]=="N2monoatomic"){
				gastype=3;
				a.mass=28.02;
				a.name="N2";
				a.coeff1=0.14397;
				a.coeff2=3.798;
			}
			if(readings[1]=="Ar"){
				gastype=4;
				a.mass=28.02;
				a.name="Ar";
				a.coeff1=0.14397;
				a.coeff2=3.798;
			}
			con->Nof_around_gas=stoi(readings[2]);
			vars->atypes.push_back(a);
			cout<<"Gastype\t\t\t"<<readings[1]<<endl;
			cout<<"Number of gases\t\t"<<con->Nof_around_gas<<endl;
		}
		if(readings[0]=="Temperature"){
			pp->T=stod(readings[1]);
			cout<<"Temperature\t\t"<<pp->T<<" K"<<endl;
		}
		if(readings[0]=="Pressure"){
			pp->p=stod(readings[1]);
			cout<<"Pressure\t\t"<<pp->p<<" Pa"<<endl;
		}
		if (readings[0]=="dt") {
			dt=stod(readings[1]);
			cout<<"Time step\t\t"<<dt<<" fs"<<endl;
		}
		if (readings[0]=="TotalSteps") {
			con->Noftimestep=stod(readings[1]);
			cout<<"Total steps\t\t"<<float(con->Noftimestep)<<endl;
		}
		if (readings[0]=="RelaxSteps") {
			con->step_relax=stod(readings[1]);
			cout<<"Relax steps\t\t"<<float(con->step_relax)<<endl;
		}
		if (readings[0]=="CutOff") {
			con->CUTOFF=stod(readings[1]);
			cout<<"Cutoff\t\t\t"<<con->CUTOFF<<" ang."<<endl;
		}
		if (readings[0]=="Margin") {
			con->MARGIN=stod(readings[1]);
			cout<<"Margin size\t\t"<<con->MARGIN<<" ang."<<endl;
		}
		if (readings[0]=="Output") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dump";
			string tmp2=ss.str();
			obs->fileDump=new char[tmp2.length()+1];
			strcpy(obs->fileDump,tmp2.c_str());
			obs->OBSERVE=stoi(readings[2]);
			FILE*f=fopen(obs->fileDump, "w");
			fclose(f);
			cout<<"Dump file -->\t\t"<<tmp2<<endl;
		}
		if (readings[0]=="NVTion") {
			if (readings[1]=="OFF") {
				cout<<"Nose-Hoover for ion --> OFF"<<endl;
			}
			else if(readings[1]=="scale") {
				delete thermo;
				thermo = new ThermostatVscale(vars, obs, stod(readings[1]));
				cout<<"Velocity scaling for ion --> ON --> "<<readings[1]<<" K"<<endl;
			}
			else {
				delete thermo;
				thermo = new ThermostatNH(vars, obs, stod(readings[1]));
				cout<<"Nose-Hoover for ion --> ON --> "<<readings[1]<<" K"<<endl;
			}
		}

		if (readings[0]=="Interactions") {
			con->CL2 = con->CUTOFF*con->CUTOFF;
			con->ML2 = (con->CUTOFF+con->MARGIN)*(con->CUTOFF+con->MARGIN);
			continue;
		}
		if (readings[1]=="gg") {
			if (readings[2]=="LJ") InterInter.push_back(new PotentialGasGas(con->ML2));
			else if (readings[2]=="OFF") ;
			else printf("**************Uknown gas gas parameter was found**************\n");
		}
		if (readings[1]=="gi"||readings[1]=="ig") {
			if (readings[2]=="LJ") InterInter.push_back(new PotentialGasIon(con->ML2));
			else if (readings[2]=="ion dipole") InterInter.push_back(new PotentialIonDipole());
			else if (readings[2]=="OFF") ;
			else printf("**************Uknown gas ion parameter was found**************\n");
		}
		if (readings[1]=="ion") {
			if (readings[2]=="AMBER") IntraInter.push_back(new PotentialAMBER());
			else if (readings[2]=="Stilinger-Weber") IntraInter.push_back(new PotentialSW());
			else if (readings[2]=="Tersoff") IntraInter.push_back(new PotentialTersoff());
			else if (readings[2]=="Born-Mayer-Huggins-NaCl") IntraInter.push_back(new PotentialBorn());
			else printf("**************Uknown ion parameter was found**************\n");
		}
		if (readings[1]=="vi"||readings[1]=="iv") {
			if (readings[2]=="LJcoul") InterInter.push_back(new PotentialVaporIon(con->ML2));
			else if (readings[2]=="OFF") ;
			else printf("**************Uknown vapor ion parameter was found**************\n");
		}
		if (readings[1]=="vv"||readings[1]=="vv") {
			if (readings[2]=="LJcoul") InterInter.push_back(new PotentialVaporVapor(con->ML2));
			else if (readings[2]=="OFF") ;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		if (readings[1]=="gv"||readings[1]=="vg") {
			if (readings[2]=="LJ") InterInter.push_back(new PotentialVaporGas(con->ML2));
			else if (readings[2]=="OFF") ;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		if (readings[1]=="Efield") InterInter.push_back(new PotentialEfield(stod(readings[2]),stod(readings[3]),stod(readings[4])));
		if (readings[0]=="Gyration") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dat";
			string tmp=ss.str();
			pp->gyration_path=new char[tmp.length()+10];
			strcpy(pp->gyration_path,tmp.c_str());
			flags->gyration=1;
			FILE*f=fopen(pp->gyration_path, "w");
			fclose(f);
			cout<<"Gyration --> ON -->\t"<<tmp<<endl;
		}
	}
	con->L=pow(con->Nof_around_gas*kb*pp->T/pp->p,1/3.0)*1e10;//pow(28.0855*8/6.02e23/(2.218e-24),1/3.0)*5;
	con->HL=con->L*0.5;
	con->V=con->L*con->L*con->L;
	IntraInter.push_back(new PotentialVaporIntra());
	cout<<"Cut off length\t\t"<<con->CUTOFF<<" ang."<<endl;
	cout<<"Margin length\t\t"<<con->MARGIN<<" ang."<<endl;
	cout<<"Domain size\t\t"<<con->L<<" ang."<<endl;
}
