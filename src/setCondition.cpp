#include "md.hpp"

void
MD::setCondition(char* condfile, char* file){
	readCondFile(condfile, file);
	pp->readAtomsFile(file);
	pp->setPhysicalProp(gastype,vaportype,T,p);
}

void
MD::readCondFile(char* condfile, char* file){
	ifstream stream(condfile);
	string str;
	int iflag=0;
	cout<<"**************************************************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"************  Calculation conditions  ************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"**************************************************"<<endl;
	while(getline(stream,str)) {
		cout<<iflag<<endl;
		if(str.length()==0) continue;
		if (str=="Gas") {iflag=1; continue;}
		if (str=="Calculation number") {iflag=2; continue;}
		if (str=="Temperature"||str=="T") {iflag=3; continue;}
		if (str=="Pressure"||str=="P") {iflag=4; continue;}
		if (str=="Time step"||str=="dt") {iflag=5; continue;}
		if (str=="Total number of steps") {iflag=6; continue;}
		if (str=="Number of steps for relax") {iflag=7; continue;}
		if (str=="Vapor") {iflag=8; continue;}
		if (str=="Cut off length") {iflag=9; continue;}
		if (str=="Margin size") {iflag=10; continue;}
		if (str=="Gyration path") {iflag=16; continue;}
		if (str=="RDF path") {iflag=17; continue;}
		if (str=="gas gas interaction") {iflag=18; continue;}
		if (str=="gas ion interaction") {iflag=19; continue;}
		if (str=="ion ion interaction") {iflag=20; continue;}
		if (str=="Nose-Hoover ion") {iflag=21; continue;}
		if (str=="Nose-Hoover gas") {iflag=22; continue;}
		if (str=="Output") {iflag=23; continue;}
		if (str=="Electric field") {iflag=25; continue;}
		if (str=="Input") {iflag=26; continue;}
		if (iflag==1) {
			istringstream stream(str);
			string tmp;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0){
					if(tmp=="He"){
						gastype=1;
						cout<<"Gastype\t\t\tHe"<<endl;
					}
					if(tmp=="N2"){
						gastype=2;
						cout<<"Gastype\t\t\tN2(diatomic)"<<endl;
					}
					if(tmp=="N2monoatomic"){
						gastype=3;
						cout<<"Gastype\t\t\tN2(monoatomic)"<<endl;
					}
					if(tmp=="Ar"){
						gastype=4;
						cout<<"Gastype\t\t\tAr"<<endl;
					}
				}
				if (loop==1) {
					Nof_around_gas=stoi(tmp);
					cout<<"Number of gases\t\t"<<Nof_around_gas<<endl;
				}
				loop++;
			}
		}
		if (iflag==2) {
			calculation_number=stoi(str);
		}
		if (iflag==3) {
			T=stod(str);
			cout<<"Temperature\t\t"<<T<<" K"<<endl;
		}
		if (iflag==4) {
			p=stod(str);
			cout<<"Pressure\t\t"<<p<<" Pa"<<endl;
		}
		if (iflag==5) {
			dt=stod(str);
			cout<<"Time step\t\t"<<dt<<" fs"<<endl;
		}
		if (iflag==6) {
			double inter=stod(str);
			Noftimestep=inter;
			cout<<"Total steps\t\t"<<float(Noftimestep)<<endl;
		}
		if (iflag==7) {
			double inter=stod(str);
			step_relax=inter;
			cout<<"Relax steps\t\t"<<float(step_relax)<<endl;
		}
		if (iflag==8) {
			istringstream stream(str);
			string tmp;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0){
					if(tmp=="MeOH"){
						vaportype=1;
						cout<<"Vapor type\t\tMeOH"<<endl;
					}
					if(tmp=="H2O"){
						vaportype=2;
						cout<<"Vapor type\t\tH2O"<<endl;
					}
					if(tmp=="EtOH"){
						vaportype=3;
						cout<<"Vapor type\t\tEtOH"<<endl;
					}
				}
				if (loop==1) {
					Nof_around_vapor=stoi(tmp);
					cout<<"Number of vapors\t"<<Nof_around_vapor<<endl;
				}
				loop++;
			}
		}
		if (iflag==9) {
			CUTOFF=stod(str);
			cout<<"Cutoff\t\t\t"<<CUTOFF<<" ang."<<endl;
		}
		if (iflag==10) {
			MARGIN=stod(str);
			cout<<"Margin size\t\t"<<MARGIN<<" ang."<<endl;
		}
		if (iflag==16) {
			ostringstream ss;
			ss<<str<<"_"<<calculation_number<<".dat";
			string tmp=ss.str();
			pp->gyration_path=new char[tmp.length()+10];
			strcpy(pp->gyration_path,tmp.c_str());
			flags->gyration=1;
			FILE*f=fopen(pp->gyration_path, "w");
			fclose(f);
			cout<<"Gyration --> ON -->\t"<<tmp<<endl;
		}
		if (iflag==17) {
			pp->RDF_path=new char[str.length()+1];
			strcpy(pp->RDF_path,str.c_str());
			flags->RDF=1;
			cout<<"RDF\tON\tPath="<<str<<endl;
		}
		if (iflag==18) {
			if (str=="LJ") {
				flags->inter_gg=1;
				cout<<"gas-gas interaction -->\tON"<<endl;
			}
			else if (str=="LJsemiNVT") {
				flags->inter_gg=1;
				flags->semi_NVT_gasgas=1;
				cout<<"gas-gas interaction -->\tON"<<endl;
				cout<<"Boundary rescaling -->\tON"<<endl;
			}
			else if (str=="OFF") {
				flags->inter_gg=0;
				flags->semi_NVT_gasgas=0;
				cout<<"gas-gas interaction -->\tOFF"<<endl;
				cout<<"Boundary rescaling -->\tOFF"<<endl;
			}
			else printf("**************Uknown gas gas parameter was found**************\n");
		}
		if (iflag==19) {
			if (str=="LJ") flags->force_lj=1;
			else if (str=="ion dipole") flags->force_ion_dipole=1;
			else if (str=="OFF") {
				flags->force_ion_dipole=0;
				flags->force_lj=0;
			}
			else printf("**************Uknown gas ion parameter was found**************\n");
		}
		if (iflag==20) {
			if (str=="AMBER") flags->intra_AMBER=1;
			else if (str=="Stilinger-Weber") flags->force_sw=1;
			else if (str=="Tersoff") flags->force_ters=1;
			else if (str=="Born-Mayer-Huggins-NaCl") flags->force_born=1;
			else printf("**************Uknown ion ion parameter was found**************\n");
		}
		if (iflag==21) {
			if (str=="OFF") {flags->nose_hoover_ion=0; cout<<"Nose-Hoover for ion --> OFF"<<endl;}
            else{
                if(str=="scale") {flags->nose_hoover_ion=0;flags->velocity_scaling=1; cout<<"Nose-Hoover for ion --> OFF\nVelocity scaling for ion --> ON"<<endl;}
                else {
                    flags->nose_hoover_ion=1;
                    pp->Tnh_ion=stod(str);
                    cout<<"Nose-Hoover for ion --> ON --> "<<pp->Tnh_ion<<" K"<<endl;
                }
            }
		}
		if (iflag==22) {
			if (str=="OFF") {flags->nose_hoover_gas=0; cout<<"Nose-Hoover for gas --> OFF"<<endl;}
			else {
				flags->nose_hoover_gas=1;
				pp->Tnh_gas=stod(str);
				cout<<"Nose-Hoover for gas --> ON --> "<<pp->Tnh_gas<<" K"<<endl;
			}
		}
		if (iflag==23) {
			istringstream stream(str);
			string tmp;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) {
					ostringstream ss;
					ss<<tmp<<"_"<<calculation_number<<".dump";
					string tmp2=ss.str();
					pp->dump_path=new char[tmp2.length()+1];
					strcpy(pp->dump_path,tmp2.c_str());
					cout<<"Dump file -->\t\t"<<tmp2<<endl;
				}
				if (loop==1) {
					OBSERVE=stoi(tmp);
				}
				if (tmp=="w") {
					FILE*f=fopen(pp->dump_path, "w");
					fclose(f);
				}
				loop++;
			}

		}
		if (iflag==25) {
			istringstream stream(str);
			string tmp;
			int loop=0;
			double mass;
			flags->efield=1;
			while(getline(stream,tmp,'\t')) {
				Ecoeff[loop]=stod(tmp);
				loop++;
			}
		}
		if (iflag==26) {
			strcpy(file,str.c_str());
			cout<<"Atom file -->\t\t"<<str<<endl;
		}
	}
	d_size=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;//pow(28.0855*8/6.02e23/(2.218e-24),1/3.0)*5;
	V=d_size*d_size*d_size;
	CL2 = (CUTOFF)*(CUTOFF);
	ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
	cout<<"Cut off length\t\t"<<CUTOFF<<" ang."<<endl;
	cout<<"Margin length\t\t"<<MARGIN<<" ang."<<endl;
	cout<<"Domain size\t\t"<<d_size<<" ang."<<endl;
}
