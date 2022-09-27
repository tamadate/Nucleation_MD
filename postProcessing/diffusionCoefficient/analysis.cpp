#include "analysis.hpp"


int main(int argc,char *argv[]){
	int calculationNumber=stoi(argv[1]);	// calculation number
	int pl=stoi(argv[2]);									// number of replica MD simulations
	double dt_block=stod(argv[6]);							// dt in the block file

	double skipTime=stod(argv[3]);								// skip steps
	double tEND=stod(argv[8]);							// dt in the block file
	int skip=int(skipTime/dt_block);				//
	int stop=int(tEND/dt_block);

	double startTime_block=stod(argv[4]);				// t_inf
	double endTime_block=stod(argv[5]);					// time of a block
	int blockStart=int(startTime_block/dt_block);		//
	int blockSize=int(endTime_block/dt_block);				//

	for(int j=0;j<blockSize;j++){
		MSDx.push_back(0);
		MSDy.push_back(0);
		MSDz.push_back(0);
		MSD.push_back(0);
		VAF.push_back(0);
	}

	int Nblocks=0;
	char dump[100];
	for(int i=0;i<pl;i++){
		sprintf(dump, "%s%s%d%s", argv[7], "ion_300_", i+calculationNumber, ".dat"); // time x y z vx vy vz
		ifstream stream(dump);
		string str;
		int lineIndex=0;
		std::vector<data> datas;
		// loop of ion file
		while(getline(stream,str)) {
			istringstream stream(str);
			string tmp;
			if(lineIndex > skip && lineIndex <= stop){
				data d;
				while(getline(stream,tmp,' ')) d.MD.push_back(stod(tmp));
				datas.push_back(d);
			}
			lineIndex++;
		}
		int N=datas.size()-blockSize-1;
		Nblocks+=N;

		for(int k=0;k<N;k++){
			for(int j=0;j<blockSize;j++){
				int testIndex=k+j;
				double dx=datas[testIndex].MD[1]-datas[k].MD[1];
				double dy=datas[testIndex].MD[2]-datas[k].MD[2];
				double dz=datas[testIndex].MD[3]-datas[k].MD[3];
				MSDx[j]+=dx*dx;
				MSDy[j]+=dy*dy;
				MSDz[j]+=dz*dz;
				MSD[j]+=dx*dx+dy*dy+dz*dz;

				double vxvx=datas[testIndex].MD[4]*datas[k].MD[4];
				double vyvy=datas[testIndex].MD[5]*datas[k].MD[5];
				double vzvz=datas[testIndex].MD[6]*datas[k].MD[6];
				VAF[j]+=vxvx+vyvy+vzvz;
				/*velocityAve[0][j]+=dx;
				velocityAve[1][j]+=dy;
				velocityAve[2][j]+=dz;*/
			}
		}
	}

	for(int i=0;i<blockSize;i++)	{
		VAF[i]/=Nblocks;
		MSD[i]/=Nblocks;
		MSDx[i]/=Nblocks;
		MSDy[i]/=Nblocks;
		MSDz[i]/=Nblocks;
	}

	sprintf(dump, "%s%s%d", argv[7], "TIME_MSD_VAF.", calculationNumber);
	FILE*f=fopen(dump, "w");
	fclose(f);
	for(int i=0;i<blockSize;i++){
		f=fopen(dump, "a");
		fprintf(f, "%f\t%f\t%e\n", dt_block*i*1e9,MSD[i],VAF[i]);
		fclose(f);
	}

	double Dvaf=dt_block/9.0*integral(VAF,blockSize)*1e14;
	double Dmsd=slope(MSD,blockStart,blockSize-1,dt_block)/6.0*1e-20*1e4;
	double Dmsdx=slope(MSDx,blockStart,blockSize-1,dt_block)/2.0*1e-20*1e4;
	double Dmsdy=slope(MSDy,blockStart,blockSize-1,dt_block)/2.0*1e-20*1e4;
	double Dmsdz=slope(MSDz,blockStart,blockSize-1,dt_block)/2.0*1e-20*1e4;
	double coeff=e/kb/T;
	double Kvaf=Dvaf*coeff;
	double Kmsd=Dmsd*coeff;
	double Kmsdx=Dmsdx*coeff;
	double Kmsdy=Dmsdy*coeff;
	double Kmsdz=Dmsdz*coeff;

	sprintf(dump, "%s%s%d", argv[7], "DiffusionCoefficients.", calculationNumber);
	f=fopen(dump, "w");
	fprintf(f, "Dmsdx\tDmdsy\tDmsdz\tDmsd\tDvaf\te/kb/T\tCCScoeff\n");
	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%e\n", Dmsdx,Dmsdy,Dmsdz,Dmsd,Dvaf,coeff,3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)*1e4*1e20);
	fclose(f);
}
