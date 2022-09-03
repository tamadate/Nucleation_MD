#include "analysis.hpp"


int main(int argc,char *argv[]){
	int calculationNumber=stoi(argv[1]);	// calculation number
	int pl=stoi(argv[2]);									// number of replica MD simulations
	double skipTime=stod(argv[3]);								// skip steps
	double startTime=stod(argv[4]);				// t_inf
	double endTime=stod(argv[5]);					// time of a block
	double dt=stod(argv[6]);							// dt in the block file
	int skip=int(skipTime/dt);				//
	int blockSize=int(endTime/dt);				//
	int skipInBlock=int(startTime/endTime);		//

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
		sprintf(dump, argv[7], i+calculationNumber); // time x y z vx vy vz
		ifstream stream(dump);
		string str;
		int lineIndex=0;
		// loop of ion file
		while(getline(stream,str)) {
			istringstream stream(str);
			string tmp;
			if(lineIndex > skip){
				data d;
				while(getline(stream,tmp,' ')) d.MD.push_back(stod(tmp));
				datas.push_back(d);
			}
			lineIndex++;
		}
		int dataSize=datas.size();
		int N=(dataSize-blockSize)/dstep-1;
		Nblocks+=N;

		for(int k=0;k<N;k++){
			double X=datas[k*dstep].MD[1];	// initial position in this block
			double Y=datas[k*dstep].MD[2];
			double Z=datas[k*dstep].MD[3];

			double vX=datas[k*dstep].MD[4];
			double vY=datas[k*dstep].MD[5];
			double vZ=datas[k*dstep].MD[6];

			for(int j=0;j<blockSize;j++){
				int jj=k*dstep+j;

				double dx=datas[jj].MD[1]-X;
				double dy=datas[jj].MD[2]-Y;
				double dz=datas[jj].MD[3]-Z;
				MSDx[j]+=dx*dx;
				MSDy[j]+=dy*dy;
				MSDz[j]+=dz*dz;
				MSD[j]+=dx*dx+dy*dy+dz*dz;

				/*velocityAve[0][j]+=dx;
				velocityAve[1][j]+=dy;
				velocityAve[2][j]+=dz;*/

				VAF[j]+=vX*(datas[jj].MD[4])+vY*(datas[jj].MD[5])+vZ*(datas[jj].MD[6]);
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

	sprintf(dump, "TIME_MSD_VAF.%d", calculationNumber);
	FILE*f=fopen(dump, "w");
	fclose(f);
	for(int i=0;i<blockSize;i++){
		f=fopen(dump, "a");
		fprintf(f, "%f\t%f\t%e\n", dt*i*1e9,MSD[i],VAF[i]);
		fclose(f);
	}


	double Dvaf=dt/9.0*integral(VAF,blockSize)*1e14;
	double Dmsd=slope(MSD,skipInBlock,blockSize-1,dt)/6.0*1e-20*1e4;
	double Dmsdx=slope(MSDx,skipInBlock,blockSize-1,dt)/2.0*1e-20*1e4;
	double Dmsdy=slope(MSDy,skipInBlock,blockSize-1,dt)/2.0*1e-20*1e4;
	double Dmsdz=slope(MSDz,skipInBlock,blockSize-1,dt)/2.0*1e-20*1e4;
	double coeff=e/kb/T;
	double Kvaf=Dvaf*coeff;
	double Kmsd=Dmsd*coeff;
	double Kmsdx=Dmsdx*coeff;
	double Kmsdy=Dmsdy*coeff;
	double Kmsdz=Dmsdz*coeff;

	sprintf(dump, "DiffusionCoefficients.%d", int(calculationNumber));
	f=fopen(dump, "w");
	fprintf(f, "Dmsdx\tDmdsy\tDmsdz\tDmsd\tDvaf\te/kb/T\tCCScoeff\n");
	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%e\n", Dmsdx,Dmsdy,Dmsdz,Dmsd,Dvaf,coeff,3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)*1e4*1e20);
	fclose(f);

	coeff=3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)*1e4*1e20;
	cout << Dvaf << endl;
	cout << Dmsd << endl;
	cout << Dmsdx << endl;
	cout << Dmsdy << endl;
	cout << Dmsdz << endl;
	cout <<"K="<<Kmsd<<endl;
}
