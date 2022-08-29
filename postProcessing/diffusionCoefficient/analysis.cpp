#include "analysis.hpp"


int main(int argc,char *argv[]){
	double x,y,z,X,Y,Z,dx,dy,dz,vX,vY,vZ;
	int i,j,k,jj;
	int calculationNumber=stoi(argv[1]);
	cout<<calculationNumber<<endl;
	char dump[100];
	char dump_gas[100];
	for(i=0;i<pl;i++){
		sprintf(dump, "ion_300_%d.dat", i+calculationNumber);
		X=Y=Z=0;
		ifstream stream(dump);
		string str;
		j=k=0;
		while(getline(stream,str)) {
			k=0;
			istringstream stream(str);
			string tmp;
			if(j>skip){
				while(getline(stream,tmp,' ')) {
					data[j-skip][k]=stod(tmp);
					k++;
				}
			}
			j++;
			if(j>total_step+skip+1) break;
		}
		ifstream stream_gas(dump_gas);
		j=k=0;

		for(k=0;k<N;k++){
			X=data[k*dstep][1];
			Y=data[k*dstep][2];
			Z=data[k*dstep][3];
			for(j=0;j<sample_step;j++){
				jj=k*dstep+j;
				dx=data[jj][1]-X;
				dy=data[jj][2]-Y;
				dz=data[jj][3]-Z;
				MSDx[j]+=dx*dx;
				MSDy[j]+=dy*dy;
				MSDz[j]+=dz*dz;
				MSD[j]+=dx*dx+dy*dy+dz*dz;
				velocityAve[0][j]+=dx;
				velocityAve[1][j]+=dy;
				velocityAve[2][j]+=dz;
			}
		}
		vX=vY=vZ=0;
		for(k=0;k<N;k++){
			vX=data[k*dstep][4]-data_gas[k*dstep][4];
			vY=data[k*dstep][5]-data_gas[k*dstep][5];
			vZ=data[k*dstep][6]-data_gas[k*dstep][6];
			for(j=0;j<sample_step;j++){
				jj=k*dstep+j;
				VAF[j]+=vX*(data[jj][4]-data_gas[jj][4])+vY*(data[jj][5]-data_gas[jj][5])+vZ*(data[jj][6]-data_gas[jj][6]);
			}
		}
	}

	for(k=0;k<sample_step;k++)	VAF[k]/=double(N*pl);
	double sim2=0.0;
	for(k=1;k<(sample_step)/2-1;k++) sim2+=VAF[2*k];
	double sim3=0.0;
	for(k=1;k<(sample_step)/2;k++) sim3+=VAF[2*k-1];
	sprintf(dump, "TIME_MSD_VAF.%d", calculationNumber);
	FILE*f=fopen(dump, "w");
	fclose(f);
	for(j=0;j<sample_step;j++){
		f=fopen(dump, "a");
		fprintf(f, "%f\t%f\t%e\n", data[j][0],MSD[j]/N/pl,VAF[j]);
		fclose(f);
	}
	double Dvaf=dt/9.0*(VAF[0]+VAF[sample_step-1]+2*sim2+4*sim3)*1e14;
	double Dmsd=(MSD[start]-MSD[sample_step-1])/((start-sample_step)*dt)/(N*pl)/6.0*1e-20*1e4;
	double Dmsdx=(MSDx[start]-MSDx[sample_step-1])/((start-sample_step)*dt)/(N*pl)/2.0*1e-20*1e4;
	double Dmsdy=(MSDy[start]-MSDy[sample_step-1])/((start-sample_step)*dt)/(N*pl)/2.0*1e-20*1e4;
	double Dmsdz=(MSDz[start]-MSDz[sample_step-1])/((start-sample_step)*dt)/(N*pl)/2.0*1e-20*1e4;
	double Kvaf=Dvaf/kb/T*e;
	double Kmsd=Dmsd/kb/T*e;
	double Kmsdx=Dmsdx/kb/T*e;
	double Kmsdy=Dmsdy/kb/T*e;
	double Kmsdz=Dmsdz/kb/T*e;
	sprintf(dump, "Dx_Dy_Dz_DMSD_DVAF_Kcoff_dpcoff.%d", int(T));
	f=fopen(dump, "a");
	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%e\n", Dmsdx,Dmsdy,Dmsdz,Dmsd,Dvaf,e/kb/T,3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)*1e4*1e20);
	fclose(f);

	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)/Kvaf*1e4*1e20 << endl;
	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)/Kmsd*1e4*1e20 << endl;
	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)/Kmsdx*1e4*1e20 << endl;
	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)/Kmsdy*1e4*1e20 << endl;
	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)/Kmsdz*1e4*1e20 << endl;
	cout << 3.0*e/16.0/(P/kb/T)*sqrt(2.0*M_PI/mred/kb/T)*1e4*1e20 << endl;
	cout <<"K="<<Kmsd<<endl;
}
