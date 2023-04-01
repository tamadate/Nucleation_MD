#include "observer.hpp"

void
Observer::display(void){
	computeGasProps();
	computeIonProps();
	computeVaporProps();
    double virial=0;//ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
    int Ngas=vars->gases.size();
    double gaspress=(kb*Ngas*T_g + vars->totalVirial/3.0*6.95e-21)/(con->V*1e-30);
    double U = vars->Usum();
	double Kin=Kion+Kin_g+Kin_v;
	double Kout=Kout_g+Kout_v;
    std::cout << "----------------------TIME = " << vars->time/1000.0 << " ps-------------------------" << endl;
    cout<<"Inside propeties"<<endl;
    printf("  Kion = %1.2e  Kgas = %1.2e  Kvap = %1.2e\n  Tion = %1.2f  Tgas = %1.2f  Tvap = %1.2f\n  Uion = %1.2e  Ugas = %1.2e  Uvap = %1.2e\n  Ugi = %1.2e  Ugg = %1.2e  Uvi = %1.2e	\n  Uvg = %1.2e	Uvv= %1.2e  \n",
		Kion, Kin_g, Kin_v, Tion, Tin_g, Tin_v, vars->U.Uion, vars->U.Ugas, vars->U.Uvap, vars->U.Ugi, vars->U.Ugg,	vars->U.Uvi, vars->U.Uvg, vars->U.Uvv);
    cout<<"Out side propeties"<<endl;
    printf("  Kgas = %1.2e  Tgas = %1.2f  Ugas = %1.2e	\n  Kvap = %1.2e  Tvap = %1.2f  Uvap = %1.2e	\n  Kout = %1.2e    Uout = %1.2e	\n",
		Kout_g, Tout_g, 0.0, Kout_v, Tout_v, 0.0, Kout, 0.0);
    cout<<"System propeties"<<endl;
    printf("  K = %1.2e	U = %1.2e	Press = %f\n",Kin+Kout, U, gaspress/101300.0);
    cout<<"Times"<<endl;
    printf("  tion  = %1.1f s	tgas = %1.1f s 	tvap = %1.1f s\n",vars->times.tion,vars->times.tgas,vars->times.tvap);
    printf("  tvi   = %1.1f s	tgi  = %1.1f s	tvg  = %1.1f s	tvv  = %1.1f s\n",vars->times.tvi,vars->times.tgi,vars->times.tvg,vars->times.tvv);
    printf("  tpair = %1.1f s	tpos = %1.1f s	tvel = %1.1f s	tetc = %1.1f s\n",vars->times.tpair,vars->times.tpos,vars->times.tvel,vars->times.tetc);
    printf("  tpot  = %1.1f s	ttot = %1.1f s\n",(vars->times.tvi+vars->times.tgi+vars->times.tvv+vars->times.tvg+vars->times.tion+vars->times.tgas+vars->times.tvap), omp_get_wtime()-startTime);

    cout <<endl;

    FILE*f=fopen(fileKinetic, "a");
    fprintf(f,"%e %e %e %e %e %e\n",vars->time,Kion,Kin_g,Kin_v,Kout_g,Kout_v);
    fclose(f);

    f=fopen(filePotential, "a");
    fprintf(f,"%e %e %e %e %e %e %e %e %e\n",vars->time,vars->U.Uion,vars->U.Ugas,vars->U.Uvap,vars->U.Ugi,vars->U.Ugg,vars->U.Uvg,vars->U.Uvi,vars->U.Uvv);
    fclose(f);
}
