#pragma once
#include "rigid.hpp"
#define tolerance 1e-20

void
Rigid::updatePosition(Variables *vars,double dt){
    vars->times.tpos-=omp_get_wtime();
    if(rigid_pairs_ion.size()) recordOldPosition(vars->ions);
    for (auto &a : vars->ions) {
        a.qx += a.px * dt;
        a.qy += a.py * dt;
        a.qz += a.pz * dt;
        a.fx=a.fy=a.fz=0.0;
    }
    if(rigid_pairs_ion.size()) shake_ion(vars);

    for (auto &i : vars->vapor_in) {
        if(rigid_pairs_vapor[i].size()) recordOldPosition(vars->vapors[i].inAtoms);
        for (auto &a : vars->vapors[i].inAtoms){
            a.qx += a.px * dt;
            a.qy += a.py * dt;
            a.qz += a.pz * dt;
            a.fx=a.fy=a.fz=0.0;
        }
        if(rigid_pairs_vapor[i].size()) shake_vapor(vars,i);
    }

    for (auto &i : vars->gas_in) {
        if(rigid_pairs_gas[i].size()) recordOldPosition(vars->gases[i].inAtoms);
        for (auto &a : vars->gases[i].inAtoms){
            a.qx += a.px * dt;
            a.qy += a.py * dt;
            a.qz += a.pz * dt;
            a.fx=a.fy=a.fz=0.0;
        }
        if(rigid_pairs_gas[i].size()) shake_gas(vars,i);
    }
    vars->times.tpos+=omp_get_wtime();

}

void 
Rigid::recordOldPosition(std::vector<Atom> &target){
    r.clear();
    for (auto &a : target) {
        std::array<double,3> pos;
        pos[0]=a.qx;
        pos[1]=a.qx;
        pos[2]=a.qx;
        r.push_back(pos);
    }
}

void
Rigid::shake_ion(Variables *vars){
    int l=rigid_pairs_ion.size();
    bool converge=false;
    while(converge){
        converge=true;
        for(int k=0;k<l;k++){
            int i=rigid_pairs_ion[k].atom1;
            int j=rigid_pairs_ion[k].atom2;
            Atom &ai=vars->ions[i];
            Atom &aj=vars->ions[j];
            int dk=vars->btypes[rigid_pairs_ion[k].type].coeff[1];
            double rkx=r[i][0]-r[j][0];
            double rky=r[i][1]-r[j][1];
            double rkz=r[i][2]-r[j][2];
            double _rkx=ai.qx-aj.qx;
            double _rky=ai.qy-aj.qy;
            double _rkz=ai.qz-aj.qz;
            double _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dot=rkx*_rkx+rky*_rky+rkz*_rkz;
            double ramk;
            if(dot>1e-10){
                double mred_inv=1/ai.mass+1/aj.mass;
                double deno=2*mred_inv*dot;
                ramk=(dk*dk-_rk2)/deno;
            }
            double coeffi=ramk/ai.mass;
            double coeffj=ramk/aj.mass;
            ai.qx+=r[i][0]*coeffi;
            ai.qy+=r[i][1]*coeffi;
            ai.qz+=r[i][2]*coeffi;
            aj.qx-=r[j][0]*coeffj;
            aj.qy-=r[j][1]*coeffj;
            aj.qz-=r[j][2]*coeffj;
            double _rkx=ai.qx-aj.qx;
            _rky=ai.qy-aj.qy;
            _rkz=ai.qz-aj.qz;
            _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dr=_rk2-dk*dk;
            if(dr*dr>tolerance) converge=false;
        }
    }
}

void
Rigid::shake_vapor(Variables *vars, int vaporID){
    int l=rigid_pairs_vapor[vaporID].size();
    bool converge=false;
    while(converge){
        converge=true;
        for(int k=0;k<l;k++){
            int i=rigid_pairs_vapor[vaporID][k].atom1;
            int j=rigid_pairs_vapor[vaporID][k].atom2;
            Atom &ai=vars->vapors[vaporID].inAtoms[i];
            Atom &aj=vars->vapors[vaporID].inAtoms[j];
            int dk=vars->btypes[rigid_pairs_ion[k].type].coeff[1];
            double rkx=r[i][0]-r[j][0];
            double rky=r[i][1]-r[j][1];
            double rkz=r[i][2]-r[j][2];
            double _rkx=ai.qx-aj.qx;
            double _rky=ai.qy-aj.qy;
            double _rkz=ai.qz-aj.qz;
            double _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dot=rkx*_rkx+rky*_rky+rkz*_rkz;
            double ramk;
            if(dot>1e-10){
                double mred_inv=1/ai.mass+1/aj.mass;
                double deno=2*mred_inv*dot;
                ramk=(dk*dk-_rk2)/deno;
            }
            double coeffi=ramk/ai.mass;
            double coeffj=ramk/aj.mass;
            ai.qx+=r[i][0]*coeffi;
            ai.qy+=r[i][1]*coeffi;
            ai.qz+=r[i][2]*coeffi;
            aj.qx-=r[j][0]*coeffj;
            aj.qy-=r[j][1]*coeffj;
            aj.qz-=r[j][2]*coeffj;
            double _rkx=ai.qx-aj.qx;
            _rky=ai.qy-aj.qy;
            _rkz=ai.qz-aj.qz;
            _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dr=_rk2-dk*dk;
            if(dr*dr>tolerance) converge=false;
        }
    }

}

void
Rigid::shake_gas(Variables *vars, int gasID){
    int l=rigid_pairs_vapor[gasID].size();
    bool converge=false;
    while(converge){
        converge=true;
        for(int k=0;k<l;k++){
            int i=rigid_pairs_vapor[gasID][k].atom1;
            int j=rigid_pairs_vapor[gasID][k].atom2;
            Atom &ai=vars->vapors[gasID].inAtoms[i];
            Atom &aj=vars->vapors[gasID].inAtoms[j];
            int dk=vars->btypes[rigid_pairs_ion[k].type].coeff[1];
            double rkx=r[i][0]-r[j][0];
            double rky=r[i][1]-r[j][1];
            double rkz=r[i][2]-r[j][2];
            double _rkx=ai.qx-aj.qx;
            double _rky=ai.qy-aj.qy;
            double _rkz=ai.qz-aj.qz;
            double _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dot=rkx*_rkx+rky*_rky+rkz*_rkz;
            double ramk;
            if(dot>1e-10){
                double mred_inv=1/ai.mass+1/aj.mass;
                double deno=2*mred_inv*dot;
                ramk=(dk*dk-_rk2)/deno;
            }
            double coeffi=ramk/ai.mass;
            double coeffj=ramk/aj.mass;
            ai.qx+=r[i][0]*coeffi;
            ai.qy+=r[i][1]*coeffi;
            ai.qz+=r[i][2]*coeffi;
            aj.qx-=r[j][0]*coeffj;
            aj.qy-=r[j][1]*coeffj;
            aj.qz-=r[j][2]*coeffj;
            double _rkx=ai.qx-aj.qx;
            _rky=ai.qy-aj.qy;
            _rkz=ai.qz-aj.qz;
            _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dr=_rk2-dk*dk;
            if(dr*dr>tolerance) converge=false;
        }
    }

}