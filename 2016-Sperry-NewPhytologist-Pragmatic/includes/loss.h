#ifndef LOSS_H
#define LOSS_H

#include "supply.h"

double loss_p_reg_max(double data[3][101]){
	int i;
	double p_reg,p_reg_max=0.0;
	for(i=0; i<101; i++){
		p_reg = data[2][i]/data[2][0] * (data[0][i]-data[0][0]) + data[0][0];
		if(p_reg >= p_reg_max){
			p_reg_max = p_reg;
		}
		else{
			break;
		}
	}
	return p_reg_max;
}

double loss_e_reg_max(double data[3][101], double weibulls[6], double kmaxs[3], double k_soil, double p_soil, double bulk[3]){
	double e_reg_max;
	int i;
	double p_reg,p_reg_max=0.0;
	for(i=0; i<101; i++){
		p_reg = data[2][i]/data[2][0] * (data[0][i]-data[0][0]) + data[0][0];
		if(p_reg >= p_reg_max){
			p_reg_max = p_reg;
		}
		else{
			break;
		}
	}
	e_reg_max = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_reg_max,bulk);
	return e_reg_max;
}

double loss_e_reg_leaf_stem(double data[3][101], double weibulls[6], double kmaxs[3], double p_soil){
	double e_reg_max;
	int i;
	double p_reg,p_reg_max=0.0;
	for(i=0; i<101; i++){
		p_reg = data[2][i]/data[2][0] * (data[0][i]-data[0][0]) + data[0][0];
		if(p_reg >= p_reg_max){
			p_reg_max = p_reg;
		}
		else{
			break;
		}
	}
	e_reg_max = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_soil,p_reg_max);
	return e_reg_max;
}

double loss_point_one(double gmax, double vpd, double data[3][101], double weibulls[6], double kmaxs[3], double k_soil, double p_soil, double bulk[3]){
	double e_reg,e_reg_bef,e_reg_max,s_max;
	int i,max_location;
	double p_reg,p_reg_max=0.0;
	for(i=0; i<101; i++){
		p_reg = data[2][i]/data[2][0] * (data[0][i]-data[0][0]) + data[0][0];
		if(p_reg >= p_reg_max){
			p_reg_max = p_reg;
			//p_reg_bef = data[0][i];
			max_location = i;
		}
		else{
			break;
		}
	}
	e_reg_bef = data[1][max_location];
	e_reg_max = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_reg_max,bulk);
	s_max = data[2][0];
	double g_vpd;
	if(vpd == 0.0){
		g_vpd = 0.0;
	}
	else if(gmax*vpd/100.0 >= e_reg_bef){
		g_vpd = e_reg_max / (vpd/100.0);
	}
	else{
		double e1, e2, e_targ = gmax * vpd/100.0;
		double p_leaf=p_soil, dp=1E-4,slope;
		int i=0;
		while(1){
			e1 = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf,bulk) - e_targ;
			if(absolute(e1) < 1E-6){
				break;
			}
			e2 = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf+dp,bulk) - e_targ;
			slope = (e2 - e1) / dp;
			p_leaf = p_leaf - e1/slope;
			if(p_leaf < 0){
				p_leaf = 0.0;
			}
			i++;
			if(i > 100){
				//printf("There is something wrong with the loss_point_one function!\n");
				break;
			}
		}
		double e_tree, e_tree_dp, e_dp = 0.001, s_reg;
		e_tree = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf,bulk);
		e_tree_dp = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf+e_dp,bulk);
		s_reg = (e_tree_dp - e_tree) / e_dp;
		p_reg = s_reg/s_max * (p_leaf-p_soil) + p_soil;
		e_reg = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_reg,bulk);
		g_vpd = e_reg / (vpd/100.0);
	}
	return g_vpd;
}

double loss_point_multi(double gmax, double vpd, double data[3][101], double weibulls[6], double kmaxs[3], double k_root[20], double k_soil[20], double p_soil[20], double bulk[3], double p_predawn, int root_layers){
	double e_reg,e_reg_bef,e_reg_max,s_max;
	int i,max_location;
	double p_reg,p_reg_max=0.0;
	for(i=0; i<101; i++){
		p_reg = data[2][i]/data[2][0] * (data[0][i]-data[0][0]) + data[0][0];
		if(p_reg >= p_reg_max){
			p_reg_max = p_reg;
			//p_reg_bef = data[0][i];
			max_location = i;
		}
		else{
			break;
		}
	}
	e_reg_bef = data[1][max_location];
	e_reg_max = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_reg_max,bulk,root_layers);
	s_max = data[2][0];
	double g_vpd;
	if(vpd == 0.0){
		g_vpd = 0.0;
	}
	else if(gmax*vpd/100.0 >= e_reg_bef){
		g_vpd = e_reg_max / (vpd/100.0);
	}
	else{
		double e1, e2, e_targ = gmax * vpd/100.0;
		double p_leaf, dp=1E-6,slope;
		p_leaf = p_predawn;
		int i=0;
		while(1){
			e1 = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_leaf,bulk,root_layers) - e_targ;
			if(absolute(e1) < 1E-6){
				break;
			}
			e2 = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_leaf+dp,bulk,root_layers) - e_targ;
			slope = (e2 - e1) / dp;
			p_leaf = p_leaf - e1/slope;
			if(p_leaf < 0){
				p_leaf = 0.0;
			}
			i++;
			if(i > 100){
				//printf("There is something wrong with the loss_point_one function!\n");
				break;
			}
		}
		double e_tree, e_tree_dp, e_dp = 0.001, s_reg;
		e_tree = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_leaf,bulk,root_layers);
		e_tree_dp = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_leaf+e_dp,bulk,root_layers);
		s_reg = (e_tree_dp - e_tree) / e_dp;
		p_reg = s_reg/s_max * (p_leaf-p_predawn) + p_predawn;
		e_reg = supply_point_multi(weibulls,kmaxs,k_root,k_soil,p_soil,p_reg,bulk,root_layers);
		g_vpd = e_reg / (vpd/100.0);
	}
	return g_vpd;
}

#endif
