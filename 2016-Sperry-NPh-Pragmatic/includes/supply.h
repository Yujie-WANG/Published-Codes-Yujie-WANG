#ifndef SUPPLY_H
#define SUPPLY_H

#include "math.h"
#include "soil_root.h"
#include "stem_leaf.h"

double root_depth(double percent, double beta){
	double depth;
	depth = log(1.0-percent) / log(beta);
	return depth;
}

void root_layers_fraction(int root_layers, double root_beta, double root_alpha, double root_fraction[20]){
	int i;
	double layer_boundary[20], layer_depth[20], root_length[20];
	layer_boundary[0] = 0.0;
	double sum_k=0;
	for(i=0; i<root_layers; i++){
		layer_boundary[i+1] = root_depth(0.995/root_layers*(i+1),root_beta);
		layer_depth[i] = layer_boundary[i+1] - layer_boundary[i];
		double tmp_v, tmp_h;
		tmp_v = root_depth(0.995/root_layers*(0.5+i),root_beta);
		tmp_h = 0.995/root_layers / layer_depth[i] / (2.0 * 3.141593) * root_alpha;
		root_length[i] = tmp_v + tmp_h;
		sum_k = sum_k + 1.0/root_length[i];
	}
	for(i=0; i<root_layers; i++){
		root_fraction[i] = 1.0/root_length[i] / sum_k;
	}
}

double supply_point_one(double weibulls[6], double kmaxs[3], double k_soil, double p_soil, double p_leaf, double bulk[3]){
	double p_root, dp=1E-6, e_above, e_below, e_tree, e_above_dp, e_below_dp, e1, e2, slope;
	p_root = p_soil;
	int i=0;
	while(1){
		e_above = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_root,p_leaf);
		e_below = soil_root_flow_one(weibulls[0],weibulls[1],kmaxs[0],k_soil,p_soil,p_root,bulk);
		e1 = e_above - e_below;
		if(absolute(e1) < 1E-6){
			e_tree = (e_above+e_below) / 2.0;
			break;
		}
		e_above_dp = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_root+dp,p_leaf);
		e_below_dp = soil_root_flow_one(weibulls[0],weibulls[1],kmaxs[0],k_soil,p_soil,p_root+dp,bulk);
		e2 = e_above_dp - e_below_dp;
		slope = (e2 - e1) / dp;
		p_root = p_root - e1/slope;
		if(p_root < 0){
			p_root = 0.0;
		}
		i++;
		if(i > 100){
			//printf("There is something wrong in supply_point_one function!\n");
			break;
		}
	}
	return e_tree;
}

void supply_curve_stem_leaf(double data[3][101], double weibulls[6], double kmaxs[3], double p_soil){
	double bl,cl,p_crit,p_leaf,dp=0.001;
	double e_tree,e_tree_dp,s_tree;
	bl = weibulls[4];
	cl = weibulls[5];
	p_crit = bl * pow(log(1000.0),1.0/cl);
	int i;
	for(i=0; i<101; i++){
		//printf("%d ",i);
		p_leaf = p_soil + (p_crit-p_soil)*i/100;
		e_tree = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_soil,p_leaf);
		e_tree_dp = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_soil,p_leaf+dp);
		s_tree = (e_tree_dp-e_tree) / dp;
		//printf("P,E,S are %lf,%lf,%lf\n",p_leaf,e_tree,s_tree);
		data[0][i] = p_leaf;
		data[1][i] = e_tree;
		data[2][i] = s_tree;
	}
	//printf("Supply Curve Finished!\n");
}

void supply_curve_one(double data[3][101], double weibulls[6], double kmaxs[3], double k_soil, double p_soil, double bulk[3]){
	double bl,cl,p_crit,p_leaf,dp;
	double e_tree,e_tree_dp,s_tree;
	bl = weibulls[4];
	cl = weibulls[5];
	p_crit = bl * pow(log(1000.0),1.0/cl);
	dp = p_crit / 1000.0;
	int i;
	for(i=0; i<101; i++){
		//printf("%d ",i);
		p_leaf = p_soil + (p_crit-p_soil)*i/100;
		e_tree = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf,bulk);
		e_tree_dp = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_leaf+dp,bulk);
		s_tree = (e_tree_dp-e_tree) / dp;
		//printf("P,E,S are %lf,%lf,%lf\n",p_leaf,e_tree,s_tree);
		data[0][i] = p_leaf;
		data[1][i] = e_tree;
		data[2][i] = s_tree;
	}
	//printf("\nSupply Curve Finished!\n\n\n");
}

double supply_point_multi(double weibulls[6], double kmaxs[3], double k_root[20], double k_soil[20], double p_soil[20], double p_leaf, double bulk[3], int root_layers){
	double p_root, dp=1E-6, e_above, e_below, e_tree, e_above_dp, e_below_dp, e1, e2, slope;
	p_root = minimal(p_soil,root_layers);
	int i=0;
	while(1){
		e_above = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_root,p_leaf);
		e_below = soil_root_flow_multi(weibulls[0],weibulls[1],k_root,k_soil,p_soil,p_root,bulk,root_layers);
		e1 = e_above - e_below;
		if(absolute(e1) < 1E-6){
			e_tree = (e_above+e_below) / 2.0;
			break;
		}
		e_above_dp = stem_leaf_flow(weibulls[2],weibulls[3],kmaxs[1],weibulls[4],weibulls[5],kmaxs[2],p_root+dp,p_leaf);
		e_below_dp = soil_root_flow_multi(weibulls[0],weibulls[1],k_root,k_soil,p_soil,p_root+dp,bulk,root_layers);
		e2 = e_above_dp - e_below_dp;
		slope = (e2 - e1) / dp;
		p_root = p_root - e1/slope;
		if(p_root < 0){
			p_root = 0.0;
		}
		i++;
		if(i > 100){
			//printf("There is something wrong in supply_point_one function!\n");
			break;
		}
	}
	return e_tree;
}

void supply_curve_multi(double data[3][101], double weibulls[6], double kmaxs[3], double k_soil, double p_soil[20], double bulk[3], int root_layers, double root_beta, double root_alpha){
	double bl,cl,p_crit,p_leaf,p_predawn,dp;
	double e_tree,e_tree_dp,s_tree;
	double root_fraction[20],kr_list[20],ks_list[20];
	int i;
	root_layers_fraction(root_layers,root_beta,root_alpha,root_fraction);
	for(i=0; i<root_layers; i++){
		kr_list[i] = root_fraction[i] * kmaxs[0];
		ks_list[i] = 1.0/root_layers * k_soil;
	}
	for(i=0; i<root_layers; i++){
		printf("%lf\t%lf\t%lf\n",kr_list[i],ks_list[i],p_soil[i]);
	}
	p_predawn = get_p_predawn(weibulls[0],weibulls[1],kr_list,ks_list,p_soil,bulk,root_layers);
	printf("\nP predawn is %lf\n\n",p_predawn);
	bl = weibulls[4];
	cl = weibulls[5];
	p_crit = bl * pow(log(1000.0),1.0/cl);
	dp = (p_crit-p_predawn) *4E-3;
	for(i=0; i<101; i++){
		//printf("%d\n",i);
		p_leaf = p_predawn + (p_crit-p_predawn)*i/100;
		e_tree = supply_point_multi(weibulls,kmaxs,kr_list,ks_list,p_soil,p_leaf,bulk,root_layers);
		e_tree_dp = supply_point_multi(weibulls,kmaxs,kr_list,ks_list,p_soil,p_leaf+dp,bulk,root_layers);
		s_tree = (e_tree_dp-e_tree) / dp;
		printf("P,E,S are\t%lf\t%lf\t%lf\n",p_leaf,e_tree,s_tree);
		data[0][i] = p_leaf;
		data[1][i] = e_tree;
		data[2][i] = s_tree;
	}
	//printf("Supply Curve Finished!\n");
}

double find_p_vpd_one(double weibulls[6], double kmaxs[3], double k_soil, double p_soil, double bulk[3], double e_targ){
	double p_reg, dp=1E-6, e_reg, e_reg_dp, slope;
	p_reg = p_soil;
	double bl, cl, p_crit;
	bl = weibulls[4];
	cl = weibulls[5];
	p_crit = bl * pow(log(1000.0),1.0/cl);
	int i=0;
	while(1){
		e_reg = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_reg,bulk) - e_targ;
		if(absolute(e_reg) < 1E-6){
			break;
		}
		e_reg_dp = supply_point_one(weibulls,kmaxs,k_soil,p_soil,p_reg+dp,bulk) - e_targ;
		slope = (e_reg_dp - e_reg) / dp;
		p_reg = p_reg - e_reg/slope;
		if(p_reg < 0){
			p_reg = 0.0;
		}
		if(p_reg > p_crit){
			p_reg = p_crit;
		}
		i++;
		if(i > 100){
			//printf("There is something wrong in find_p_vpd_one function!\n");
			//printf("Psoil: %lf\tP_crit %lf\tp_reg %lf\t Slope %lf\t%d\n",p_soil,p_crit,p_reg,slope,i);
			break;
		}
	}
	return p_reg;
}

#endif
