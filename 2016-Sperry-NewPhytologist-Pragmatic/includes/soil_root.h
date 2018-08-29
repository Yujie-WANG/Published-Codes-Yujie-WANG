#ifndef SOIL_ROOT_H
#define SOIL_ROOT_H

#include "math.h"

#include "weibull.h"

double get_soil_r_fraction(double kmaxs[3], double weibulls[6], double k_soil, double bulk[3]){
	double max_r_fraction;
	double p_crit, p_soil, f_soil, f_list[101];
	double b_alpha,b_n,b_m,b_theta,b_fract, k_b, k_r, k_s, k_l;
	b_alpha = bulk[0]; b_n = bulk[1]; b_m = bulk[2];
	p_crit = weibulls[4] * pow(log(1000.0),1.0/weibulls[5]);
	int i;
	for(i=0; i<101; i++){
		p_soil = p_crit/100.0 * i;
		b_theta = pow((1.0/(1.0 + pow(b_alpha*p_soil,b_n))),b_m);
		b_fract = sqrt(b_theta) * pow((1.0 - pow(1.0-pow(b_theta,1.0/b_m),b_m)),2.0);
		k_b = k_soil * b_fract;
		k_r = kmaxs[0] * weibull_k(weibulls[0], weibulls[1], p_soil);
		k_s = kmaxs[1] * weibull_k(weibulls[2], weibulls[3], p_soil);
		k_l = kmaxs[2] * weibull_k(weibulls[4], weibulls[5], p_soil);
		f_soil = (1.0/k_b) / (1.0/k_r + 1.0/k_s + 1.0/k_l + 1.0/k_b);
		f_list[i] = f_soil;
		//printf("%lf\t",f_soil);
	}
	max_r_fraction = maximal(f_list,101);
	return max_r_fraction;
}

double soil_pressure(double ks, double flow, double p_soil, double bulk[3]){
	double b_alpha,b_n,b_m,b_theta,b_fract,b_k;
	double dp = 0.0;
	b_alpha = bulk[0];
	b_n = bulk[1];
	b_m = bulk[2];
	double r_root=1.0,r_bulk=10.0,r_layers[11];
	int i;
	for(i=0; i<11; i++){
		r_layers[i] = r_bulk - (r_bulk-r_root)*i/10;
	}
	for(i=0; i<10; i++){
		b_theta = pow((1.0/(1.0 + pow(b_alpha*p_soil,b_n))),b_m);
		b_fract = sqrt(b_theta) * pow((1.0 - pow(1.0-pow(b_theta,1.0/b_m),b_m)),2.0);
		b_k = ks * b_fract * log(r_bulk/r_root) / log(r_layers[i]/r_layers[i+1]);
		dp = dp + flow/b_k;
	}
	return p_soil + dp;
}

double soil_root_flow_one(double br, double cr, double kr, double ks, double p_soil, double p_root, double bulk[3]){
	double e_below;
	if(p_root <= p_soil){
		e_below = 0.0;
	}
	else{
		double p_rizo, tmp_p_rizo, dp, e_below_dp, p1, p2, slope;
		p_rizo = p_soil;
		dp = (p_root-p_soil) * 1E-6;
		int i = 0;
		while(1){
			e_below = weibull_to_e(br,cr,p_rizo,p_root) * kr;
			tmp_p_rizo = soil_pressure(ks,e_below,p_soil,bulk);
			p1 = p_rizo - tmp_p_rizo;
			//printf("%d\t%lf\t%lf\t%lf\n",i,p_rizo,tmp_p_rizo,p1);
			if(absolute(p1) < 1E-6){
				break;
			}
			e_below_dp = weibull_to_e(br,cr,p_rizo+dp,p_root) * kr;
			tmp_p_rizo = soil_pressure(ks,e_below_dp,p_soil,bulk);
			p2 = p_rizo + dp - tmp_p_rizo;
			slope = (p2 - p1) / dp;
			p_rizo = p_rizo - p1/slope;
			if(p_rizo > p_root){
				p_rizo = p_root;
			}
			else if(p_rizo < p_soil){
				p_rizo = p_soil;
			}
			i++;
			if(i > 100){
				break;
				//printf("There is something wrong with the soil_root_flow_one function!\n");
			}
		}
	}
	return e_below;
}

double soil_root_flow_multi(double br, double cr, double kr[20], double ks[20], double p_soil[20], double p_root, double bulk[3], int root_layers){
	double e_below=0, e_layer;
	int i;
	for(i=0; i<root_layers; i++){
		double p_rizo, tmp_p_rizo, dp, e_layer_dp, p1, p2, slope;
		p_rizo = p_soil[i];
		dp = 1E-6;
		int j = 0;
		while(1){
			e_layer = weibull_to_e(br,cr,p_rizo,p_root) * kr[i];
			tmp_p_rizo = soil_pressure(ks[i],e_layer,p_soil[i],bulk);
			p1 = p_rizo - tmp_p_rizo;
			//printf("%d\t%lf\t%lf\t%lf\n",j,p_rizo,tmp_p_rizo,p1);
			if(absolute(p1) < 1E-6){
				break;
			}
			e_layer_dp = weibull_to_e(br,cr,p_rizo+dp,p_root) * kr[i];
			tmp_p_rizo = soil_pressure(ks[i],e_layer_dp,p_soil[i],bulk);
			p2 = p_rizo + dp - tmp_p_rizo;
			slope = (p2 - p1) / dp;
			p_rizo = p_rizo - p1/slope;
			if(p_rizo < 0.0){
				p_rizo = 0.0;
			}
			j++;
			if(j > 100){
				//printf("There is something wrong with the soil_root_flow_one function!\n");
				break;
			}
			e_below = e_below + e_layer;
		}
	}
	return e_below;
}

double get_p_predawn(double br, double cr, double kr[20], double ks[20], double p_soil[20], double bulk[3], int root_layers){
	double p_predawn, dp, e_below, e_below_dp, slope;
	p_predawn = minimal(p_soil,root_layers);
	dp = 1E-6;
	int i=0;
	while(1){
		e_below = soil_root_flow_multi(br,cr,kr,ks,p_soil,p_predawn,bulk,root_layers);
		//printf("P predawn: %d\t%lf\n",i,e_below);
		if(absolute(e_below) < 1E-6){
			break;
		}
		e_below_dp = soil_root_flow_multi(br,cr,kr,ks,p_soil,p_predawn+dp,bulk,root_layers);
		slope = (e_below_dp - e_below) / dp;
		p_predawn = p_predawn - e_below/slope;
		i++;
		if(i > 100){
			//printf("There is something wrong when computing p_predawn!\n");
			break;
		}
	}
	return p_predawn;
}

#endif
