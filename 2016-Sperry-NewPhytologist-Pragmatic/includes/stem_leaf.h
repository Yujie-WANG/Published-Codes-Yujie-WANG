#ifndef STEM_LEAF_H
#define STEM_LEAF_H

#include "unistd.h"
#include "weibull.h"

double stem_leaf_flow(double bs, double cs, double ks, double bl, double cl, double kl, double p_root, double p_leaf){
	double e_above;
	if(p_leaf <= p_root){
		e_above = 0.0;
	}
	else{
		double p_stem, dp, e_stem, e_leaf, e1, e2, de, slope;
		p_stem = p_root;
		dp = 1E-6;
		int i=0;
		while(1){
			e_stem = weibull_to_e(bs,cs,p_root,p_stem) * ks;
			e_leaf = weibull_to_e(bl,cl,p_stem,p_leaf) * kl;
			e1 = e_stem - e_leaf;
			if(absolute(e1) < 1E-6){
				break;
			}
			e2 = weibull_to_e(bs,cs,p_root,p_stem+dp) * ks - weibull_to_e(bl,cl,p_stem+dp,p_leaf) * kl;
			de = e2 - e1;
			slope = de / dp;
			//printf("P: %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",p_stem,e1,e2,de,dp,slope,i); sleep(3);
			p_stem = p_stem - e1/slope;
			if(p_stem < p_root){
				p_stem = p_root;
			}
			if(p_stem > p_leaf){
				p_stem = p_leaf;
			}
			i++;
			if(i > 100){
				//printf("There is something wrong with the stem_leaf_flow function!\n");
				//printf("P: %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",p_root,p_leaf,e1,e2,de,dp,slope,i); sleep(300);
				break;
			}
		}
		e_above = (e_stem + e_leaf) / 2.0;
	}
	return e_above;
}

#endif
