#include "stdio.h"

#include "./includes/loss.h"
#include "./includes/soil_root.h"
#include "./includes/stem_leaf.h"
#include "./includes/supply.h"
#include "./includes/weibull.h"

int main(){
	double weibulls[6] = {1.7,3.0,2.0,3.0,1.3,3.0};
	double kmaxs[3] = {20.0,40.0,40.0};
	double k_soil = 2.07e12;
	double bulk[3] = {602.0419,1.48,0.324324};
	double data[3][101];
	double bl,cl,p_crit,p_predawn, e_reg_max;
	bl = weibulls[4];
	cl = weibulls[5];
	p_crit = bl * pow(log(1000.0),1.0/cl);
	//printf("%lf\n",p_crit);
	int i;
	double p_list[50], e_list[50];
	for(i=0; i< 50; i++){
		p_predawn = p_crit/50.0 * i;
		p_list[i] = p_predawn;
		//supply_curve_one(data,weibulls,kmaxs,k_soil,p_predawn,bulk);
		supply_curve_stem_leaf(data,weibulls,kmaxs,p_predawn);
		//e_reg_max = loss_e_reg_max(data,weibulls,kmaxs,k_soil,p_predawn,bulk);
		e_reg_max = loss_e_reg_leaf_stem(data,weibulls,kmaxs,p_predawn);
		e_list[i] = e_reg_max/2.0; /*WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW*/
	}
	FILE *save_file;
	save_file = fopen("./output/multi_vc_e.txt","w");
	for(i=0; i<50; i++){
		fprintf(save_file,"%lf\t%lf\n",p_list[i],e_list[i]);
	}
	fclose(save_file);
	return 0;
}
