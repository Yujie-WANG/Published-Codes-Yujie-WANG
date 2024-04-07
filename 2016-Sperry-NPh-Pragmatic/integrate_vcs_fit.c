#include "stdio.h"
#include <stdlib.h>

#include "./includes/loss.h"
#include "./includes/soil_root.h"
#include "./includes/stem_leaf.h"
#include "./includes/supply.h"
#include "./includes/weibull.h"

int main(int agvc, char **agvs){
	if(agvc == 5){
		printf("B and values are inputed!\n");
		double fb,fc,lb,lc;
		fb = atof(agvs[1]);
		fc = atof(agvs[2]);
		lb = atof(agvs[3]);
		lc = atof(agvs[4]);
		double weibulls[6] = {fb,fc,fb,fc,fb,fc};
		double kmaxs[3] = {20.0,40.0,40.0};
		double k_soil = 2.07e12;
		double bulk[3] = {602.0419,1.48,0.324324};
		double data[3][101];
		double p_crit,p_predawn, e_reg_max;
		p_crit = lb * pow(log(1000.0),1.0/lc);
		//printf("%lf\n",p_crit);
		int i;
		double p_list[50], e_list[50];
		for(i=0; i< 50; i++){
			p_predawn = p_crit/50.0 * i;
			printf("%d ",i);
			p_list[i] = p_predawn;
			supply_curve_one(data,weibulls,kmaxs,k_soil,p_predawn,bulk);
			e_reg_max = loss_e_reg_max(data,weibulls,kmaxs,k_soil,p_predawn,bulk);
			e_list[i] = e_reg_max;
		}
		printf("\n");
		FILE *save_file;
		save_file = fopen("./output/single_vc_e.txt","w");
		for(i=0; i<50; i++){
			fprintf(save_file,"%lf\t%lf\n",p_list[i],e_list[i]);
		}
		fclose(save_file);
	}
	return 0;
}
