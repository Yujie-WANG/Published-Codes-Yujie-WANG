#ifndef WEIBULL_H
#define WEIBULL_H

#include "math.h"

double absolute(double a){
	if(a >= 0){
		return a;
	}
	else{
		return -a;
	}
}

double median(double list[200], int length){
	double result, copy[200], tmp;
	int i, j;
	for(i=0; i<length; i++){
		copy[i] = list[i];
	}
	for(i=0; i<length-1; i++){
		for(j=0; j<length-1-i; j++){
			if(copy[j] > copy[j+1]){
				tmp = copy[j];
				copy[j] = copy[j+1];
				copy[j+1] = tmp;
			}
		}
	}
	if(length%2 == 0){
		result = (copy[length/2] + copy[length/2 - 1]) /2.0;
	}
	else{
		result = copy[length/2];
	}
	return result;
}

double minimal(double list[200], int length){
	double min = list[0];
	int i;
	for(i=0; i<length; i++){
		if(list[i] < min){
			min = list[i];
		}
	}
	return min;
}

int minimal_site(double list[200], int length){
	double min = list[0];
	int i, min_site=0;
	for(i=0; i<length; i++){
		if(list[i] < min){
			min = list[i];
			min_site = i;
		}
	}
	return min_site;
}

double maximal(double list[200], int length){
	double max = list[0];
	int i;
	for(i=0; i<length; i++){
		if(list[i] > max){
			max = list[i];
		}
	}
	return max;
}

int maximal_site(double list[200], int length){
	double max = list[0];
	int i, max_site=0;
	for(i=0; i<length; i++){
		if(list[i] > max){
			max = list[i];
			max_site = i;
		}
	}
	return max_site;
}

double norm_k(double mean, double sd, double targ){
	double f;
	f = 1.0/(sd*sqrt(2.0*3.141592654)) * exp(-pow(targ-mean,2)/(2.0*pow(sd,2)));
	return f;
}

double weibull_k(double b, double c, double t){
	double ratio;
	ratio = exp(-pow(t/b,c));
	return ratio;
}

double weibull_to_e(double b, double c, double p1, double p2){
	int i,N=30;
	double sum=0.0,tmp,dp,list[N+1];
	dp = (p2-p1) / N;
	for(i=0; i<N+1; i++){
		tmp = p1 + (p2-p1)*i/N;
		list[i] = weibull_k(b,c,tmp);
		if(i==0 || i==N){
			sum = sum + list[i]*0.5*dp;
		}
		else{
			sum = sum + list[i]*dp;
		}
	}
	return sum;
}

#endif
