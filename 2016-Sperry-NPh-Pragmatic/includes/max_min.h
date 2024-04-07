#ifndef MAX_MIN_H
#define MAX_MIN_H

#include "stdio.h"

void insert(double curve[2][200], int length, double point[2]){
	int i,j,inserted=0;
	for(i=0; i<length; i++){
		if(point[0]<curve[0][i]){
			for(j=length; j>i; j--){
				curve[0][j] = curve[0][j-1];
				curve[1][j] = curve[1][j-1];
			}
			curve[0][i] = point[0];
			curve[1][i] = point[1];
			inserted = 1;
			break;
		}
	}
	if(inserted == 0){
		curve[0][i] = point[0];
		curve[1][i] = point[1];
		inserted = 1;
	}
}

int judge_convex(double list[50], int length){
	int i,judge=0;
	double multi;
	if(length>=3){
		for(i=1; i<length-1;i++){
			multi = (list[i]-list[i-1]) * (list[i+1]-list[i]);
			if(multi < 0){
				judge = 1;
				break;
			}
		}
	}
	return judge;
}

int get_optimal_gmax(double curve[2][200], int length){
	int i,max=0;
	for(i=1; i<length; i++){
		if(curve[1][i] > curve[1][max]){
			max = i;
		}
	}
	return max;
}

void get_max_2(double curve[2][200], int length, double maxs[2]){
	int max1,max2;
	max1 = 0; max2=0;
	int i;
	for(i=1; i<length; i++){
		if(curve[1][i] > curve[1][max1]){
			max1 = i;
		}
	}
	if(max1 == 0){
		max2 = 1;
	}
	for(i=1; i<length; i++){
		if(curve[1][i] > curve[1][max2] && i != max1){
			max2 = i;
		}
	}
	maxs[0] = curve[0][max1];
	maxs[1] = curve[0][max2];
}

void get_max_3(double curve[2][200], int length, double maxs[3]){
	int i,max=0;
	for(i=1; i<length; i++){
		if(curve[1][i] > curve[1][max]){
			max = i;
		}
	}
	maxs[0] = curve[0][max-1];
	maxs[1] = curve[0][max];
	maxs[2] = curve[0][max+1];
}

void get_min_2(double curve[2][200], int length, double mins[2]){
	int min1,min2;
	min1 = 0; min2=0;
	int i;
	for(i=1; i<length; i++){
		if(curve[1][i] < curve[1][min1]){
			min1 = i;
		}
	}
	if(min1 == 0){
		min2 = 1;
	}
	for(i=1; i<length; i++){
		if(curve[1][i] < curve[1][min2] && i != min1){
			min2 = i;
		}
	}
	mins[0] = curve[0][min1];
	mins[1] = curve[0][min2];
}

void get_min_3(double curve[2][200], int length, double mins[3]){
	int i,min=0;
	for(i=1; i<length; i++){
		if(curve[1][i] < curve[1][min]){
			min = i;
		}
	}
	mins[0] = curve[0][min-1];
	mins[1] = curve[0][min];
	mins[2] = curve[0][min+1];
}

#endif
