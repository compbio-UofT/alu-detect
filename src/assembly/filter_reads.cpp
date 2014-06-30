#include<iostream>
#include<cstdio>
#include<string>
#include<cstring>
#define MAX 1000
using  namespace std;

char sample[MAX],name[MAX],rest[MAX];
int num,j;
bool rfail;

inline bool readlines(char *lines) {
	for (j=0;j<4;j++) {
		rfail = fgets(lines + j*MAX,MAX,stdin) == NULL && !feof(stdin);
		if (rfail) { //if nothing read
			break;
		}
	}
	if (rfail) {
		cerr << j << "\n";
		throw 1;
	}
	return !feof(stdin);
}

int main(int argc, char* argv[]) {
	if (argc > 2) {
		cerr << "Only takes one argument, the prefix of the output files\n";
		return 1;
	}
	string NAME=argv[1];
	char lines[MAX * 4];

	FILE *out[2]={fopen((NAME + "_1.fq").c_str(),"w"),fopen((NAME + "_2.fq").c_str(),"w")};
	int i[2]={-1,-1};
	while (readlines(lines)) {
		sscanf(lines,"%[^':']:%[^':']:%d:%s",sample,name,&num,rest);
		i[--num]++;
		//append new read IDs
		fprintf(out[num],"%s:%011d,%s:%d:%s\n%s%s%s",sample,i[num],name,num+1,rest,lines + MAX,lines + 2*MAX,lines + 3*MAX);
/*		fputs(lines + MAX,out[num]);
		fputs(lines + 2*MAX,out[num]);
		fputs(lines + 3*MAX,out[num]);
		memmove(lines + MAX + strlen(lines + MAX),lines + 2*MAX,strlen(lines + 2*MAX) + 1);
		memmove(lines + MAX + strlen(lines + MAX),lines + 3*MAX,strlen(lines + 3*MAX) + 1);
		fputs(lines + MAX,out[num]);*/
		
		//reset input
		lines[0]=0;
	}
	fclose(out[0]);
	fclose(out[1]);
	return 0;
}
