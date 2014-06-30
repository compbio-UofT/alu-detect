#include<iostream>
#include<string>
#include<cstdio>
#include<cstring>
#include<sstream>
#include<fstream>
#include<map>
#define MAX 10000
using namespace std;

int main(int argc,char* argv[]) {
	//generate a map from read ID to read label
	map<string,string> filter;
	ifstream novel(argv[1]);
	string line,chr,first,second,calls,strand,reads,temp,label;
	ostringstream lab;
	istringstream reads_s;
	int i;
	for (i=1;novel >> chr >> first >> second >> calls >> temp >> strand >> temp >>  temp >> temp >> temp >> temp >> temp >> temp;i++) {
//		ostringstream lab;
		lab.clear();
		lab.str("");
		lab << ":" << argv[2] << ":" << i << ":" << chr << ":" << first << ":" << second << ":" << calls << ":" << strand;
		label = lab.str();
		novel >> reads;
		reads_s.clear();
		reads_s.str(reads);
//		istringstream reads_s(reads);
		while (getline(reads_s,temp,',')) {
			filter[temp] = label;
		}
	}
	novel.close();

	map<string,string>::iterator check;
	char id1[MAX],id2[MAX],l1[MAX],l2[MAX],l3[MAX],l4[MAX];

	string pname = "pv ";
	for (i=4;i<argc;i++) {
		pname += string(argv[i]) + " ";
	}
	pname += "| zcat";
	
	//read read files, output those which are in read ID map, add label
	FILE *outf = fopen((string(argv[3]) + ".reads." + string(argv[2]) + ".fq").c_str(),"w");
	FILE *inf = popen(pname.c_str(),"r");
	while (fgets(l1,MAX,inf) && fgets(l2,MAX,inf) && fgets(l3,MAX,inf) && fgets(l4,MAX,inf)) {
		sscanf(l1,"@%[^':']:%[^':']",id1,id2);
		l1[strlen(l1) - 1] = 0;

		check = filter.find(id2);

		if (check != filter.end()) {
			fprintf(outf,"%s%s\n%s%s%s",l1,check->second.c_str(),l2,l3,l4);
		} else {
			check = filter.find(id1);
			if (check != filter.end()) {
				fprintf(outf,"@0:%s%s\n%s%s%s",l1,check->second.c_str(),l2,l3,l4);
			}
		}
	}
	fclose(outf);
	pclose(inf);
	return 0;
}
