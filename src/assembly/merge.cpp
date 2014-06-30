#include<iostream>
#include<string>
#include<cstring>
#include<fstream>
#include<utility>

#include "CallTree.h"

int main(int argc, char* argv[]) {
	std::ifstream calls_file;
	CallTree call_tree[16]; 
	Call call;
	std::string line;
	int i=1,j=0,str;
	int debug = (!strcmp(argv[1],"-d") || !strcmp(argv[1],"--debug")) + ((!strcmp(argv[1],"-D") || !strcmp(argv[1],"--debug2")) << 1);
	std::ostream& out = debug ? i++,std::cout : std::cout;
	
	std::set<Call>::iterator it,jt;

	for (;i<argc;i++) {
		calls_file.open(argv[i]);
		while (getline(calls_file,line)) {
			//merge everything with same support kind first
			call = Call(line,argv[i],++j);
			
			//+8: TS deletion, +4: negative, +2: no left coord, +1: no right coord
			str = ((call.COORD.get(4)<0) << 3) + ((!call.COORD.get(3)) << 2) + ((!call.COUNTS[2]) << 1) + (!call.COUNTS[1]);

			call_tree[str].add(call,1);

		}
		calls_file.close();
		j=0;
	}

	//now merge everything into everything
	std::pair<std::set<Call>::iterator,bool> insert;
	for (i=0;i<13;i+=4) {
		//merge 0bp into everything
		std::set<Call> temp;
		if (!debug) {
			for (j=i;j<i+3;j++) {
				call_tree[j].merge_calls(call_tree[i+3],0,0);
			}
			for (it=call_tree[i+2].call_tree.begin();it!=call_tree[i+2].call_tree.end();it++) {
				temp.clear();
				while (jt = call_tree[i+1].call_tree.find(*it), jt != call_tree[i+1].call_tree.end()) {
					call_tree[i].add(Call::merge(*it,*jt),1);
					temp.insert(*jt);
					call_tree[i+1].call_tree.erase(jt);
					it->toprint=0;
					jt->toprint=0;
				}
				for (jt=temp.begin();jt!=temp.end();jt++) {
					call_tree[i+1].add(*jt,1);
				}
			}
		}

		for (j=i;j<i+4;j++) {
			out << call_tree[j];
		}
	}
	return 0;
}

