#include "Call.h"

#include<sstream>
#include<cstdio>
#include<cstdlib>
#include<map>

#define DIFF 10 //if on polyA
#define SDIFF 2 //if not on polyA, stricter diff
#define OVERLAP 0.3

template <typename T> T max(const T& a, const T& b) {
	return a>b ? a : b;
}

template <typename T> T min(const T& a, const T& b) {
	return a>b ? b : a;
}
template <class T> void set_union(T& a, const T& b) {
	for (typename T::iterator it=b.begin();it!=b.end();it++) {
		a.insert(*it);
	}
}
Call::Call(){
	COORD[1]=NUM=-1;
	toprint=1;
}

int Call::get_sig() {
	return !COUNTS[1] + (!COUNTS[2] << 1);
}

Call::Call(std::string& line,std::string f,int lnum) {
	std::istringstream fields_s;
	int ind,i;
	char temp[10];
	toprint=1;
	fields_s.clear();
	fields_s.str(line);
	std::string fields[15],s;
	
	for (i=0;getline(fields_s,s,'\t');i++) {
		fields[i] = s;
	}
	if (i==15) {
		fields_s.clear();
		fields_s.str(fields[13]);
		fields_s >> this->NUM;
		fields[13] = fields[14];
	} else {
		this->NUM = 1;
	}

	ind = f.rfind('/');
	ind = ind == -1 ? 0 : ind + 1;
	read_suffix = f.substr(ind);
	ind = read_suffix.find('.');
	ind = ind == -1 ? read_suffix.length() : ind;

	sprintf(temp,"%d",lnum);
	sample = read_suffix.substr(0,ind);
	read_suffix = sample + ":" + std::string(temp);

	fields_s.clear();
	fields_s.str(fields[13]);
	while (getline(fields_s,s,',')) {
		READS.insert(s + ":" + read_suffix);
	}
	if (READS.empty()) {
		READS.insert("0:" + read_suffix);
	}
	COORD[3] = fields[5] == "+" ? 1 : 0;
	COORD.chr = fields[0];

	sscanf((fields[1]+" "+
			fields[2]+" "+
			fields[4]+" "+
			fields[8]+" "+
			fields[9]+" "+
			fields[10]+" "+
			fields[11]).c_str(),
			"%d%d%d%d%d%d%d",
			&COORD[1],
			&COORD[2],
			&score,
			&COUNTS[0],
			&COUNTS[1],
			&COUNTS[2],
			&COORD[4]
			);
	fields_s.clear();
	fields_s.str(fields[3]);
	while (getline(fields_s,s,',')) {
		CALLS.insert(AluCall(s,score));
	}
	evidence[get_sig()].insert(*this);
}

bool Call::less(const Call& that) const {
	if (this->COORD < that.COORD) {
		return 1;
	} else if (that.COORD < this->COORD) {
		return 0;
	} else {
		for (int i=1;i<=2;i++) {//check left coord, then right
			if (this->COUNTS[3-i] && that.COUNTS[3-i]) {//if has evidence
				if ((this->COORD.get(3) == 2-i && abs(this->COORD.get(i) - that.COORD.get(i)) < DIFF) || abs(this->COORD.get(i) - that.COORD.get(i)) < SDIFF) {
					continue;
				}
				if (this->COORD.get(i) < that.COORD.get(i)) {
					return 1;
				} else if (that.COORD.get(i) < this->COORD.get(i)) {
					return 0;
				}
			}
		}
//		int overlap = min(this->COORD.get(2),that.COORD.get(2)) - max(this->COORD.get(1),that.COORD.get(1)) + 1;
//		return ((float)overlap / (float)(that.COORD.get(2) - that.COORD.get(1)+1) < OVERLAP && (float)overlap / (float)(this->COORD.get(2) - this->COORD.get(1)+1) < OVERLAP); 
		return 0;
	}
}

bool Call::operator<(const Call& that) const {
	if(!this->less(that)) {
		return 0;
	}
	std::set<Call>::iterator it,jt;
	for (int i=0;i<3;i++) {
		for (it=this->evidence[i].begin();it!=this->evidence[i].end();it++) {
			if (!it->less(that)) {
				return 0;
			}
			for (jt=that.evidence[i].begin();jt!=that.evidence[i].end();jt++) {
				if (!it->less(*jt)) {
					return 0;
				}
			}
		}
		if (this->evidence[i].empty()) {
			for (jt=that.evidence[i].begin();jt!=that.evidence[i].end();jt++) {
				if (!this->less(*jt)) {
					return 0;
				}
				for (it=this->evidence[i].begin();it!=this->evidence[i].end();it++) {
					if (!it->less(*jt)) {
						return 0;
					}
				}
			}
		}
	}
	return 1;
}

bool Call::operator==(const Call& that) const {
	return !(*this < that) && !(that < *this);
}

void Call::update_info(const Call& B) {
	COUNTS[0] += B.COUNTS[0];
	COUNTS[1] += B.COUNTS[1];
	COUNTS[2] += B.COUNTS[2];
	NUM += B.NUM;
	CALLS.insert(B.CALLS);
	set_union(READS,B.READS);
	score += B.score;
}

//a is in tree, b is from file
Call Call::merge(const Call& a, const Call& b) {
	bool check = ((!a.COUNTS[1] + !a.COUNTS[2]) <= (!b.COUNTS[1] + !b.COUNTS[2]));
	Call A=check ? a : b,B=check ? b : a; //give preference to those with both or at least some breakpoint detected
	if (B.COUNTS[1] || !A.COUNTS[1]) {
		A.COORD[2] = max(A.COORD[2],B.COORD[2]);
	}
	if (B.COUNTS[2] || !A.COUNTS[2]) {
		A.COORD[1] = min(A.COORD[1],B.COORD[1]);
	}
	//add to the pile of evidence and think about it later, flattening the evidence
	for (int i=0;i<4;i++) {
		set_union(A.evidence[i],B.evidence[i]);
		B.evidence[i].clear();
	}
//	A.evidence[!B.COUNTS[1] + (!B.COUNTS[2] << 1)].insert(B);
	A.evidence[B.get_sig()].insert(B);
	if (check) {
		b.toprint=0;
	} else {
		a.toprint=0;
	}
	return A;
}

Call Call::flatten() const {
	Call call = *this;
	std::set<Call>::iterator it;

	call.READS.clear();
	call.CALLS.clear();
	std::map<unsigned int,unsigned int> coord[3];
	call.NUM = 0;
	bool single_coord = !call.COUNTS[1] || !call.COUNTS[2] || call.COORD[4];
	call.COUNTS[0]=0;
	call.COUNTS[1]=0;
	call.COUNTS[2]=0;
	call.score=0;

	unsigned int last,n,coords[] = {0,call.COORD.get(1),call.COORD.get(2)},scores[]={0,0,0};
	for (int j=0;j<4;j++) {
		for (it=call.evidence[j].begin();it!=call.evidence[j].end();it++) {
			for (int i=1;i<=2;i++) {
				if (it->COUNTS[3-i] && (!j || i==j)) {
					single_coord &= (!it->COUNTS[1] || !it->COUNTS[2] || it->COORD.get(4));
					last = coord[i].size();
					n = it->COORD.get(i);
					coord[i][n];
					if (coord[i].size() > last) {
						coord[i][n] = it->score;
					} else {
						coord[i][n] += it->score;
					}
					if (coord[i][n] > scores[i]) {
						coords[i] = n;
						scores[i] = coord[i][n];
					} else if (coord[i][n] == scores[i]) {
						if (i==1) {
							coords[i] = max(coords[i],n);
						} else {
							coords[i] = min(coords[i],n);
						}
					}
				}
			}
			call.NUM += it->NUM;
			call.COUNTS[0] += it->COUNTS[0];
			call.COUNTS[1] += it->COUNTS[1];
			call.COUNTS[2] += it->COUNTS[2];
			set_union(call.READS,it->READS);
			call.CALLS.insert(it->CALLS);
			call.score += it->score;
		}
	}
	if (call.COORD[4] >= 0) {
		if (coord[1].size()-1 || coord[2].size()-1) {
			call.COORD[4] = 0;
		} else {
			call.COORD[4] = (single_coord && call.COUNTS[1] && call.COUNTS[2] ? coords[2] - coords[1] - 2 : 0);
		}
	} else {
		call.COORD[4] = coords[1] - coords[2] + 2;
	}
	call.COORD[1] = coords[1];
	call.COORD[2] = coords[2];
	return call;
}

std::ostream& operator<<(std::ostream& os, const Call& cll) {
	Call call = cll.flatten();

	os	<< call.COORD.chr << "\t"
		<< call.COORD[1] << "\t"
		<< call.COORD[2] << "\t"
		<< call.CALLS << "\t"
		<< call.score/call.NUM << "\t"
		<< (call.COORD[3] ? "+" : "-") << "\t"
		<< ".\t"
		<< ".\t"
		<< call.COUNTS[0] << "\t"
		<< call.COUNTS[1] << "\t"
		<< call.COUNTS[2] << "\t"
		<< call.COORD[4] << "\t"
		<< ".\t"
		<< call.NUM;
	if (call.READS.size()) {
		os << "\t" << flatten(call.READS,",");
	}
//	os << "\n";

	return os;
}

