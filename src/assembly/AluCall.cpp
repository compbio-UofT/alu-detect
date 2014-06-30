#include "AluCall.h"

AluCall::AluCall() {
	call="";
	num=-1;
}

AluCall::AluCall(std::string call, int num) {
	this->call = call;
	this->num = num;
}

bool AluCall::operator<(const AluCall& that) const {
	return call < that.call;
}

std::ostream& operator<<(std::ostream& os, const AluCall call) {
	os << call.call;
	return os;
}

void AluCall::add(const AluCall& that) const {
	if (that.call == call) {
		this->num += that.num;
	}
}

AluCallSet::AluCallSet() {
}

void AluCallSet::insert(const AluCall& call) {
	std::pair<std::set<AluCall>::iterator,bool> ins = calls.insert(call);
	if (!ins.second) {
		ins.first->add(call);
	}
}

void AluCallSet::insert(const std::string& call, const int& score) {
	insert(AluCall(call,score));
}

void AluCallSet::insert(const AluCallSet& call_set) {
	for (std::set<AluCall>::iterator it=call_set.begin();it!=call_set.end();it++) {
		insert(*it);
	}
}

void AluCallSet::clear() {
	calls.clear();
}

std::set<AluCall>::iterator AluCallSet::begin() const {
	return calls.begin();
}

std::set<AluCall>::iterator AluCallSet::end() const {
	return calls.end();
}

std::ostream& operator<<(std::ostream& os, const AluCallSet calls) {
	AluCallSet best;
	std::set<AluCall>::iterator it = calls.begin();
	best.insert(*it++);
	for (;it!=calls.end();it++) {
		if (it->num >= best.begin()->num) {
			if (it->num > best.begin()->num) {
				best.clear();
			}
			best.insert(*it);
		}
	}
	os << flatten(best.calls,",");

	return os;
}
