#include<string>
#include<ostream>
#include<sstream>
#include<set>

class AluCall {
	public:
		std::string call;
		mutable int num;
		AluCall();
		AluCall(std::string call,int num);
		bool operator<(const AluCall& that) const;
		bool rsort_num(const AluCall& that) const;
		void add(const AluCall& that) const;
};

class AluCallSet {
	public:
		std::set<AluCall> calls;

		AluCallSet();
		void insert(const std::string& call, const int& score);
		void insert(const AluCall& call);
		void insert(const AluCallSet& call_set);
		void clear();
		std::set<AluCall>::iterator begin() const;
		std::set<AluCall>::iterator end() const;
};

std::ostream& operator<<(std::ostream& os, const AluCall call);
std::ostream& operator<<(std::ostream& os, const AluCallSet calls);

template <class T> std::string flatten(T& list, std::string delim="\t") {
	std::ostringstream test;
	for (typename T::iterator i=list.begin(); i!=list.end();i++) {
		test << *i << delim;
	}
	std::string t = test.str();
	return t.substr(0,t.length()-1);
}
