#include<iostream>
#include<set>
#include<string>

#include "Coord.h"
//#include "Evidence.h"
#include "AluCall.h"


class Call {
	private:
		std::string read_suffix;
		struct CallCmp {
			bool operator()(const Call& a, const Call& b) {
				return a.read_suffix < b.read_suffix;
			}
		};
	public:
		AluCallSet CALLS;
		mutable std::set<Call,CallCmp> evidence[4]; //+1 no left bp, +2 no right bp
		int NUM;
		Coord COORD;
		int COUNTS[3];
		int score;
		mutable bool toprint;
		std::string sample;

		Call();
		Call(std::string& line, std::string f, int lnum);
		mutable std::set<std::string> READS;
		bool less(const Call& that) const;
		bool operator<(const Call& that) const;
		bool operator==(const Call& that) const;
		void update_info(const Call& B);
		void add_evidence(const Call& b);
		int get_sig();
		static Call merge(const Call& a, const Call& b);
		Call flatten() const;
};


std::ostream& operator<<(std::ostream& os, const Call& call);

