#include<set>
#include<ostream>

#include "Call.h"

class CallTree {
	public:
		std::set<Call> call_tree;
		
		CallTree();
//		std::set<Call>::iterator add(const Call& cll, bool mrg);
		bool add(const Call& cll, bool mrg);
		void merge_calls(CallTree& that_tree, bool toprint, bool todelete);
};

std::ostream& operator<<(std::ostream& os, const CallTree& tree);
