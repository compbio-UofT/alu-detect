#include<utility>

#include "CallTree.h"

CallTree::CallTree() {
}

bool CallTree::add(const Call& cll, bool mrg) { //used to add stuff
	Call tcall, call = cll;
	std::set<Call> temp;
	std::pair<std::set<Call>::iterator,bool> ins;
	std::set<Call>::iterator cur;
	cur = call_tree.find(call);
	if (cur != call_tree.end()) {
		//if found, adds it
		temp.clear();
		temp.insert(Call::merge(*cur,call));
		call_tree.erase(cur);
		//sees if the newly added thing overlaps with anything else
		while (cur = call_tree.find(call), cur!=call_tree.end()) {
			tcall = Call::merge(*cur,call);
			call_tree.erase(cur);
			//if it does, merge it with everything else found previously
			while (ins = temp.insert(tcall), !ins.second) {
				tcall = Call::merge(tcall,*ins.first);
				temp.erase(ins.first);
			}
		}
		//readds merged stuff to tree
		for (cur=temp.begin();cur!=temp.end();cur++) {
			tcall = *cur;
			while (ins = call_tree.insert(tcall), !ins.second) {
				tcall = Call::merge(tcall,*ins.first);
				call_tree.erase(ins.first);
			}
			call_tree.insert(tcall);
		}
		return 1;
	} else {
		if (mrg) {
			call_tree.insert(call);
		}
		return 0;
	}
}
void CallTree::merge_calls(CallTree& that_tree, bool toprint, bool todelete) {
	std::set<Call> calls;
	bool inserted;
	std::set<Call>::iterator it,jt;//,insert;
	std::pair<std::set<Call>::iterator, bool> insert_try;
	for (it=that_tree.call_tree.begin();it!=that_tree.call_tree.end();) {
		inserted = add(*it,0);
		if (inserted) {
			it->toprint=it->toprint && toprint;
			if (todelete) {
				that_tree.call_tree.erase(it++);
			} else {
				it++;
			}
		} else {
			it++;
		}
	}
}

std::ostream& operator<<(std::ostream& os, const CallTree& call_tree) {
	for (std::set<Call>::iterator cur=call_tree.call_tree.begin();cur!=call_tree.call_tree.end();cur++) {
		if (cur->toprint) {
			os << *cur << "\n";
		}
	}
	return os;
}
