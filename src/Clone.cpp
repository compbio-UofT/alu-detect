#include "Clone.hpp"

#include <iostream>
#include <algorithm>

#include "Mapping.hpp"
#include "globals.hpp"
#include "deep_size.hpp"


template <class T>
Interval<T>
minEnclosingInterval(const Interval<T>& i1, const Interval<T>& i2)
{
  Interval<T> result;
  result[0] = min(i1[0], i2[0]);
  result[1] = max(i1[1], i2[1]);
  return result;
}

template<class T>
Interval<T>
minEnclosingInterval(vector<Interval<T> >& v)
{
  Interval<T> result;
  result[0] = v[0][0];
  result[1] = v[0][1];
  for_iterable(typename vector<Interval<T> >, v, it) {
    result[0] = min(result[0], (*it)[0]);
    result[1] = max(result[1], (*it)[1]);
  }
  return result;
}

Interval<long long>
minEnclosingInterval(vector<BP>& v)
{
  vector<Interval<long long> > u;
  for_iterable(vector<BP>, v, it) {
    u.push_back(it->pos);
  }
  return minEnclosingInterval(u);
}

BP
reduceBP(BP& bp1, BP& bp2)
{
  BP result;

  if (bp1.solid and bp2.solid) {
    result.solid = true;
    result.refAnchor[0] = bp1.refAnchor[0] or bp2.refAnchor[0];
    result.refAnchor[1] = bp1.refAnchor[1] or bp2.refAnchor[1];
    result.pos = minEnclosingInterval(bp1.pos, bp2.pos);
  } else if (bp1.solid or bp2.solid) {
    if (bp1.solid) {
      result = bp1;
    } else {
      result = bp2;
    }
  } else {
    result.solid = false;
    result.refAnchor[0] = bp1.refAnchor[0] or bp2.refAnchor[0];
    result.refAnchor[1] = bp1.refAnchor[1] or bp2.refAnchor[1];
    result.pos = minEnclosingInterval(bp1.pos, bp2.pos);
  }

  return result;
}

void
reduceIntersectingBP(vector<BP>& v)
{
  for (size_t i = 0; i < v.size(); ++i) {
    size_t j = i + 1;
    while (j < v.size()) {
      if (v[i].pos.intersects(v[j].pos)) {
	v[i] = reduceBP(v[i], v[j]);
	v.erase(v.begin() + j);
      } else {
	++j;
      }
    }
  }
}


Interval<long long int>
getFragPosFromMapping(const Mapping& mapping)
{
  Interval<long long int> result;

  result[0] = max(0ll, mapping.dbPos[0] - mapping.qrPos[0] - 100);
  result[1] = min(mapping.db->len - 1,
		  mapping.dbPos[1] + (mapping.qr->len - 1 - mapping.qrPos[1]) + 100);

  //cerr << "getFragPosFromMapping(" << mapping << ") = " << result << endl;

  return result;
}

vector<BP>
getBPListFromMapping(const Read& r, int refMappingIdx)
{
  const Mapping& mapping = r.mapping[refMappingIdx];
  vector<BP> result;
  BP tmp;

  if (mapping.qrPos[0] > 0) {
    tmp.pos[0] = max(0ll, mapping.dbPos[0] - mapping.qrPos[0] - 50);
    tmp.pos[1] = min(mapping.db->len - 1, mapping.dbPos[0] + 50);
    tmp.solid = (mapping.mqv >= 4 and refMappingIdx > 0);
    // reverse solidity preference because of tsds
    tmp.refAnchor[0] = true;
    tmp.refAnchor[1] = false;
    result.push_back(tmp);
  }

  if (mapping.qrPos[1] < mapping.qr->len - 1) {
    tmp.pos[0] = max(0ll, mapping.dbPos[1] - 50);
    tmp.pos[1] = min(mapping.db->len - 1,
		     mapping.dbPos[1] + (mapping.qr->len - 1 - mapping.qrPos[1]) + 50);
    tmp.solid = (mapping.mqv >= 4 and refMappingIdx < int(r.mapping.size()) - 1);
    tmp.refAnchor[0] = false;
    tmp.refAnchor[1] = true;
    result.push_back(tmp);
  }

  /*
  cerr << "getBPPosListFromMapping(" << mapping << ") = [";
  if (result.size() > 0)
    cerr << result[0];
  if (result.size() > 1)
    cerr << result[1];
  cerr << "]" << endl;
  */

  return result;
}

BP
getBPFromMissingMapping(const Mapping& mapping, int r_st, const Pairing& pairing)
{
  BP result;
  vector<Interval<long long int> > mp_pos = pairing.get_mp_pos(mapping, r_st);

  result.pos[0] = max(0ll, mp_pos[0][0] - 50);
  result.pos[0] = min(mapping.db->len - 1, result.pos[0]);
  result.pos[1] = min(mapping.db->len - 1, mp_pos[1][1] + 50);
  result.pos[1] = max(0ll, result.pos[1]);
  if (pairing.is_mp_downstream(mapping, r_st)) {
    result.pos[0] = max(result.pos[0], mapping.dbPos[1] - 10);
  } else {
    result.pos[1] = min(result.pos[1], mapping.dbPos[0] + 10);
  }
  result.solid = false;
  result.refAnchor[0] = false;
  result.refAnchor[1] = false;

  //cerr << "getBPPosFromMissingMapping(" << mapping << ") = " << result << endl;

  return result;
}

void
Clone::computePosition(bool remove_discordant, bool collapse_bp)
{
  int nip;

  //cerr << "compute position on: " << *this << endl;

  if (pairing == NULL) {
    cerr << "error: pairing unset" << endl;
    exit(1);
  }

  int refMappingIdx[2];
  refMappingIdx[0] = read[0].getRefMappingIdx();
  refMappingIdx[1] = read[1].getRefMappingIdx();

  if (!read[0].usable() and !read[1].usable())
    {
      ref = NULL;
    }
  else if (!read[0].usable() or !read[1].usable())
    {
      nip = (read[0].usable()? 0 : 1);
      if (refMappingIdx[nip] < 0) {
	ref = NULL;
      } else {
	ref = read[nip].mapping[refMappingIdx[nip]].db;
	if (read[nip].mapping[refMappingIdx[nip]].st == 1) {
	  read[nip].flip_st();
	}
	fragPos = getFragPosFromMapping(read[nip].mapping[refMappingIdx[nip]]);
	bp = getBPListFromMapping(read[nip], refMappingIdx[nip]);
	mappedToRepeatSt[0] = read[nip].mappedToRepeatSt[0];
	mappedToRepeatSt[1] = read[nip].mappedToRepeatSt[1];
      }
    }
  // now both usable
  else if (refMappingIdx[0] < 0 and refMappingIdx[1] < 0)
    {
      ref = NULL;
    }
  else
    {
      // if both reads are mapped to ref but discordant, remove mapping with lower mqv
      if (remove_discordant and refMappingIdx[0] >= 0 and refMappingIdx[1] >= 0
	  and !pairing->pair_concordant(read[0].mapping[refMappingIdx[0]], read[0].st,
					read[1].mapping[refMappingIdx[1]], read[1].st)) {
	if (read[0].mapping[refMappingIdx[0]].mqv < read[1].mapping[refMappingIdx[1]].mqv) {
	  read[0].mapping.erase(read[0].mapping.begin() + refMappingIdx[0]);
	  refMappingIdx[0] = read[0].getRefMappingIdx();
	} else {
	  read[1].mapping.erase(read[1].mapping.begin() + refMappingIdx[1]);
	  refMappingIdx[1] = read[1].getRefMappingIdx();
	}
      }

      if (refMappingIdx[0] >= 0 and refMappingIdx[1] >= 0)
	{
	  // now we know the pair is concordant
	  ref = read[0].mapping[refMappingIdx[0]].db;
	  if (read[0].mapping[refMappingIdx[0]].st == 1) {
	    read[0].flip_st();
	  }
	  if (read[1].mapping[refMappingIdx[1]].st == 1) {
	    read[1].flip_st();
	  }
	  fragPos = minEnclosingInterval(getFragPosFromMapping(read[0].mapping[refMappingIdx[0]]),
					 getFragPosFromMapping(read[1].mapping[refMappingIdx[1]]));
	  bp = getBPListFromMapping(read[0], refMappingIdx[0]);
	  vector<BP> tmp = 
	    getBPListFromMapping(read[1], refMappingIdx[1]);
	  bp.insert(bp.end(), tmp.begin(), tmp.end());
	}
      else
	{
	  // only one side mapped
	  nip = (refMappingIdx[0] >= 0? 0 : 1);
	  ref = read[nip].mapping[refMappingIdx[nip]].db;
	  if (read[nip].mapping[refMappingIdx[nip]].st == 1) {
	    read[nip].flip_st();
	  }
	  if (pairing->get_mp_st(read[nip].mapping[refMappingIdx[nip]], read[nip].st)
	      != read[1 - nip].st) {
	    read[1 - nip].flip_st();
	  }
	  bp = getBPListFromMapping(read[nip], refMappingIdx[nip]);
	  BP tmp =
	    getBPFromMissingMapping(read[nip].mapping[refMappingIdx[nip]], read[nip].st, *pairing);
	  bp.push_back(tmp);
	  Interval<long long int> tmp2 = minEnclosingInterval(bp);
	  //cerr << "tmp2=" << tmp2 << endl;
	  fragPos = minEnclosingInterval(getFragPosFromMapping(read[nip].mapping[refMappingIdx[nip]]),
					 tmp2);
	}
      mappedToRepeatSt[0] = read[0].mappedToRepeatSt[0] || read[1].mappedToRepeatSt[0];
      mappedToRepeatSt[1] = read[0].mappedToRepeatSt[1] || read[1].mappedToRepeatSt[1];
    }

  reduceIntersectingBP(bp);

  if (collapse_bp) {
    while (bp.size() > 1) {
      bp[0] = reduceBP(bp[0], bp[1]);
      bp.erase(bp.begin() + 1);
    }
  }

  size_t i = 0;
  while (i < bp.size()) {
    if (bp[i].pos[0] >= bp[i].pos[1]) {
      bp.erase(bp.begin() + i);
    } else {
      ++i;
    }
  }
}


ostream&
operator <<(ostream& ostr, const Clone& clone)
{
  ostr << "[clone_name=" << clone.name << ", ";
  if (clone.ref != NULL) {
    ostr << "ref=" << clone.ref->name << ", fragpos=" << clone.fragPos
	 << ", mappedToRepeatSt=[" << clone.mappedToRepeatSt[0] << "," << clone.mappedToRepeatSt[1]
	 << "], bpPos=[";
    if (clone.bp.size() > 0) {
      ostr << clone.bp[0].pos << ",solid=" << clone.bp[0].solid
	   << ",left=" << clone.bp[0].refAnchor[0] << ",right=" << clone.bp[0].refAnchor[1];
    }
    ostr << "], ";
  }
  ostr << clone.read[0] << ", " << clone.read[1] << "]";
  return ostr;
}

pair<int,Clone>*
readFastq(istream& istr, void (*fullNameParser)(const string&, Clone&, int&))
{
  int nip;
  Clone c;
  string s;

  // name
  getline(istr, s);
  if (istr.bad()) {
    cerr << "error reading Fastq entry" << endl;
    exit(1);
  }
  if (istr.eof()) {
    return NULL;
  }
  if (s[0] != '@') {
    cerr << "error: fastq read name does not being with '@'" << endl;
    exit(1);
  }
  fullNameParser(s.substr(1), c, nip);

  // seq
  getline(istr, s);
  if (istr.bad() or istr.eof()) {
    cerr << "error: could not read sequence from Fastq file" << endl;
    exit(1);
  }
  if (c.read[nip].seq == "") {
    c.read[nip].seq = s;
    if ((int)c.read[nip].seq.length() != c.read[nip].len) {
      cerr << "error: conflicting read length for read [" << c.read[nip].name << "]" << endl;
      exit(1);
    }
  } else if (c.read[nip].seq != s) {
    cerr << "error: conflicting sequences for read [" << c.read[nip].name << "]" << endl;
    exit(1);
  }

  // comment
  getline(istr, s);
  if (istr.bad() or istr.eof() or s[0] != '+') {
    cerr << "error: could not read qvString from Fastq file" << endl;
    exit(1);
  }
  // just ignore

  // qvString
  getline(istr, s);
  if (istr.bad() or istr.eof()) {
    cerr << "error: could not read qvString from Fastq file" << endl;
    exit(1);
  }
  if (c.read[nip].qvString == "") {
    c.read[nip].qvString = s;
    if ((int)c.read[nip].qvString.length() != c.read[nip].len) {
      cerr << "error: conflicting qvString length for read [" << c.read[nip].name << "]" << endl;
      exit(1);
    }
  } else if (c.read[nip].qvString != s) {
    cerr << "error: conflicting qvStrings for read [" << c.read[nip].name << "]" << endl;
    exit(1);
  }

  //cerr << "got read: " << c << endl;

  return new pair<int,Clone>(nip, c);
}


vector<Clone>
readAllFastq(istream& istr, 
	     string (*cloneNameParser)(const string&),
	     void (*fullNameParser)(const string&, Clone&, int&))
{
  vector<Clone> result;
  pair<int,Clone>* p[2];

  while (true) {
    p[0] = readFastq(istr, fullNameParser);
    if (p[0] == NULL) {
      break;
    }
    p[0]->second.name = cloneNameParser(p[0]->second.read[p[0]->first].name);
    p[1] = readFastq(istr, fullNameParser);
    if (p[1] == NULL) {
      cerr << "error: unmatched read entry in Fastq file" << endl;
      exit(1);
    }
    p[1]->second.name = cloneNameParser(p[1]->second.read[p[1]->first].name);
    if (p[0]->second.name != p[1]->second.name) {
      cerr << "error: consecutive entries in Fastq file mention different clones: "
	   << p[0]->second.name << " vs " << p[1]->second.name << endl;
      exit(1);
    }
    if (p[0]->first == p[1]->first) {
      cerr << "error: consecutive entries in Fastq file about the same clone side: "
	   << p[0]->second.name << ":" << p[0]->first << endl;
      exit(1);
    }
    // copy the read sequences and qvstrings
    Clone& c = p[0]->second;
    int& nip = p[1]->first;
    c.read[nip].st = p[1]->second.read[nip].st;
    c.read[nip].name = p[1]->second.read[nip].name;
    c.read[nip].seq = p[1]->second.read[nip].seq;
    c.read[nip].qvString = p[1]->second.read[nip].qvString;
    cerr << "got clone: " << c << endl;
    result.push_back(c);
    delete p[0];
    delete p[1];
  }

  return result;
}


size_t
size_below(const Clone& c)
{
  return size_below(c.name)
    + size_below(c.read[0])
    + size_below(c.read[1]);
}
