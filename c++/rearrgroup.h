/**   LICENCE
* Copyright (c) 2014-2017 Genome Research Ltd.
*
* Author: Cancer Genome Project <cgpit@sanger.ac.uk>
*
* This file is part of BRASS.
*
* BRASS is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*
* 1. The usage of a range of years within a copyright statement contained within
* this distribution should be interpreted as being equivalent to a list of years
* including the first and last year specified and all consecutive years between
* them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
* 2009, 2011-2012’ should be interpreted as being identical to a statement that
* reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
* statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
* identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
* 2009, 2010, 2011, 2012’."
*/

// Author: John Marshall

// rearrgroup.h -- Classes for rearrangement groups and sets thereof.

#ifndef REARRGROUP_H
#define REARRGROUP_H

#include <algorithm>
#include <iterator>
#include <list>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"

using sam::coord_t;
using sam::scoord_t;
using sam::alignment;
using sam::collection;

// Information extracted from @RG headers

struct readgroup_info {
  std::string sample;
  int sample_index;
  bool is_main;
  scoord_t max_insert;

  readgroup_info() { }
  readgroup_info(const std::string& s, bool m, scoord_t mi)
    : sample(s), sample_index(0), is_main(m), max_insert(mi) { }
};

class readgroup_set {
public:
  readgroup_set(const collection& headers, scoord_t default_max,
		std::string default_sample, std::string main_sample);
  ~readgroup_set() { }

  const readgroup_info& find(const alignment& aln) const;

  const std::vector<std::string>& samples() const { return samples_; }

private:
  struct c_str_less {
    bool operator() (const char* a, const char* b) const
      { return strcmp(a, b) < 0; }
  };

  typedef std::map<const char*, readgroup_info, c_str_less> readgroup_map;
  readgroup_map readgroups;
  std::vector<std::string> samples_;

  static const char* const NO_RG;
};


struct interval {
  scoord_t pos5, pos3;

  interval() { }
  interval(scoord_t p5, scoord_t p3) : pos5(p5), pos3(p3) { }
  interval(const interval& rhs) : pos5(rhs.pos5), pos3(rhs.pos3) { }
  ~interval() { }

  interval& assign(scoord_t p5, scoord_t p3)
    { pos5 = p5; pos3 = p3; return *this; }

  interval(const alignment& aln, coord_t pos, int strand,
	   coord_t ref_length, const readgroup_info& info);

  scoord_t length() const { return pos3 - pos5; }

  // Assigns the result of set intersection with RHS.  If this and RHS
  // do not intersect, the result is an empty interval with pos5 > pos3.
  interval& operator*= (const interval& rhs) {
    if (rhs.pos5 > pos5)  pos5 = rhs.pos5;
    if (rhs.pos3 < pos3)  pos3 = rhs.pos3;
    return *this;
  }

  // Assigns an interval that covers the result of set union with RHS.
  // If this and RHS intersect, this is just union; otherwise it's more.
  interval& operator+= (const interval& rhs) {
    if (rhs.pos5 < pos5)  pos5 = rhs.pos5;
    if (rhs.pos3 > pos3)  pos3 = rhs.pos3;
    return *this;
  }
};

// Returns whether the set intersection of LHS and RHS is non-empty.
inline bool intersect(const interval& lhs, const interval& rhs) {
  return ! (lhs.pos3 < rhs.pos5 || lhs.pos5 > rhs.pos3);
}


class rearr_group {
public:
  rearr_group() { }
  rearr_group(alignment& aln, const interval& alnL, const interval& alnH,
	      const readgroup_info& info, const readgroup_set& readgroups);
  ~rearr_group() { }

  friend std::ostream& operator<< (std::ostream&, const rearr_group&);

  typedef std::map<std::string, int>::iterator mate_iterator;

  void insert(const alignment& aln, const interval& alnL, const interval& alnH,
	      const readgroup_info& info);
  void insert_mate(mate_iterator hint, const alignment& aln);

  // Reconstruct afresh: assign()+reinsert()xN is similar to construct+insert()xN
  // but keeps the existing alnlist, mate_primary_count, and higher_count.
  void assign(alignment& aln, const interval& alnL, const interval& alnH,
	      const readgroup_info& info, const readgroup_set& readgroups);
  void reinsert(const alignment& aln, const interval& alnL,
		const interval& alnH, const readgroup_info& info);

  // Update alnlist based on mate data from insert_mate() calls
  bool synchronise();

  bool
  matches(const alignment& aln, const interval& alnL, const interval& alnH) {
    return aln.strand() == canonical.strand() &&
	   aln.mate_strand() == canonical.mate_strand() &&
	   intersect(overlapL, alnL) && intersect(overlapH, alnH) &&
	   1 /* complicated thing involving $dir */;
  }

  bool matches_mate(const alignment& aln, const interval& alnL,
		    const interval& alnH, mate_iterator& hint) {
    static std::string buf;
    return aln.mate_strand() == canonical.strand() &&
	   aln.strand() == canonical.mate_strand() &&
	   intersect(overlapL, alnL) && intersect(overlapH, alnH) &&
	   1 /* complicated thing involving $dir */ &&
	   (hint = higher_count.find(aln.qname(buf))) != higher_count.end();
  }

  int rindex() const { return canonical.rindex(); }
  int mate_rindex() const { return canonical.mate_rindex(); }
  std::string rname() const { return canonical.rname(); }
  std::string mate_rname() const { return canonical.mate_rname(); }

  int lower_total_count() const { return total_count; }
  int lower_primary_count() const { return primary_count; }
  int higher_primary_count() const { return mate_primary_count; }
  int higher_total_count() const;
  int sample_count(int sindex) const { return samples[sindex].count; }
  const std::set<std::string>& sample_reads(int sindex) const
    { return samples[sindex].readnames; }

  bool either_side_stacked() const
    { return readL.length() <= max_read_length + 1 ||
	     readH.length() <= max_read_length + 1; }

  void set_notes(const std::string& str) { notes = str; }

  // Iterator and begin()/end() for iterating over alignments in the group
  typedef std::list<alignment>::iterator iterator;
  iterator begin() { return alnlist.begin(); }
  iterator end() { return alnlist.end(); }

//private: // FIXME  Decide on encapsulation of these fields
  alignment canonical;
  interval overlapL, overlapH, readL, readH;
  std::string notes;
  scoord_t max_insert;
private:
  int max_read_length;

  struct per_sample {
    int count;
    std::set<std::string> readnames;
    per_sample() : count(0) { }
  };
  std::vector<per_sample> samples;
  int total_count, primary_count, mate_primary_count;
  std::list<alignment> alnlist;
  std::map<std::string, int> higher_count;

  typedef std::vector<per_sample>::iterator sample_iterator;
  typedef std::vector<per_sample>::const_iterator const_sample_iterator;
};

std::ostream& operator<< (std::ostream& stream, const rearr_group& group);


class rearr_group_set {
public:
  typedef std::list<rearr_group>::iterator iterator;
  typedef std::pair<iterator, iterator> iterator_pair;

  rearr_group_set(const collection& refseqs)
    : ref_size(refseqs.ref_size()), lists(I(ref_size-1, ref_size-1) + 2),
      list_sentinel(lists.end() - 1) {
    // Vector sentinel needs to be a non-empty list
    list_sentinel->push_front(rearr_group());
  }
  ~rearr_group_set() { }

  // Fetch groups matching ALN (given ALN itself is less than its mate)
  iterator_pair range_lower(const alignment& aln)
    { std::list<rearr_group>& list = lists[I(aln.rindex(), aln.mate_rindex())];
      return make_pair(list.begin(), list.end()); }

  // Fetch groups matching ALN (given ALN's mate is less than ALN itself)
  iterator_pair range_higher(const alignment& aln)
    { std::list<rearr_group>& list = lists[I(aln.mate_rindex(), aln.rindex())];
      return make_pair(list.begin(), list.end()); }

  void insert(const rearr_group& group) { lists[I(group)].push_back(group); }
  iterator erase(iterator pos) { return lists[I(*pos)].erase(pos); }

  iterator transfer(std::list<rearr_group>& dest, iterator pos)
    { iterator succ = pos; ++succ;
      dest.splice(dest.end(), lists[I(*pos)], pos);
      return succ; }

  // Iterator and begin()/end() for iterating over all groups
  class full_iterator :
    public std::iterator<std::forward_iterator_tag, rearr_group> {
  public:
    typedef std::vector<std::list<rearr_group> > list_vector;

    full_iterator() { }
    full_iterator(const full_iterator& rhs) : i(rhs.i), it(rhs.it) { }
    ~full_iterator() { }
    full_iterator& operator= (const full_iterator& rhs)
      { i = rhs.i; it = rhs.it; return *this; }

    bool operator== (const full_iterator& rhs) const
      { return i == rhs.i && it == rhs.it; }
    bool operator!= (const full_iterator& rhs) const
      { return i != rhs.i || it != rhs.it; }

    rearr_group& operator* () const { return *it; }
    rearr_group_set::iterator operator-> () const { return it; }

    full_iterator& operator++ ()
      { ++it;
	while (it == i->end())  it = (++i)->begin();
	return *this; }

  private:
    friend class rearr_group_set;
    full_iterator(list_vector::iterator i0, rearr_group_set::iterator it0)
      : i(i0), it(it0) { while (it == i->end())  it = (++i)->begin(); }

    list_vector::iterator i;
    rearr_group_set::iterator it;
  };

  full_iterator begin()
    { return full_iterator(lists.begin(), lists.front().begin()); }
  full_iterator end()
    { return full_iterator(list_sentinel, list_sentinel->begin()); }

  full_iterator erase(const full_iterator& pos)
    { return full_iterator(pos.i, pos.i->erase(pos.it)); }

  full_iterator transfer(std::list<rearr_group>& dest, full_iterator pos)
    { full_iterator succ = pos; ++succ;
      dest.splice(dest.end(), *(pos.i), pos.it);
      return succ; }

private:
  // Compute 1D index into a 2D triangular half matrix (given i <= j)
  size_t I(size_t i, size_t j) const { return i*(2*ref_size - (i+1))/2 + j; }
  size_t I(const rearr_group& g) const { return I(g.rindex(),g.mate_rindex()); }

  const size_t ref_size;
  std::vector<std::list<rearr_group> > lists;
  full_iterator::list_vector::iterator list_sentinel;
};

#endif
