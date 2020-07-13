/**   LICENCE
* Copyright (c) 2014-2017 Genome Research Ltd.
* Portions copyright (C) 2018-2020 University of Glasgow.
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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cerrno>
#include <cstdlib>

#include <time.h>
#include <unistd.h>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "cansam/exception.h"
#include "cansam/intervalmap.h"

#include "feature.h"
#include "imergestream.h"
#include "rearrgroup.h"
#include "version.inc"

using std::string;
using namespace sam;


// Expands occurrences of "%XY" in TEXT into the value of REF's "XY" field.
void expand_ref(string& text, const header& ref) {
  size_t pos = 0;
  while ((pos = text.find('%', pos)) != string::npos && pos+2 < text.length()) {
    string value = ref.field<string>(&text[pos+1]);
    text.replace(pos, 3, value);
    pos += value.length();
  }
}

// Opens FILENAME in FILE and returns it, or returns std::cout when given "-".
std::ostream& open_or_cout(std::ofstream& file, const string& filename) {
  if (filename == "-")  return std::cout;

  file.open(filename.c_str());
  if (! file.is_open())
    throw sam::system_error("can't write to ", filename, errno);

  return file;
}

// Returns whether ALN itself is less than its mate, either by location or,
// in the (unlikely) case that both are mapped to the same position, by being
// the one that has the FIRST flag set.
inline bool less_than_mate(const alignment& aln) {
  if (aln.rindex() != aln.mate_rindex())
    return aln.rindex() < aln.mate_rindex();
  else if (aln.pos() != aln.mate_pos())
    return aln.pos() < aln.mate_pos();
  else
    return (aln.flags() & FIRST_IN_PAIR) != 0;
}

// Returns whether ALN and its mate are in the natural orientation for the
// library type.  FIXME  For us, this is really a property of the rearr_group.
inline bool natural_orientation(const alignment& aln) {
  // FIXME Not implemented for other than short-insert-solexa
  int natural = MATE_REVERSE_STRAND;

  return (aln.flags() & (REVERSE_STRAND | MATE_REVERSE_STRAND)) == natural;
}

// Returns whether ALN appears to span a small intrachromosomal insertion,
// i.e., the two reads are on the same chromosome in the natural orientation
// for the library type, with too short an insert.
inline bool
apparent_insertion(const alignment& aln, const readgroup_info& info) {
  // FIXME This is only right for short-insert-solexa
  int natural = less_than_mate(aln)? MATE_REVERSE_STRAND : REVERSE_STRAND;

  return aln.rindex() == aln.mate_rindex() &&
      (aln.flags() & (REVERSE_STRAND | MATE_REVERSE_STRAND)) == natural &&
      std::abs(aln.isize()) <= info.max_insert;
}

// Record and print out the time taken for various computations.
class stopwatch {
public:
  stopwatch(bool started = true) : total(0), starttime(started? now() : 0) { }
  void restart() { total = 0; starttime = now(); }
  void start() { starttime = now(); }
  void stop() { total += now() - starttime; starttime = 0; }

  unsigned long long elapsed() const
    { return starttime? total + now() - starttime : total; }

private:
  static unsigned long long now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000000000ULL + ts.tv_nsec;
  }

  unsigned long long total, starttime;
};

std::ostream& operator<< (std::ostream& s, const stopwatch& timer) {
  unsigned long long elapsed = timer.elapsed();
  char fill = s.fill();
  return s << (elapsed / 1000000000ULL) << '.'
	   << std::setfill('0') << std::setw(9) << (elapsed % 1000000000ULL)
	   << std::setfill(fill);
}

// Construct a seqinterval from a (rearrgroup.h) interval.
inline seqinterval
make_seqinterval(const string& rname, const ::interval& ival) {
  return seqinterval(rname, ival.pos5 - 1, ival.pos3);
}


// Command-line options.
struct options {
  string output_filename;
  std::vector<string> ignores;
  std::vector<string> ignore_filenames;
  std::vector<string> feature_filenames;
  std::vector<string> anchor_filenames;
  std::vector<string> retrotransposon_filenames;
  std::map<string, bool> discards;
  string default_sample;
  int max_insert;
  int min_count;
  int min_quality;
};


class rearrangement_grouper {
public:
  rearrangement_grouper(const options& opt, const collection& headers);
  ~rearrangement_grouper() { }

  template <typename InputSamStream>
  void group_alignments(InputSamStream& in);
  void recalculate_group(rearr_group& group);

  void print_preamble(const collection& headers, const string& preamble);
  void print_trailer() { out << "#\n"; print_statistics(out, "# ", "#\n"); }
  void log_statistics(std::ostream& s) { print_statistics(s, "", "\n"); }

private:
  readgroup_set readgroups;
  interval_multimap<feature> ignores;
  interval_multimap<feature> filters;
  interval_multimap<feature> transposons;
  interval_multimap<feature> anchors;
  interval_multimap<feature> retrotransposons;
  rearr_group_set active;
  interval_multimap<rearr_group*> active_by_readL;
  interval_multimap<rearr_group*> active_by_readH;
  std::vector<coord_t> ref_length;
  bool discard_apparent_insertions;
  bool discard_within_repeats;
  bool discard_repeat_mapped;
  int min_count;
  int min_quality;
  int min_clipped_length;
  double max_polytract_frac;
  double max_samebase_frac;

  std::ofstream outfile;
  std::ostream& out;

  struct {
    unsigned long total, proper, unmapped, low_quality, low_mate_quality,
		  clipped_short, polytract, samebase,
		  repeats, repetitive, ignored, insertion, near_mate;
    void clear()
      { total = proper = unmapped = low_quality = low_mate_quality =
	  clipped_short = polytract = samebase =
	  repeats = repetitive = ignored = insertion = near_mate = 0; }
  } read_stats;

  struct {
    unsigned long total, small, stacked, supponly, emitted, unanchored,
		  swamped;
    void clear()
      { total = small = stacked = supponly = emitted = unanchored =
	  swamped = 0; }
  } group_stats, pass1_group_stats;

  // Methods that update GROUP by adding notes and annotations.
  void annotate_intrachromosomal_deletions(rearr_group& group);

  // Returns true iff ALN (as an interval) is covered by intervals in FEATURES,
  // with up to MAX_UNCOVERED positions remaining uncovered.
  // FIXME Once interval_multimap has const iterators etc, revert to const
  bool within(interval_multimap<feature>& features, const seqinterval& aln,
	      int max_uncovered = 0);

  // Returns true iff ALN (as an interval) is covered by repeat features.
  // FIXME Once interval_multimap has const iterators etc, revert to const
  bool within_repeat(const seqinterval& aln)
    { return within(filters, aln, 10); }

  // Returns the number of distinct reads from other groups that intersect
  // with GROUP's read windows (either rname/readL or mate_rname/readH).
  int count_swamping_reads(rearr_group* group, const string& group_rname,
			   const ::interval& group_read_window);

  void print_statistics(std::ostream& s, const char* prefix, const char* blank);
};

rearrangement_grouper::rearrangement_grouper(const options& opt,
					     const collection& headers)
  : readgroups(headers, opt.max_insert, opt.default_sample),
    active(headers),
    ref_length(headers.ref_size()),
    discard_apparent_insertions(opt.discards.find("insertion")->second),
    discard_within_repeats(opt.discards.find("repeat")->second),
    discard_repeat_mapped(opt.discards.find("repetitive")->second),
    min_count(opt.min_count),
    min_quality(opt.min_quality),
    min_clipped_length(35),
    max_polytract_frac(0.5),
    max_samebase_frac(0.9),
    outfile(),
    out(open_or_cout(outfile, opt.output_filename)) {

  // Ensure that we hear about any write failures.
  out.exceptions(std::ios::failbit | std::ios::badbit);

  // We need to know, in advance, how many reference sequences there are.
  if (headers.ref_empty())
    throw sam::exception("input contains no reference sequence headers");

  const header& ref = *(headers.ref_begin());

  for (std::vector<string>::const_iterator it = opt.ignores.begin();
       it != opt.ignores.end(); ++it)
    ignores.insert(std::make_pair(seqinterval(*it), string()));

  for (std::vector<string>::const_iterator it = opt.ignore_filenames.begin();
       it != opt.ignore_filenames.end(); ++it) {
    string filename = *it;
    expand_ref(filename, ref);
    insert(ignores, ignores, ignores, filename, ignore_reads);
  }

  for (std::vector<string>::const_iterator it = opt.feature_filenames.begin();
       it != opt.feature_filenames.end(); ++it) {
    string filename = *it;
    expand_ref(filename, ref);
    if (insert(filters, transposons, ignores, filename) == 0) {
      sam::exception e("no features selected (missing track metadata?)");
      e.set_filename(filename);
      throw e;
    }
  }

  for (std::vector<string>::const_iterator it = opt.anchor_filenames.begin();
       it != opt.anchor_filenames.end(); ++it) {
    string filename = *it;
    expand_ref(filename, ref);
    insert(anchors, anchors, anchors, filename, filter_reads);
  }

  for (std::vector<string>::const_iterator it = opt.retrotransposon_filenames.begin();
       it != opt.retrotransposon_filenames.end(); ++it) {
    string filename = *it;
    expand_ref(filename, ref);
    insert(retrotransposons, retrotransposons, retrotransposons, filename, filter_reads);
  }

  for (size_t i = 0; i < headers.ref_size(); i++)
    ref_length[i] = headers.findseq(i).length();

  read_stats.clear();
  group_stats.clear();
  pass1_group_stats.clear();
}

void rearrangement_grouper::print_preamble(const collection& headers,
					   const string& preamble) {
  out << preamble;

  const refsequence& ref = *(headers.ref_begin());
  if (ref.find("SP") != ref.end()) {
    out << "#\n#REFERENCE\tSP:" << ref.species();
    if (ref.find("AS") != ref.end())  out << "\tAS:" << ref.assembly();
    if (ref.find("UR") != ref.end())  out << "\tUR:" << ref.uri();
    if (ref.find("M5") != ref.end())  out << "\tM5:" << ref.checksum();
    out << "\n#\n";
  }

  const std::vector<string>& samples = readgroups.samples();
  out << "#NSAMPLES\t" << samples.size() << '\n';
  for (size_t i = 0; i < samples.size(); i++)
    out << "#SAMPLE\t" << i+1 << '\t' << samples[i] << '\n';
}

void rearrangement_grouper::print_statistics(std::ostream& s, const char* p,
					     const char* blank) {
  s << p << "Total reads scanned:\t" << read_stats.total << '\n'
    << p << "Reads discarded due to being\n"
    << p << "  Properly paired:\t" << read_stats.proper << '\n'
    << p << "  (Half-)unmapped:\t" << read_stats.unmapped << '\n'
    << p << "  Near mate:\t\t" << read_stats.near_mate << '\n';

  if (min_quality > 0)
    s << p << "  Low quality:\t" << (*p? "" : "\t") << read_stats.low_quality << '\n'
      << p << "  Low mate quality:\t" << read_stats.low_mate_quality << '\n';
  if (discard_apparent_insertions)
    s << p << "  Small insertion:\t" << read_stats.insertion << '\n';
  if (discard_within_repeats)
    s << p << "  Repeat features:\t" << read_stats.repeats << '\n';
  if (discard_repeat_mapped)
    s << p << "  Repeat-mapped:\t" << read_stats.repetitive << '\n';
  if (! ignores.empty())
    s << p << "  In ignored regions:\t" << read_stats.ignored << '\n';

  s << p << "  Clipped to < " << min_clipped_length << "bp:\t"
	 << read_stats.clipped_short << '\n'
    << p << "  > " << max_polytract_frac * 100.0 << "% poly. tract:\t"
	 << read_stats.polytract << '\n'
    << p << "  > " << max_samebase_frac * 100.0 << "% same base:\t"
	 << read_stats.samebase << '\n';

  s << blank
    << p << "Total groups found:\t" << group_stats.total << '\n';

  s << p << "Rearrangement groups omitted due to being\n"
    << p << "  Supplementary-only:\t" << group_stats.supponly
	 << "\t  (first pass: " << pass1_group_stats.supponly << ")\n";
  if (min_count >= 2)
    s << p << "  < " << min_count << " read pairs:\t" << group_stats.small
	   << "\t  (first pass: " << pass1_group_stats.small << ")\n";
  s << p << "  Stacked narrowly:\t" << group_stats.stacked << '\n';
  if (! anchors.empty())
    s << p << "  Neither anchored:\t" << group_stats.unanchored << '\n';
  s << p << "  Swamped by others:\t" << group_stats.swamped << '\n';

  s << blank
    << p <<  "Total groups emitted:\t" << group_stats.emitted << '\n';
}


bool rearrangement_grouper::
within(interval_multimap<feature>& features, const seqinterval& aln,
       int max_uncovered) {
  interval_multimap<feature>::iterator_pair range =
      features.intersecting_range(aln);

  int covered = 0;           // number of so-far covered bases within the read
  coord_t pos = aln.start(); // first as-yet uncovered position within the read

  for (interval_multimap<feature>::iterator it = range.first;
       it != range.second && pos <= aln.end(); ++it)
    if (it->first.end() >= pos) {
      coord_t newpos = it->first.limit();
      if (newpos > aln.limit())  newpos = aln.limit();

      if (it->first.start() > pos)
	covered += newpos - it->first.start();
      else
	covered += newpos - pos;

      pos = newpos;
    }

  return covered + max_uncovered >= aln.length();
}

void trim_clipped(const alignment& aln, int& start, int& end) {
  size_t cigar_len = aln.cigar_length();
  size_t i, j;

  start = end = 0;

  for (i = 0; i < cigar_len; i++) {
    cigar_op op = aln.cigar(i);
    if (op.opcode() == SOFT_CLIP)  start += op.length();
    else if (op.opcode() == HARD_CLIP) { }
    else  break;
  }

  for (j = cigar_len - 1; j >= i; j--) {
    cigar_op op = aln.cigar(j);
    if (op.opcode() == SOFT_CLIP)  end += op.length();
    else if (op.opcode() == HARD_CLIP) { }
    else  break;
  }
}

class base_count {
public:
  base_count() { clear(); }
  void clear() { for (int i = 0; i < 27; i++)  count[i] = 0; }
  size_t& operator[] (char base) { return count[index(base)]; }
  size_t operator[] (char base) const { return count[index(base)]; }
  size_t max_count() const { return (*std::max_element(count, &count[27])); }

private:
  static size_t index(char base) { return (base != '=')? (base - 'A') : 26; }
  size_t count[27];
};

static bool intersect(interval_multimap<feature>& map,
		      const string& name, const ::interval &region) {
    interval_multimap<feature>::iterator_pair range =
	map.intersecting_range(make_seqinterval(name, region));

    for (interval_multimap<feature>::iterator it = range.first;
	 it != range.second; ++it)
      if (it->first.end() >= region.pos5 && it->first.start() <= region.pos3)
	return true;

    return false;
}

inline bool nearly_proper(const alignment& aln) {
  return aln.rindex() == aln.mate_rindex() && natural_orientation(aln) &&
	 std::abs(aln.isize()) <= 500;
}

int rearrangement_grouper::
count_swamping_reads(rearr_group* group, const string& group_rname,
		     const ::interval& group_read_window) {
  typedef interval_multimap<rearr_group*>::iterator range_iterator;

  seqinterval group_seqival = make_seqinterval(group_rname, group_read_window);
  std::set<string> qnames;
  string buf;

  interval_multimap<rearr_group*>::iterator_pair
      range = active_by_readL.intersecting_range(group_seqival);

  for (range_iterator git = range.first; git != range.second; ++git) {
    rearr_group* g = git->second;
    if (intersect(g->readL, group_read_window) && g != group)
      for (rearr_group::iterator it = g->begin(); it != g->end(); ++it) {
	::interval rwindow(it->pos(), it->pos() + it->length());
	if (intersect(rwindow, group_read_window) && ! nearly_proper(*it))
	  qnames.insert(it->qname(buf));
      }
  }

  range = active_by_readH.intersecting_range(group_seqival);
  for (range_iterator git = range.first; git != range.second; ++git) {
    rearr_group* g = git->second;
    if (intersect(g->readH, group_read_window) && g != group)
      for (rearr_group::iterator it = g->begin(); it != g->end(); ++it) {
	::interval rwindow(it->mate_pos(), it->mate_pos() + it->length());
	if (intersect(rwindow, group_read_window) && ! nearly_proper(*it))
	  qnames.insert(it->qname(buf));
      }
  }

  return qnames.size();
}

void
rearrangement_grouper::annotate_intrachromosomal_deletions(rearr_group& group) {
  const alignment& aln = group.canonical;

  // Annotate intrachromosomal deletions that span repeats.
  if (aln.rindex() == aln.mate_rindex() && natural_orientation(aln)) {
    std::map<string, int> spanned;
    coord_t min_length =
	      group.overlapH.pos3 - group.overlapL.pos5 - group.max_insert;

    interval_multimap<feature>::iterator_pair range =
	transposons.intersecting_range(seqinterval(aln.rname(),
		group.overlapL.pos5-1, group.overlapL.pos5 + group.max_insert));

    for (interval_multimap<feature>::iterator it = range.first;
	 it != range.second; ++it)
      if (it->first.start() >= group.overlapL.pos5 &&
	  it->first.end()   <= group.overlapH.pos3 &&
	  it->first.length() >= min_length)
	spanned[it->second.name()]++;

    if (! spanned.empty()) {
      std::ostringstream note;
      for (std::map<string, int>::const_iterator it = spanned.begin();
	   it != spanned.end(); ++it) {
	if (it != spanned.begin())  note << ' ';
	note << it->first;
	if (it->second > 1)  note << '*' << it->second;
      }

      group.set_notes(note.str());
    }
  }
}

template <typename InputSamStream>
void rearrangement_grouper::group_alignments(InputSamStream& in) {
  alignment aln;
  seqinterval aln_ival, mate_ival;
  string seq;
  stopwatch elapsed;

  while (in >> aln) {
    read_stats.total++;

    // Apply symmetric read filters that discard read pairs equivalently
    // irrespective of whether they're considering a read or its mate.

    int flags = aln.flags();
    if (flags & PROPER_PAIRED) { read_stats.proper++; continue; }
    if (flags & (UNMAPPED|MATE_UNMAPPED)) { read_stats.unmapped++; continue; }
    if (aln.mapq() < min_quality) { read_stats.low_quality++; continue; }
    if (aln.aux("MQ", 255) < min_quality) { read_stats.low_mate_quality++; continue; }
    if (aln.rindex() == aln.mate_rindex() && std::abs(aln.zpos() - aln.mate_zpos()) < 10) { read_stats.near_mate++; continue; }

    aln_ival.assign(aln.rname(), aln.zpos(), aln.zpos() + aln.length());
    mate_ival.assign(aln.mate_rname(),
		     aln.mate_zpos(), aln.mate_zpos() + aln.length());

    if (discard_within_repeats &&
        (within_repeat(aln_ival) || within_repeat(mate_ival)))
      { read_stats.repeats++; continue; }

    if (discard_repeat_mapped && aln.aux("XT", '.') == 'R')
      { read_stats.repetitive++; continue; }

    if (! ignores.empty() &&
	(within(ignores, aln_ival) || within(ignores, mate_ival)))
      { read_stats.ignored++; continue; }

    const readgroup_info& info = readgroups.find(aln);

    // We may want to ignore apparent small intrachromosomal insertions,
    // which are likely to be artifacts.
    if (discard_apparent_insertions && apparent_insertion(aln, info))
      { read_stats.insertion++; continue; }

    int start, end;
    trim_clipped(aln, start, end);
    int mapped_length = aln.length() - start - end;

    if (start >= 2 && end >= 2 && mapped_length < min_clipped_length)
      { read_stats.clipped_short++; continue; }

    aln.seq(seq);
    base_count bases;
    size_t polytract_start = start, max_polytract = 0;
    for (size_t i = start; i < seq.length() - end; i++) {
      bases[seq[i]]++;
      if (i+1 >= seq.length() - end || seq[i+1] != seq[i]) {
	max_polytract = std::max(max_polytract, i - polytract_start + 1);
	polytract_start = i+1;
      }
    }

    if (max_polytract > max_polytract_frac * mapped_length)
      { read_stats.polytract++; continue; }
    if (bases.max_count() > max_samebase_frac * mapped_length)
      { read_stats.samebase++; continue; }

    if (less_than_mate(aln)) {
      // Process the lesser (by location) read in each pair.
      // Find the intervals subtended by the read pair, and add it
      // to any rearrangement groups that match those intervals.
      ::interval alnL(aln, aln.pos(), aln.strand(),
		    ref_length[aln.rindex()], info);
      ::interval alnH(aln, aln.mate_pos(), aln.mate_strand(),
		    ref_length[aln.mate_rindex()], info);

      int matched = 0;
      rearr_group_set::iterator_pair range = active.range_lower(aln);
      rearr_group_set::iterator group = range.first;
      while (group != range.second)
	if (group->overlapL.pos3 <= alnL.pos5 &&
	    group->canonical.strand() == aln.strand()) {
	  // Remove any clearly too small groups immediately.
	  if (group->lower_primary_count() == 0)
	    { group_stats.supponly++; group = active.erase(group); }
	  else if (group->lower_total_count() < min_count)
	    { group_stats.small++; group = active.erase(group); }
	  else
	    ++group;
	}
	else if (group->matches(aln, alnL, alnH)) {
	  group->insert(aln, alnL, alnH, info);
	  ++matched;
	  ++group;
	}
	else
	  ++group;

      if (matched == 0) {
	// If the pair is matched by no groups, create a new group.
	active.insert(rearr_group(aln, alnL, alnH, info, readgroups));
	group_stats.total++;
      }
    }
    else {
      // Process the greater (by location) read in each pair.
      // Record that the mate also passed filters and add mate information
      // to any groups that contain the read pair.

      ::interval alnL(aln, aln.mate_pos(), aln.mate_strand(),
		    ref_length[aln.mate_rindex()], info);
      ::interval alnH(aln, aln.pos(), aln.strand(),
		    ref_length[aln.rindex()], info);

      rearr_group::mate_iterator hint;
      rearr_group_set::iterator_pair range = active.range_higher(aln);
      for (rearr_group_set::iterator group = range.first;
	   group != range.second; ++group)
	if (group->matches_mate(aln, alnL, alnH, hint))
	  group->insert_mate(hint, aln);
    }
  }

  pass1_group_stats = group_stats;
  std::clog << "Input & read filtering:\t"
	    << std::setw(5) << elapsed << " seconds\n";
  elapsed.restart();

  // Apply all individual group filters to remaining groups, and recalculate
  // group intervals (now that mate information has been supplied).

  rearr_group_set::full_iterator it = active.begin();
  while (it != active.end()) {
    if (it->higher_primary_count() == 0) { group_stats.supponly++; goto erase; }
    if (it->higher_total_count() < min_count) {group_stats.small++; goto erase;}

    // Recalculate if the group changed during sychronisation.
    if (it->synchronise())  recalculate_group(*it);

    if (it->either_side_stacked()) { group_stats.stacked++; goto erase; }

    if (! anchors.empty()) {
      if (! intersect(anchors, it->rname(), it->overlapL) &&
	  ! intersect(anchors, it->mate_rname(), it->overlapH))
	{ group_stats.unanchored++; goto erase; }
    }

    ++it;
    continue;

  erase:
    it = active.erase(it);
  }

  std::clog << "Group filtering:\t" << std::setw(5) << elapsed << " seconds\n";
  elapsed.restart();

  for (it = active.begin(); it != active.end(); ++it) {
    active_by_readL.insert(std::make_pair(
	make_seqinterval(it->rname(), it->readL), &*it));
    active_by_readH.insert(std::make_pair(
	make_seqinterval(it->mate_rname(), it->readH), &*it));
  }

  // Annotate and print out groups, skipping those that are swamped on either
  // side by other groups based in intersecting read intervals.

  std::set<string> readnamesL, readnamesH;
  for (it = active.begin(); it != active.end(); ++it) {
    int max_swamping =
	std::max(count_swamping_reads(&*it, it->rname(), it->readL),
		 count_swamping_reads(&*it, it->mate_rname(), it->readH));
    if (max_swamping >= 2 * it->lower_total_count() &&
	! (intersect(retrotransposons, it->rname(), it->overlapL) ||
	   intersect(retrotransposons, it->mate_rname(), it->overlapH)))
      { group_stats.swamped++; continue; }

    // Update the group by adding any relevant notes and annotations.
    annotate_intrachromosomal_deletions(*it);

    out << *it << '\n';
    group_stats.emitted++;
  }

  std::clog << "Swamped filter: \t" << std::setw(5) << elapsed << " seconds\n";
  std::clog << '\n';
}

void rearrangement_grouper::recalculate_group(rearr_group& group) {
  bool first = true;
  for (rearr_group::iterator it = group.begin(); it != group.end(); ++it) {
    const readgroup_info& info = readgroups.find(*it);

    ::interval alnL(*it, it->pos(), it->strand(),
		    ref_length[it->rindex()], info);
    ::interval alnH(*it, it->mate_pos(), it->mate_strand(),
		    ref_length[it->mate_rindex()], info);

    if (first)  group.assign(*it, alnL, alnH, info, readgroups);
    else  group.reinsert(*it, alnL, alnH, info);
    first = false;
  }
}


int main(int argc, char** argv)
try {
  static const char version[] = "brass-group (Brass) " BRASS_VERSION;

  static const char copyright[] =
"Copyright (C) 2014-2017 Genome Research Ltd.\n"
"Portions copyright (C) 2018-2020 University of Glasgow.\n"
"This is free software; you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n"
"";

  static const char usage[] =
"Usage: brass-group [OPTION]... FILE [FILE] [FILE]\n"
"Options:\n"
"  -A FILE    Keep only groups with either side intersecting regions in FILE\n"
"  -d TYPE    Discard read pairs or groups matching condition TYPE\n"
"  -F FILE    Read annotation features from FILE (in BED or range format)\n"
"  -i RANGE   Omit groups in or near the locations encompassed by RANGE\n"
"  -I FILE      ...or locations encompassed by ranges listed in FILE\n"
"  -k TYPE    Keep read pairs or groups matching condition TYPE\n"
"  -m NUM     Use maximum insert size NUM unless specified by the library\n"
"  -n NUM     Omit groups containing fewer than NUM read pairs (default 2)\n"
"  -o FILE    Write rearrangement groups to FILE rather than standard output\n"
"  -q NUM     Discard read pairs with mapping quality less than NUM (default 1)\n"
"  -R FILE    Read retrotransposon features from FILE (in BED or range format)\n"
"  -s NAME    Use sample NAME for read pairs that are not in any read group\n"
"Conditions:\n"
"  insertion  Intrachromosomal insertions smaller than the insert (discarded)\n"
"  repeat     Groups touching listed repeat features (discarded)\n"
"  repetitive Read pairs marked as repetitively mapped (kept)\n"
"";

  if (argc >= 2) {
    string arg = argv[1];
    if (arg == "--help") { std::cout << usage; return EXIT_SUCCESS; }
    else if (arg == "--version")
      { std::cout << version << '\n' << copyright; return EXIT_SUCCESS; }
  }

  options opt;
  opt.output_filename = "-";  // standard output by default
  opt.max_insert = -1;  // an error if a read group needs to fall back to this
  opt.min_count = 2;
  opt.min_quality = 1;
  opt.discards["insertion"] = true;
  opt.discards["repeat"] = true;
  opt.discards["repetitive"] = false;

  int c;
  while ((c = getopt(argc, argv, ":A:d:F:i:I:k:m:n:o:q:R:s:")) >= 0)
    switch (c) {
    case 'A':  opt.anchor_filenames.push_back(optarg);  break;
    case 'F':  opt.feature_filenames.push_back(optarg);  break;
    case 'i':  opt.ignores.push_back(optarg);  break;
    case 'I':  opt.ignore_filenames.push_back(optarg);  break;
    case 'm':  opt.max_insert = atoi(optarg);  break;
    case 'n':  opt.min_count = atoi(optarg);  break;
    case 'o':  opt.output_filename = optarg;  break;
    case 'q':  opt.min_quality = atoi(optarg);  break;
    case 'R':  opt.retrotransposon_filenames.push_back(optarg);  break;
    case 's':  opt.default_sample = optarg;  break;

    case 'd':
    case 'k':
      if (opt.discards.find(optarg) == opt.discards.end())
	{ std::cerr << usage; return EXIT_FAILURE; }
      opt.discards[optarg] = (c == 'd');
      break;

    default:   std::cerr << usage; return EXIT_FAILURE;
    }

  int nfiles = argc - optind;
  if (nfiles == 0 || nfiles > 3) { std::cerr << usage; return EXIT_FAILURE; }

  isamstream in(argv[optind]);
  if (! in.is_open())
    throw sam::system_error("can't open ", argv[optind], errno);
  in.exceptions(std::ios::failbit | std::ios::badbit);

  std::ostringstream preamble;
  preamble << "# Rearrangement groupings, generated by " << version << "\n#\n"
	   << "# Input files:\n";

  for (int i = optind; i < argc; i++)
    preamble << "#INPUT\t" << argv[i] << '\n';

  collection headers;
  if (nfiles == 1) {
    in >> headers;
    rearrangement_grouper grouper(opt, headers);
    grouper.print_preamble(headers, preamble.str());
    grouper.group_alignments(in);
    grouper.print_trailer();
    grouper.log_statistics(std::clog);
  }
  else if (nfiles == 2) {
    isamstream in2(argv[optind+1]);
    if (! in2.is_open())
      throw sam::system_error("can't open ", argv[optind+1], errno);
    in2.exceptions(std::ios::failbit | std::ios::badbit);

    imergestream<> merged(in, in2);
    merged >> headers;
    rearrangement_grouper grouper(opt, headers);
    grouper.print_preamble(headers, preamble.str());
    grouper.group_alignments(merged);
    grouper.print_trailer();
    grouper.log_statistics(std::clog);
  }
  else {
    isamstream in2(argv[optind+1]);
    if (! in2.is_open())
      throw sam::system_error("can't open ", argv[optind+1], errno);
    in2.exceptions(std::ios::failbit | std::ios::badbit);

    isamstream in3(argv[optind+2]);
    if (! in3.is_open())
      throw sam::system_error("can't open ", argv[optind+2], errno);
    in3.exceptions(std::ios::failbit | std::ios::badbit);

    imergestream<> merged23(in2, in3);
    imergestream<imergestream<> > merged(in, merged23);
    merged >> headers;
    rearrangement_grouper grouper(opt, headers);
    grouper.print_preamble(headers, preamble.str());
    grouper.group_alignments(merged);
    grouper.print_trailer();
    grouper.log_statistics(std::clog);
  }

  return EXIT_SUCCESS;
}
catch (const std::exception& e) {
  std::cout << std::flush;
  std::cerr << "brass-group: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
catch (const char* what) {
  std::cout << std::flush;
  std::cerr << "brass-group: " << what << std::endl;
  return EXIT_FAILURE;
}
