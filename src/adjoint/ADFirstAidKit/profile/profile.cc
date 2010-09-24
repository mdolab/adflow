#include "reader.h"
#include <fstream>
#include <string>
#include <stdexcept>
#include <map>
#include <stack>
#include <set>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <memory>
#include <iterator>
#include <limits>

using namespace std;

typedef unsigned long long int clicks;

int main(int argc, char **argv)
{
	bool remove_subtime = true;
	beginning:
	if (argc < 3) {
		cout << "Usage: profile tapenade.prof ouput.csv\n";
		exit(-1);
	}
	if (string("-nosubtime") == argv[1]) {
		remove_subtime = false;
		argv[1] = argv[2];
		argv[2] = argv[3];
		main(argc, argv);
		return 0;
	}
	const string filename(argv[1]);
	ReadingOnlyInAdjointBranch reader(filename);
	
	stack<TapenadeLogEvent> call_stack;
	int first_free_function_id = 0;
	map<string,int> function_ids;
	typedef map<string, int> functions_ids_t;
	typedef functions_ids_t::iterator ids_iter;
	vector<int> taping;
	vector<int> snapshoting;
	vector<clicks> snapshottime;
	vector<clicks> fwdtime;
	vector<clicks> totaltime;
	vector<int> totalsize;
	vector<int> call_count;
	vector<int> snap_count;
	stack<clicks> subtime;

	taping.reserve(2000);
	fwdtime.reserve(2000);
	snap_count.reserve(2000);
	snapshoting.reserve(2000);
	snapshottime.reserve(2000);
	call_count.reserve(2000);
	totalsize.reserve(2000);
	totaltime.reserve(2000);
	int event_count = 0;
	while (! reader.eof()) {
		TapenadeLogEvent event(reader.readEvent());
		int kind = event.get_kind();
		switch(kind) {
			case BEGIN:
				subtime.push(0);
				call_stack.push(event);
				break;
			case BEGINSNAP:
				call_stack.push(event);
				break;
		}
		string &name = event.get_function_name();
		clicks time = event.get_time();
		unsigned int size = event.get_size();
		// Allocate a function id if necessary
		unsigned int id;
 		if (function_ids.find(name) == function_ids.end()) {
			id = first_free_function_id++;
			function_ids[name] = id;
			
			taping[id] = 0;
			fwdtime[id] = 0;
			snap_count[id] = 0;
			snapshoting[id] = 0;
			snapshottime[id] = 0;
			call_count[id] = 0;
			totalsize[id] = 0;
			totaltime[id] = 0;
		}
		else
			id = function_ids[name];

		// ALLOCATE
		// JOINT MODE taping & time
		if (kind == ENDFWD) {
			taping[id] += size-call_stack.top().get_size();
			fwdtime[id] += time-call_stack.top().get_time();
		}
		// SNAPSHOT taping & time
		if (kind == ENDSNAP) {
			snap_count[id]++;
			snapshoting[id] += size-call_stack.top().get_size();
			snapshottime[id] += size-call_stack.top().get_size();
		}
		// ORIG, FWD, BWD taping & time
		if (kind == END) {
			call_count[id]++;
			totalsize[id] += size-call_stack.top().get_size();
			clicks t = time-call_stack.top().get_time();
			totaltime[id] += t + (remove_subtime ? - subtime.top() : 0 );
			subtime.pop();
			if (subtime.size())
				subtime.top() += t;
		}
		switch(event.kind) {
			case ENDSNAP:
				call_stack.pop();
				break;
			case END:
				call_stack.pop();
				break;
		}
		event_count ++;
	}
	ofstream output_stream(argv[2]);
	ids_iter cur = function_ids.begin(), end = function_ids.end();
	output_stream << setw(25) << "name" ; // 1
	output_stream << setw(14) << "call_count" ; // 2
	output_stream << setw(14) << "snap_count" ; // 3
	output_stream << setw(14) << "taping" ; // 4
	output_stream << setw(14) << "snapshoting" ; // 5

	output_stream << setw(14) << "fwdtime" ; // 6
	output_stream << setw(14) << "totaltime" ; // 7
	output_stream << setw(14) << "totalsize" ; // 8
	output_stream << setw(14) << "snapshottime" ; // 8
	output_stream << endl;
	for (; cur != end; ++cur) {
		const string &name = (*cur).first;
		int id = (*cur).second;
		output_stream << setw(25) << name ; // 1
		output_stream << setw(14) << call_count[id]; // 2
		output_stream << setw(14) << snap_count[id]; // 3
		output_stream << setw(14) << taping[id]; // 4
		output_stream << setw(14) << snapshoting[id]; // 5

		output_stream << setw(14) << fwdtime[id]; // 6
		output_stream << setw(14) << totaltime[id]; // 7
		output_stream << setw(14) << totalsize[id]; // 8
		output_stream << setw(14) << snapshottime[id]; // 8
		output_stream << endl;
	}
}
		
