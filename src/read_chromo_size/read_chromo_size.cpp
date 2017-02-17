/*    
 *    read_chromo_size.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;

string
basename2(string filename) {

	//cout << "bef: "<<filename<<endl;
	const string s(filename.substr(0, filename.find_last_of(".")));

	cout << s << endl;
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
	{
		//cout << "aft 1: "<<s.substr(final_slash + 1)<<endl;
		return s.substr(final_slash + 1);
	}
	else
	{
		//cout << "aft 2: "<<s<<endl;
		return s;
	}
}


void merge_sam(char* chromfile, char* chrom_size_file, char* chrom_names, map<string, pair<size_t, size_t > >& chrom_size)
{
	ifstream input_fs(chromfile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	ifstream testifs(chrom_size_file);

	ofstream output_fs;
	if (testifs.is_open())
	{
		testifs.close();
		output_fs.open(chrom_size_file,ios::app);
	}
	else
	{
		testifs.close();
		output_fs.open(chrom_size_file);
	}

	ifstream testifs2(chrom_names);
	ofstream output_fs2;
	if (testifs2.is_open())
	{
		testifs2.close();
		output_fs2.open(chrom_names,ios::app);
	}
	else
	{
		testifs2.close();
		output_fs2.open(chrom_names);
	}

	if( !output_fs ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	} 

	if( !output_fs2 ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	}

	string line;
	size_t count=0;
	size_t countline = 0;
	getline(input_fs, line);
	while(getline(input_fs, line))
	{
		if (line == "")
			continue;
		count += line.length();
		++countline;
	}

	string chromname = basename2(chromfile);
	
	if (count > 4294967290)
	{
		cerr<<"chromosome length too long: " << chromfile <<'\t' << count <<endl;

		exit(2);
	}

	chrom_size[chromname] = make_pair(count, countline);


	//output_fs<<chromname<<"\t"<<count<<endl;

	//output_fs2 << chromname << endl;

	input_fs.close();
	output_fs.close();
	output_fs2.close();
}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	} 

	map<string, pair<size_t, size_t > > chrom_size;

	for(int i=1;i<argc-4;i++)
	{
		merge_sam(argv[i],argv[argc-2], argv[argc-1], chrom_size);
	}

	map<string, pair<size_t, size_t > >::iterator chrom_size_iter;

	ofstream output_fai(argv[argc-4]);

	ofstream output_fs_head(argv[argc-3]);

	ofstream output_fs(argv[argc-2]);

	ofstream output_fs2(argv[argc-1]);

	size_t offset = 0;

	for (chrom_size_iter = chrom_size.begin(); chrom_size_iter != chrom_size.end(); ++chrom_size_iter)
	{
		size_t linelen = chrom_size_iter->second.first / chrom_size_iter->second.second;

		if (chrom_size_iter->second.first % chrom_size_iter->second.second)
			++linelen;

		output_fs<<chrom_size_iter->first <<'\t' << chrom_size_iter->second.first << endl;

		output_fs2<<chrom_size_iter->first << endl;

		output_fs_head << "@SQ\tSN:"<<chrom_size_iter->first<<"\tLN:"<<chrom_size_iter->second.first <<endl;

		output_fai << chrom_size_iter->first<<"\t"<<chrom_size_iter->second.first << '\t' << offset + chrom_size_iter->first.length() + 2<<'\t' << linelen << '\t' << linelen + 1 << endl;

		offset = offset + chrom_size_iter->first.length() + 2 + chrom_size_iter->second.second + chrom_size_iter->second.first;
	}

	//sort_sam(argv[argc-1]);
	return 0;
}
