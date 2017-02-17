#include "SamRec.h"
#include "JunctionSeed.h"
#include <algorithm>
#include <vector>

//
SamRec::SamRec(const string& line, int min_ins) : matched_id(-1),  wrong_format(false), tagidx(-1), filter_score(0), best(0), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(0), ave_junc_mis(0), canon_count(0), 
	noncanon_count(0), canon_rate(1), issmallins(false), issmalldel(false), iscanonical(false), 
	issemicanonical(false), isnoncanoical(false), isunmapped(false), isclipped(false), /*cur_line(line),*/
	/*is_insert(false),*/ mappedlen(0), intron_size(0), min_anchor(-1), max_anchor(0), pair_rate(0), m_contig_len(0), 
	fusion_prefix_len(0), fusion_suffix_len(0), is_swapped(false), mate_offset(0), mate_offset2(0), mate_diff(0), 
	mate_diff2(0), mate_match("*"), mate_match2("*"), paired_type(SINGLE), start(0), start2(0), strand_t(0), strand_t2(0),
	forward_count(0), reverse_count(0), fusion_prefix_paired(false), fusion_suffix_paired(false), need_swap(false), insertlen(0), insertlen1(0), insertlen2(0)
{
	if (line.find("seq.10001997/") != string::npos)
	{
		cout << line << endl;
	}

	char tagname[1000], chrom[100], mapped[100], seq[1000], qual_chr[1000], alters_chr[1000];

	char chrom2[100], mapped2[100];

	char mate_match_chr[1000];

	string alterstr = "";

	char strand[100];

	sscanf(line.c_str(), "%s\t%s\t%s", tagname, chrom, strand);

	if (string(strand) == "+" || string(strand) == "-")
		is_fusion = true;
	else
		is_fusion = false;

	if (line.find("ZF:Z:FUS") != string::npos)
		is_fusion_newfmt = true;
	else
		is_fusion_newfmt = false;

	if (is_fusion)
	{
		//14~R_7/2	chr1	+	11184637	54M	chr11	+	61646784	21M	GGCGCAGATCTTCATGGCCTGTTAGAAGGAAAACACACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGGATGAA	
		//283348985863448157653149274159715381432428925324118868345815774332128375427~CTGCATCACACGCTCATCCTGGCGCAGATCTTCATGGCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTG
		//CCTGAGAGATGGCCAGGATGAAGGCGGCCAGGGCACTGGGCACCCAGCCAGGACCCAGGAGGTAGATAAGGAGCC	NM:i:2

		size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tagname, chrom, &strand1, &start, mapped, 
			chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match);

		//sscanf(tagname, "%llu", &tagidx);

		mapped_seq = seq;

		qual_str = qual_chr; 

		syn_seq = qual_chr;

		size_t last_idx = syn_seq.find_last_of("~");

		qual_str = qual_str.substr(0, last_idx);

		syn_seq = syn_seq.substr(last_idx + 1, syn_seq.length() - last_idx - 1);

		tag_name = tagname;

		if (strand1 == '+')
			strand_t = 0;
		else
			strand_t = 16;

		if (strand2 == '+')
			strand_t2 = 0;
		else
			strand_t2 = 16;

		if (string(chrom) < string(chrom2) || (string(chrom) == string(chrom2) && start <= start2))
		{
			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			chrom_name = chrom;
			splice_way = mapped;
			ori_splice_way = mapped;

			splice_way2 = mapped2;
			ori_splice_way2 = mapped2;
			chrom_name2 = chrom2;
		}
		else
		{
			if (strand1 == '+')
				strand1 = '-';
			else
				strand1 = '+';

			if (strand2 == '+')
				strand2 = '-';
			else
				strand2 = '+';

			swap(strand1, strand2);

			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;

			
			swap(start, start2);
			
			chrom_name = chrom2;//chrom;
			splice_way = mapped2;//mapped;
			ori_splice_way = mapped2;//mapped;

			splice_way2 = mapped;//mapped2;
			ori_splice_way2 = mapped;//mapped2;
			chrom_name2 = chrom;//chrom2;

			is_swapped = true;
		}

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}	

		size_t index = 0;

		string flag_str = " ";

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(start, -maplen));

					insertlen1 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec.push_back(make_pair(start + maplen, 0));
				}
			}
			else if (flag_str[0] == 'M')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));

					insertlen1 += -maplen;
				}
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
				{
					spliceway_vec.back().second = -maplen;

					insertlen1 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec.back().second += -maplen;

					insertlen1 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		mappedlen1 = mappedlen;

		fusion_prefix_len = mappedlen1;

		//segment 2	

		if (mapped2[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way2.find("I") != string::npos))
			issmallins = true;
		if ((splice_way2.find("D") != string::npos))
			issmalldel = true;
		if (splice_way2.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way2.length(); ++i)
			if (splice_way2[i] == 'D')
				splice_way2[i] = 'N';

		if ((splice_way2.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way2[splice_way2.length() - 1] == 'S')
			{
				size_t m_idx = splice_way2.find_last_of("M");

				splice_way = splice_way2.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way2.find("S");

			if (s_index != string::npos)
			{
				splice_way2.substr(s_index + 1, splice_way2.length() - s_index - 1);
			}
		}

		size_t index2 = 0;

		string flag_str2 = " ";

		while (true)
		{
			if (index2 >= splice_way2.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way2.c_str() + index2, "%d%c", &maplen, &flag);

			if (flag_str2[0] == ' ')
			{
				if (flag == 'I')
				{
					spliceway_vec2.push_back(make_pair(start2, -maplen));

					insertlen2 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(start2, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec2.push_back(make_pair(start2 + maplen, 0));
				}
			}
			else if (flag_str2[0] == 'M')
			{
				if (flag == 'I')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second, -maplen));

					insertlen2 += -maplen;
				}
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec2.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'N')
			{
				if (flag == 'I')
				{
					spliceway_vec2.back().second = -maplen;

					insertlen2 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec2.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec2.back().second += -maplen;

					insertlen2 += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str2[0] = flag;

			index2 = splice_way2.find(flag_str2, index2) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end2 = spliceway_vec2.back().first + spliceway_vec2.back().second - 1;

		mappedlen2 = mappedlen - mappedlen1;

		fusion_suffix_len = mappedlen2;
		//fusion_prefix_st, fusion_prefix_end, fusion_suffix_st, fusion_suffix_end;

		if (strand1 == '+')
		{
			fusion_prefix_st = start;

			fusion_prefix_end = end;
		}
		else
		{
			fusion_prefix_st = end;

			fusion_prefix_end = start;
		}

		if (strand2 == '+')
		{
			fusion_suffix_st = start2;
			
			fusion_suffix_end = end2;
		}
		else
		{
			fusion_suffix_st = end2;
			
			fusion_suffix_end = start2;
		}
	}
	else
	{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       
		//TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		//seq.1000000/2   0       chr22   51213810        255     75M     *       0       0       TCTACTCAGTCCCAGCATGTGCCCCCTAGGGAAAGGCAGCAGGCAACAGTGGAGGGAGGCCTTGGGAGCTTGTCT     
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    NM:i:2  NM:i:2  MD:Z:16G3C55

		//seq.46/1	0	chr9	123606192	255	75M	*	0	0	GACCCAGTTAACCTTTTAGGTGGTAGGCCAGATGACAGAGGGCGGCATCTGACCCGGAAGGCCCTCTGCAACCAG	
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	MD:Z:A23GA5G45G1

		//JOHNLENNON_0006:1:1:18984:1271#0/2	435	chr7	52225290	255	50S19M	=	18546378	33678863	
		//CTCGGAGCTGATGCTTGGCGGAACCAACACACTAGTGCTGNACAACACGTGTCAGGACACGCTGCTGGC	eeeeLfffffhhhehWffffhhfghgghffhhhRh\`b^^Ba`^`bbbd_fbJdd`^`WZNX[[\^``_	NM:i:4	ZF:Z:FUS_18546378_52225308(--)

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand_t, chrom, &start, &confid, 
			mapped, mate_match_chr, &mate_offset, &mate_diff, seq, qual_chr, &mis_match, alters_chr);

		strand_t = strand_t & IS_REVERSE;

		if (mis_match)
		{
			size_t mis_idx = line.find("MD:Z:");

			if (mis_idx != string::npos)
			{
				mis_info_str = line.substr(mis_idx + 5);

				//bool was_mis = false;

				size_t current_offset = 0;

				size_t current_num = 0;

				for (size_t i = 0; i < mis_info_str.length(); ++i)
				{
					if (mis_info_str[i] >= '0' && mis_info_str[i] <= '9')
					{
						current_num = current_num * 10 + mis_info_str[i] - '0';

						if (i + 1 < mis_info_str.length() && (mis_info_str[i + 1] < '0' || mis_info_str[i + 1] > '9'))
							current_offset += current_num;
					}
					else
					{
						current_num = 0;

						mis_info.push_back(make_pair(current_offset++, mis_info_str[i]));
					}
				}
				//mis_info;
			}
		}

		if (line.find("ZF:Z:FUS") != string::npos)
		{
			is_fusion_newfmt = true;

			fusion_newfmt_str = line.substr(line.find("ZF:Z:FUS"));
		}
		else
			is_fusion_newfmt = false;

		//sscanf(tagname, "%llu", &tagidx);

		tag_name = tagname;
		chrom_name = chrom;
		splice_way = mapped;
		ori_splice_way = mapped;
		mapped_seq = seq;
		qual_str = qual_chr; 
		mate_match = mate_match_chr;

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		size_t index = 0;

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				suffixS = splice_way.substr(m_idx + 1);

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				prefixS = splice_way.substr(0, s_index + 1);

				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		if (isexonic)
		{
			filter_score = 5;

			junc_anchor_len = 500;

			ave_junc_mis = 0;

			pair_rate = 5;

			junc_hits = 0;
		}

		string flag_str = " ";

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(start, -maplen));

					insertlen += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec.push_back(make_pair(start + maplen, 0));
				}
			}
			else if (flag_str[0] == 'M')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));

					insertlen += -maplen;
				}
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
				{
					spliceway_vec.back().second = -maplen;

					insertlen += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec.back().second += -maplen;

					insertlen += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		if (mis_info.size() && mis_info.back().first > mappedlen)
		{
			cout << "mismatch position large than mapped length?"<<endl;
		}
	}
}

SamRec::SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
	const string& matematch, size_t mateoffest, int matediff, string line) : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
	mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offset(mateoffest), mate_diff(matediff), wrong_format(false), tagidx(tidx), filter_score(0), best(0), alters(alt), qual_str(qualstr), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(100000000), ave_junc_mis(100000000), canon_count(0), noncanon_count(0), canon_rate(1), issmallins(false), issmalldel(false), /*cur_line(line),*/
	/*is_insert(false),*/ mappedlen(0), mappedlen1(0), mappedlen2(0), intron_size(0), need_swap(false), insertlen(0), insertlen1(0), insertlen2(0)
{
	//sscanf(tag_name.c_str(), "%llu", &tagidx);

	size_t index = 0;

	if ((splice_way.find("I") != string::npos))
		issmallins = true;
	else if ((splice_way.find("D") != string::npos))
		issmalldel = true;
	else if (splice_way.find("N") != string::npos)
		isspliced = true;
	else
		isexonic = true;		

	string flag_str = " ";

	while (true)
	{
		if (index >= splice_way.length())
			break;

		int maplen;

		char flag;

		sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

		if (flag_str[0] == ' ')
		{
			if (flag == 'I')
				spliceway_vec.push_back(make_pair(start, -maplen));
			else if (flag == 'M')
			{
				spliceway_vec.push_back(make_pair(start, maplen));

				mappedlen += maplen;
			}
			else if (flag == 'N')
			{
				cout<<"start with N?"<<endl;
				spliceway_vec.push_back(make_pair(start + maplen, 0));
			}
		}
		else if (flag_str[0] == 'M')
		{
			if (flag == 'I')
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
			else if (flag == 'M')
			{
				cout << "continue Ms?"<<endl;
				spliceway_vec.back().second += maplen;

				mappedlen += maplen;
			}
			else if (flag == 'N')
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
		}
		else if (flag_str[0] == 'N')
		{
			if (flag == 'I')
				spliceway_vec.back().second = -maplen;
			else if (flag == 'M')
			{
				spliceway_vec.back().second = maplen;

				mappedlen += maplen;
			}
			else if (flag == 'N')
			{
				cout << "continue Ns?"<<endl;
				spliceway_vec.back().first += maplen;
			}
		}
		else if (flag_str[0] == 'I')
		{
			if (flag == 'I')
			{
				cout << "continue Is?"<<endl;
				spliceway_vec.back().second += -maplen;
			}
			else if (flag == 'M')
			{
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

				mappedlen += maplen;
			}
			else if (flag == 'N')
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
		}

		flag_str[0] = flag;

		index = splice_way.find(flag_str, index) + 1;
	}

	end = spliceway_vec.back().first + spliceway_vec.back().second - 1;
}


void SamRec::Set(const string& line, int min_ins) 
{
	//cout << line << endl;

	Clear();

	//cur_line = line;

	//if (line.find("seq.10001997/") != string::npos)
	//{
	//	cout << line << endl;
	//}

	char tagname[1000], chrom[100], mapped[100], seq[1000], qual_chr[1000], alters_chr[1000];

	char chrom2[100], mapped2[100];

	char mate_match_chr[1000];

	string alterstr = "";

	char strand[100];

	sscanf(line.c_str(), "%s\t%s\t%s", tagname, chrom, strand);

	if (string(strand) == "+" || string(strand) == "-")
		is_fusion = true;
	else
		is_fusion = false;

	if (line.find("ZF:Z:FUS") != string::npos)
		is_fusion_newfmt = true;
	else
		is_fusion_newfmt = false;

	if (is_fusion)
	{
		//14~R_7/2	chr1	+	11184637	54M	chr11	+	61646784	21M	GGCGCAGATCTTCATGGCCTGTTAGAAGGAAAACACACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGGATGAA	
		//283348985863448157653149274159715381432428925324118868345815774332128375427~CTGCATCACACGCTCATCCTGGCGCAGATCTTCATGGCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTG
		//CCTGAGAGATGGCCAGGATGAAGGCGGCCAGGGCACTGGGCACCCAGCCAGGACCCAGGAGGTAGATAAGGAGCC	NM:i:2
		
		char genestrand[100], flankseq[100];
		//\t%s\t%s , genestrand, flankseq

		//HWI-EAS217:4:1:126:890#0/1	chr5	+	13638058	1M3I73M	chr13	+	37500670	23M	CGTTTGGGAGCCTGAGCTGGAACTGAAGCTGGGGCTGCAGCCTGGGCCTTGGTTTGATCCTTGGCCTTGGCCTTGGCCTTGGCCTTGGCCTTTGGCCGGC	
		//`a^aaa``a_a`aa^a`aa]``a````\a^`\\]a````]`a][R\^^\^^VXS]TT[^]SSPP[]ZRPE\^TSRN]^URQQ`XSMKM]YWOOJM`WKPY	NM:i:5	XF:Z:CTGC,

		//size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu\t%s", tagname, chrom, &strand1, &start, mapped, 
		//	chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match, flankseq);

		////cout << flankseq<<endl;

		//if (read_count == 13)
		//{
		//	string xs_tag_str = genestrand;

		//	xf_tag = flankseq;

		//	if (xf_tag.find("XF:Z") != string::npos/* && xs_tag_str.find("XS:A") != string::npos*/)
		//	{
		//		if (xf_tag[xf_tag.length() - 1] != ',')
		//		{
		//			//cout << line << endl;
		//			xf_tag.append(",");
		//		}
		//		//sscanf(genestrand, "XS:A:%c", &xs_tag);

		//		xs_tag = '+';

		//		size_t count = (xf_tag.length() - 5) / 5;

		//		if ((xf_tag.length() - 5) % 5)
		//		{
		//			cout << line << endl;
		//		}

		//		for (size_t i = 0; i < count; ++i)
		//		{

		//			flankstrings.push_back(xf_tag.substr(i * 5 + 5, 4));
		//			//cout << flankstrings.back()<<endl;
		//		}
		//	}
		//	else
		//		xf_tag.clear();
		//}

		size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tagname, chrom, &strand1, &start, mapped, 
			chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match);

		//cout << flankseq<<endl;

		size_t xfidx = line.find("XF:Z");

		if (xfidx != string::npos)
		{
			string xs_tag_str = genestrand;

			sscanf(line.c_str() + xfidx, "%s", flankseq);

			xf_tag = flankseq;

			if (xf_tag[xf_tag.length() - 1] != ',')
			{
				//cout << line << endl;
				xf_tag.append(",");
			}
			//sscanf(genestrand, "XS:A:%c", &xs_tag);

			xs_tag = '+';

			size_t count = (xf_tag.length() - 5) / 5;

			if ((xf_tag.length() - 5) % 5)
			{
				cout << line << endl;
			}

			for (size_t i = 0; i < count; ++i)
			{

				flankstrings.push_back(xf_tag.substr(i * 5 + 5, 4));
				//cout << flankstrings.back()<<endl;
			}
		}

		//sscanf(tagname, "%llu", &tagidx);

		mapped_seq = seq;

		//cout << mapped_seq.length() << endl;

		qual_str = qual_chr; 

		syn_seq = qual_chr;

		size_t last_idx = syn_seq.find_last_of("~");

		qual_str = qual_str.substr(0, last_idx);

		syn_seq = syn_seq.substr(last_idx + 1, syn_seq.length() - last_idx - 1);

		tag_name = tagname;

		tag_base_name = tag_name.substr(0, tag_name.length() - 1);

		if (tag_name[tag_name.length() - 1] == '1' || tag_name[tag_name.length() - 1] == 'a')
			end_id = 1;
		else if (tag_name[tag_name.length() - 1] == '2' || tag_name[tag_name.length() - 1] == 'b')
			end_id = 2;

		if (strand1 == '+')
			strand_t = 0;
		else
			strand_t = 16;

		if (strand2 == '+')
			strand_t2 = 0;
		else
			strand_t2 = 16;

		if (string(chrom) < string(chrom2) || (string(chrom) == string(chrom2) && start <= start2))
		{
			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			chrom_name = chrom;
			splice_way = mapped;
			ori_splice_way = mapped;

			splice_way2 = mapped2;
			ori_splice_way2 = mapped2;
			chrom_name2 = chrom2;
		}
		else
		{
			//cout << "swap" << endl;
			if (strand1 == '+')
				strand1 = '-';
			else
				strand1 = '+';

			if (strand2 == '+')
				strand2 = '-';
			else
				strand2 = '+';

			swap(strand1, strand2);

			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			//cout << "swap1a" << endl; 
			swap(start, start2);
			
			chrom_name = chrom2;//chrom;
			splice_way = mapped2;//mapped;
			ori_splice_way = mapped2;//mapped;

			splice_way2 = mapped;//mapped2;
			ori_splice_way2 = mapped;//mapped2;
			chrom_name2 = chrom;//chrom2;
			//cout << "swap1b" << endl; 
			is_swapped = true;

			//cout << "mappedlen1:"<<mappedlen1<<endl;
			//cout << "mappedlen2:"<<mappedlen2<<endl;
			//string mapseq1 = mapped_seq.substr(0, mappedlen1);

			//string mapseq2 = mapped_seq.substr(mapped_seq.length() - mappedlen2, mappedlen2);

			//mapped_seq = mapseq2 + mapseq1;
			//cout << "swap1c" << endl; 
			//string qualstr1 = qual_str.substr(0, mappedlen1);

			//string qualstr2 = qual_str.substr(qual_str.length() - mappedlen2, mappedlen2);

			//qual_str = qualstr2 + qualstr1;

			//cout << "swap1" << endl; 
			reverse(flankstrings.begin(), flankstrings.end());

			//cout << "swap2" << endl;
			if (xs_tag == '+')
				xs_tag = '-';
			else if (xs_tag == '-')
				xs_tag = '+';
			//cout << "swap3" << endl;
			for (size_t i = 0; i < flankstrings.size();++i)
			{
				//cout << flankstrings[i]<<endl;

				flankstrings[i] = revcomp(flankstrings[i]);
			}

			//cout << "swap finished" << endl;

		}

		strand_t |= IS_PRIMARY;

		strand_t2 |= IS_PRIMARY;

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}

			//ori_splice_way = splice_way;
		}	

		int insertlen1 = 0, insertlen2 = 0;

		size_t index = 0;

		string flag_str = " ";

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(start, -maplen));
					insertlen1 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec.push_back(make_pair(start + maplen, 0));
				}
			}
			else if (flag_str[0] == 'M')
			{
				if (flag == 'I')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
					insertlen1 +=  -maplen;
				}
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
				{
					spliceway_vec.back().second = -maplen;
					insertlen1 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec.back().second += -maplen;
					insertlen1 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		mappedlen1 = mappedlen;

		fusion_prefix_len = mappedlen1;

		//segment 2	

		if (mapped2[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way2.find("I") != string::npos))
			issmallins = true;
		if ((splice_way2.find("D") != string::npos))
			issmalldel = true;
		if (splice_way2.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way2.length(); ++i)
			if (splice_way2[i] == 'D')
				splice_way2[i] = 'N';

		if ((splice_way2.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way2[splice_way2.length() - 1] == 'S')
			{
				size_t m_idx = splice_way2.find_last_of("M");

				splice_way2 = splice_way2.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way2.find("S");

			if (s_index != string::npos)
			{
				splice_way2 = splice_way2.substr(s_index + 1, splice_way2.length() - s_index - 1);
			}

			//ori_splice_way2 = splice_way2;
		}

		size_t index2 = 0;

		string flag_str2 = " ";

		while (true)
		{
			if (index2 >= splice_way2.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way2.c_str() + index2, "%d%c", &maplen, &flag);

			if (flag_str2[0] == ' ')
			{
				if (flag == 'I')
				{
					spliceway_vec2.push_back(make_pair(start2, -maplen));

					insertlen2 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(start2, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec2.push_back(make_pair(start2 + maplen, 0));
				}
			}
			else if (flag_str2[0] == 'M')
			{
				if (flag == 'I')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second, -maplen));

					insertlen2 +=  -maplen;
				}
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec2.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'N')
			{
				if (flag == 'I')
				{
					spliceway_vec2.back().second = -maplen;

					insertlen2 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec2.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec2.back().second += -maplen;
					insertlen2 +=  -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str2[0] = flag;

			index2 = splice_way2.find(flag_str2, index2) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end2 = spliceway_vec2.back().first + spliceway_vec2.back().second - 1;

		mappedlen2 = mappedlen - mappedlen1;

		fusion_suffix_len = mappedlen2;
		//fusion_prefix_st, fusion_prefix_end, fusion_suffix_st, fusion_suffix_end;

		if (strand1 == '+')
		{
			fusion_prefix_st = start;

			fusion_prefix_end = end;
		}
		else
		{
			fusion_prefix_st = end;

			fusion_prefix_end = start;
		}

		if (strand2 == '+')
		{
			fusion_suffix_st = start2;
			
			fusion_suffix_end = end2;
		}
		else
		{
			fusion_suffix_st = end2;
			
			fusion_suffix_end = start2;
		}

		
		mappedlen1 = mappedlen1 - insertlen1;

		mappedlen2 = mappedlen2 - insertlen2;

		//cout << mappedlen1 << endl;
		//cout << insertlen1 << endl;
		//cout << mappedlen2 << endl;
		//cout << insertlen2 << endl;

		//cout << mapped_seq.length() << endl;

		//cout << qual_str.length() << endl;
		if (is_swapped)
		{
			//cout << "mappedlen1:"<<mappedlen1<<endl;
			//cout << "mappedlen2:"<<mappedlen2<<endl;
			string mapseq1 = mapped_seq.substr(0, mappedlen2);

			string mapseq2 = mapped_seq.substr(mapped_seq.length() - mappedlen1, mappedlen1);

			mapped_seq = mapseq2 + mapseq1;
			
			//cout << "swap1c" << endl; 
			string qualstr1 = qual_str.substr(0, mappedlen2);

			string qualstr2 = qual_str.substr(qual_str.length() - mappedlen1, mappedlen1);

			qual_str = qualstr2 + qualstr1;
		}

		fusion_prefix_len = mappedlen1;

		fusion_suffix_len = mappedlen2;

		//cout << mapped_seq.length() << endl;

		//cout << qual_str.length() << endl;
	}
	else
	{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       
		//TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		//seq.1000000/2   0       chr22   51213810        255     75M     *       0       0       TCTACTCAGTCCCAGCATGTGCCCCCTAGGGAAAGGCAGCAGGCAACAGTGGAGGGAGGCCTTGGGAGCTTGTCT     
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    NM:i:2  NM:i:2  MD:Z:16G3C55

		//seq.46/1	0	chr9	123606192	255	75M	*	0	0	GACCCAGTTAACCTTTTAGGTGGTAGGCCAGATGACAGAGGGCGGCATCTGACCCGGAAGGCCCTCTGCAACCAG	
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	MD:Z:A23GA5G45G1

		//UNC13-SN749_0076:6:1101:2089:2682#0/1	0	chr14	93424691	138	21M3988N41M355N38M	*	0	0	
		//CTGGTTGAACACGATAGCCATCTCGTGAGAGTTGGTGCCATGAGCCACTCTGGTTTTGCAAATGAATGGGAAAGTCAAGCCGTTCTTCTCCAGCAGCCGC	
		//___c\a^cecceeg\d_fddghfg]be_dXcbeffacfgf^egfghfeghhfhadfghcffffd_dfgegedddHVZ]bbba]a`abbbb``bbbaaaa[	NM:i:0	-	CTAC,CTAC,

		char genestrand[100], flankseq[100];
		//\t%s\t%s , genestrand, flankseq

		//HWI-EAS217:4:1:36:1006#0/1	329	chr1	160236111	90	100M	*	0	0	CGCCTATGCAAGAACAGCTTAAGACCAGTCAGTGGTTGCTCCTACCCATTCAGTGGCCTGAGCAGTGGGAGCTGCAGACCAGTCTTCCGTGGCAGGCTGC	
		//aaaabbbbaababbaa_abba]aaaabbabbbab`_aba`baaZ_a`Paaba`^_^`aaa_aa^_P^USYP_\\^T`^Sa]XU`[TY\BBBBBBBBBBBB	NM:i:5	IH:i:21	HI:i:3

		//size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s\t%s", tagname, &strand_t, chrom, &start, &confid, 
		//	mapped, mate_match_chr, &mate_offset, &mate_diff, seq, qual_chr, &mis_match, genestrand, flankseq);

		//if (read_count == 14)
		//{
		//	string xs_tag_str = genestrand;

		//	xf_tag = flankseq;

		//	if (xf_tag.find("XF:Z") != string::npos && xs_tag_str.find("XS:A") != string::npos)
		//	{
		//		sscanf(genestrand, "XS:A:%c", &xs_tag);

		//		size_t count = (xf_tag.length() - 5) / 5;

		//		for (size_t i = 0; i < count; ++i)
		//		{
		//			flankstrings.push_back(xf_tag.substr(i * 5 + 5, 4));
		//		}
		//	}
		//	else
		//	{
		//		xs_tag_str.clear();

		//		xf_tag.clear();
		//	}
		//}

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu", tagname, &strand_t, chrom, &start, &confid, 
			mapped, mate_match_chr, &mate_offset, &mate_diff, seq, qual_chr, &mis_match);

		size_t xfidx = line.find("XF:Z");

		size_t xsidx = line.find("XS:A");

		if (xfidx != string::npos)
		{
			sscanf(line.c_str() + xfidx, "%s", flankseq);

			xf_tag = flankseq;

			size_t count = (xf_tag.length() - 5) / 5;

			for (size_t i = 0; i < count; ++i)
			{
				flankstrings.push_back(xf_tag.substr(i * 5 + 5, 4));
			}
		}

		if (xsidx != string::npos)
		{
			sscanf(line.c_str() + xsidx, "%s", genestrand);

			string xs_tag_str = genestrand;

			sscanf(genestrand, "XS:A:%c", &xs_tag);
		}
			//alters = alters_chr;
		strand_t = strand_t & IS_REVERSE;

		strand_t |= IS_PRIMARY;

		if (mis_match)
		{
			size_t mis_idx = line.find("MD:Z:");

			if (mis_idx != string::npos)
			{
				mis_info_str = line.substr(mis_idx + 5);

				//bool was_mis = false;

				size_t current_offset = 0;

				size_t current_num = 0;

				for (size_t i = 0; i < mis_info_str.length(); ++i)
				{
					if (mis_info_str[i] >= '0' && mis_info_str[i] <= '9')
					{
						current_num = current_num * 10 + mis_info_str[i] - '0';

						if (i + 1 < mis_info_str.length() && (mis_info_str[i + 1] < '0' || mis_info_str[i + 1] > '9'))
							current_offset += current_num;
					}
					else
					{
						current_num = 0;

						mis_info.push_back(make_pair(current_offset++, mis_info_str[i]));
					}
				}
				//mis_info;
			}
		}

		if (line.find("ZF:Z:FUS") != string::npos)
			is_fusion_newfmt = true;
		else
			is_fusion_newfmt = false;

		//sscanf(tagname, "%llu", &tagidx);

		tag_name = tagname;

		tag_base_name = tag_name.substr(0, tag_name.length() - 1);

		if (tag_name[tag_name.length() - 1] == '1' || tag_name[tag_name.length() - 1] == 'a')
			end_id = 1;
		else if (tag_name[tag_name.length() - 1] == '2' || tag_name[tag_name.length() - 1] == 'b')
			end_id = 2;

		chrom_name = chrom;
		splice_way = mapped;
		ori_splice_way = mapped;
		mapped_seq = seq;
		qual_str = qual_chr; 
		mate_match = mate_match_chr;

		

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		size_t index = 0;

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		if (isexonic && !is_fusion_newfmt)
		{
			filter_score = 5;

			junc_anchor_len = 500;

			ave_junc_mis = 0;

			junc_hits = 0;

			pair_rate = 5;
		}

		string flag_str = " ";

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
					spliceway_vec.push_back(make_pair(start, -maplen));
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec.push_back(make_pair(start + maplen, 0));
				}
			}
			else if (flag_str[0] == 'M')
			{
				if (flag == 'I')
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
					spliceway_vec.back().second = -maplen;
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec.back().second += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		if (mis_info.size() && mis_info.back().first > mappedlen)
		{
			cout << "mismatch position large than mapped length?"<<endl;
		}
	}
}

void 
SamRec::Clear()
{
	spliceway_vec.clear();

	spliceway_vec2.clear();

	output_line.clear();

	flankstrings.clear();
	//corresponding_juncs.clear();

	xf_tag.clear();

	xs_tag = 0;

	encompassed_juncs.clear();

	left_splice_ways.clear();
	
	right_splice_ways.clear();

	mis_info.clear();

	filter_type = NOT_FILTERED;

	matched_id = -1;
	
	wrong_format = false;
	
	tagidx = -1;
	
	filter_score = 0;
	
	best = 0;

	isunique = false;
	
	isexonic = false;
	
	isspliced = false;
	
	ave_intron_len = 0;
	
	ave_junc_mis = 0;
	
	canon_count = 0;

	noncanon_count = 0;
	
	canon_rate = 1;
	
	issmallins = false;
	
	issmalldel = false;
	
	iscanonical = false;

	issemicanonical = false;
	
	isnoncanoical = false; 
	
	isunmapped = false;
	
	isclipped = false; 
	/*is_insert(false),*/ 
	
	mappedlen = 0;
	
	intron_size = 0;
	
	min_anchor = -1;
	
	max_anchor = 0;
	
	pair_rate = 0;
	
	m_contig_len = 0;

	fusion_prefix_len = 0;
	
	fusion_suffix_len = 0;
	
	is_swapped = false;
	
	mate_offset = 0;
	
	mate_offset2 = 0;
	
	mate_diff = 0;

	mate_diff2 = 0;

	//m = NOT_FILTERED;
	
	mate_match = "*";
	
	mate_match2 = "*";

	paired_type = SINGLE;

	chrom_name.clear();
	splice_way.clear();
	ori_splice_way.clear();

	splice_way2.clear();
	ori_splice_way2.clear();
	chrom_name2.clear();

	start = 0;

	start2 = 0;

	strand_t = 0;
	
	strand_t2 = 0;

	forward_count = 0;
	
	reverse_count = 0;

	is_fusion_newfmt = false;

	is_fusion = false;

	fusion_prefix_paired = false;
	
	fusion_suffix_paired = false;

	fusion_mate_ptr = 0;

	fusion_mate_on_doner_side = false;

	fragment_lengths.clear();

	
}

string 
SamRec::tostring(size_t total, size_t idx) const
{
	//if (!cur_line.empty())
	//	return cur_line;

	char sam_rec_char[5000];

	if (isunmapped)
	{
		string original_mapped_seq = mapped_seq;

		string original_qual_str = qual_str;

		if (strand_t & IS_REVERSE)
		{
			original_mapped_seq = revcomp(original_mapped_seq);

			reverse(original_qual_str.begin(), original_qual_str.end());
		}

		sprintf(sam_rec_char, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tIH:i:0\tHI:i:0", 
			tag_name.c_str(), original_mapped_seq.c_str(), original_qual_str.c_str());
	}
	else if (is_fusion)
		sprintf(sam_rec_char, "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu\tXF:Z:%s", tag_name.c_str(), chrom_name.c_str(), strand1, start, ori_splice_way.c_str(), 
			chrom_name2.c_str(), strand2, start2, ori_splice_way2.c_str(), mapped_seq.c_str(), qual_str.c_str(), mis_match, toflankstring().c_str());
	else
	{
		//cout << "ispliced"<<isspliced << endl;

		//cout << "isunspliced"<<isexonic << endl;
		if (flankstrings.empty())
		{
			//cout << "xf_tag empty"<<endl;

			sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
				chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
			//if (forward_count >= reverse_count)
			
			//else 
			//	sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu\tXS:A:-"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
			//	chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
		}
		else
		{
			//cout << "xf_tag not empty:"<<xf_tag << endl;
			sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tXS:A:%c\tXF:Z:%s\tIH:i:%llu\tHI:i:%llu"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
				chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, xs_tag, toflankstring().c_str(), total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
		}
	}

	return sam_rec_char;
}

string
SamRec::tostandfusion() const 
{
	//R_112:1	115	11	61646784	255	38S37M	1	11184653	0	
	//GCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGTATGAAGGCGGCCAGGGCACTG	181166788356889826893782431399631578517864278199618228323647256272457757777	ZF:Z:FUS_11184689_446432191(++)	ZT:Z:Seed

	//R_112:1	179	1	11184653	255	38M37S	11	61646784	0	
	//GCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGTATGAAGGCGGCCAGGGCACTG	181166788356889826893782431399631578517864278199618228323647256272457757777	ZF:Z:FUS_11184689_446432191(++)	ZT:Z:Seed

	char sam_rec_char[5000];

	if (isunmapped)
	{
		//HWI-EAS217:4:100:1793:1308#0/2  4       *       0       0       *       *       0       0       
		//GTTGAAGTAGGCCGGCACGGTGATCACCGCGTTGGTCACCGGGTAGCTCACCTGCACCTTGGGCTTGTCACCGTCGTTGATCACCTGGAAAGGCCAGTGC    
		//a^`aW```a\^Z_a``aUaYKZZ`P`\T^^WY```\\Z`Za\`TR`[Z\^[Q^^Z^NZPZ\\\[```\\^\R]UZ[T]QKST]RT\^XNTR\ZPTR\T`R IH:i:0  HI:i:0

		string original_mapped_seq;// = mapped_seq;

		string original_qual_str;// = qual_str;

		string original_mapped_seq1;

		string original_mapped_seq2;

		string original_qual_str1;

		string original_qual_str2;

		original_mapped_seq1 = mapped_seq.substr(0, mappedlen1 - insertlen1);

		original_mapped_seq2 = mapped_seq.substr(mapped_seq.length() - mappedlen2 + insertlen2, mappedlen2 - insertlen2);

		original_qual_str1 = qual_str.substr(0, mappedlen1 - insertlen1);

		original_qual_str2 = qual_str.substr(mapped_seq.length() - mappedlen2 + insertlen2, mappedlen2 - insertlen2);

		if (strand_t & IS_REVERSE)
		{
			original_mapped_seq1 = revcomp(original_mapped_seq1);

			reverse(original_qual_str1.begin(), original_qual_str1.end());
		}

		if (strand_t2 & IS_REVERSE)
		{
			original_mapped_seq2 = revcomp(original_mapped_seq2);

			reverse(original_qual_str2.begin(), original_qual_str2.end());
		}

		if (is_swapped)
		{
			original_mapped_seq = original_mapped_seq2 + original_mapped_seq1;

			original_qual_str = original_qual_str2 + original_qual_str1;
		}
		else
		{
			original_mapped_seq =  original_mapped_seq1 + original_mapped_seq2;

			original_qual_str = original_qual_str1 + original_qual_str2;
		}
	
		sprintf(sam_rec_char, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tIH:i:0\tHI:i:0", 
			tag_name.c_str(), original_mapped_seq.c_str(), original_qual_str.c_str());

	}
	else if (is_fusion)
	{
		char append_chr1[1000], append_chr2[1000];

		sprintf(append_chr1, "%lluS", mappedlen2 - insertlen2);

		sprintf(append_chr2, "%lluS", mappedlen1 - insertlen1);

		string sp1, sp2;

		if (is_fusion_newfmt)
		{
			//cout << "is_fusion_newfmt" << endl;
			sp1 = ori_splice_way;

			sp2 = ori_splice_way2;
		}
		else
		{
			if (is_swapped)
			{
				//cout << "is_swapped" << endl;
				sp1 = append_chr1; sp1.append(splice_way);

				sp2 = splice_way2; sp2.append(append_chr2);
			}
			else
			{
				//cout << "is_swapped not" << endl;
				sp1 = splice_way; sp1.append(append_chr1);

				sp2 = append_chr2; sp2.append(splice_way2);
			}
		}

		string combflankstring = toflankstring();

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tXS:A:%c\tXF:Z:%s\tZF:Z:FUS_%llu_%llu(%c%c)\n%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tXS:A:%c\tXF:Z:%s\tZF:Z:FUS_%llu_%llu(%c%c)", 
			tag_name.c_str(), strand_t, chrom_name.c_str(), start, sp1.c_str(), 
			mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, xs_tag, combflankstring.c_str(), fusion_prefix_end, fusion_suffix_st, strand1, strand2,
			tag_name.c_str(), strand_t2, chrom_name2.c_str(), start2, sp2.c_str(), 
			mate_match2.c_str(), mate_offset2, mate_diff2, mapped_seq.c_str(), qual_str.c_str(), mis_match, xs_tag, combflankstring.c_str(), fusion_prefix_end, fusion_suffix_st, strand1, strand2);
	}
	else
		return "";

	return sam_rec_char;
}


void SamRec::genearte_fusion_str()
{
	/*char sam_rec_char[5000];

	if (is_fusion)
	{
		char append_chr1[1000], append_chr2[1000];

		sprintf(append_chr1, "%lluS", mappedlen2);

		sprintf(append_chr2, "%lluS", mappedlen1);

		string sp1, sp2;

		if (is_swapped)
		{
			sp1 = append_chr1; sp1.append(ori_splice_way);

			sp2 = ori_splice_way2; sp2.append(append_chr2);
		}
		else
		{
			sp1 = ori_splice_way; sp1.append(append_chr1);

			sp2 = append_chr2; sp2.append(ori_splice_way2);
		}

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)\n%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)", 
			tag_name.c_str(), strand_t, chrom_name.c_str(), start, sp1.c_str(), 
			mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2,
			tag_name.c_str(), strand_t2, chrom_name2.c_str(), start2, sp2.c_str(), 
			mate_match2.c_str(), mate_offset2, mate_diff2, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2);
	}
	else
		return "";

	return sam_rec_char;*/
}

void SamRec::generate_normal_str()
{
	//char sam_rec_char[5000];

	//if (is_fusion)
	//	sprintf(sam_rec_char, "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tag_name.c_str(), chrom_name.c_str(), strand1, start, ori_splice_way.c_str(), 
	//		chrom_name2.c_str(), strand2, start2, ori_splice_way2.c_str(), mapped_seq.c_str(), qual_str.c_str(), mis_match);
	//else
	//	sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
	//		chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);

	//return sam_rec_char;
}

void SamRec::set_unmapped()
{
	mate_match = "*";

	mate_offset = 0;

	mate_diff = 0;

	mate_match2 = "*";

	mate_offset2 = 0;

	mate_diff2 = 0;
}

bool SamRec::clip_by_small_anchor(bool add_S/*JunctionSeed* filtered_junc*/)
{
	const string& readstr = mapped_seq;

	if (max_anchor >= 50)
	{
		vector<pair<size_t, int> >::iterator vp_iter, max_vp_iter;

		int max_frame_len = 0;
	
		for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		{
			if (vp_iter->second > max_frame_len)
			{
				max_frame_len = vp_iter->second;
				max_vp_iter = vp_iter;
			}
		}

		pair<size_t, int> best = *max_vp_iter;

		if (max_frame_len >= 50)
		{
			xf_tag.clear();

			flankstrings.clear();

			size_t prev_map_len = 0;

			for (vp_iter =  spliceway_vec.begin(); vp_iter != max_vp_iter; ++vp_iter)
			{
				prev_map_len += abs(vp_iter->second);
			}

			if (!add_S)
			{
				qual_str = qual_str.substr(prev_map_len, max_vp_iter->second);

				mapped_seq = mapped_seq.substr(prev_map_len, max_vp_iter->second);
			}

			mappedlen = max_vp_iter->second;

			start = max_vp_iter->first;

			end = max_vp_iter->first + max_vp_iter->second - 1;

			min_anchor = max_vp_iter->second;

			max_anchor = max_vp_iter->second;

			spliceway_vec.clear();

			spliceway_vec.push_back(best);

			char clipped_jmp[100];

			if (true)
			{
				size_t tail_clip = mapped_seq.length() - prev_map_len - best.second;

				if (add_S && tail_clip > mapped_seq.length() )
					cerr <<" tail clip length"<<endl;

				if (!add_S)
					sprintf(clipped_jmp, "%dM", best.second);
				else if (add_S && prev_map_len > 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%lluS%dM%lluS", prev_map_len, best.second, tail_clip);
				else if (add_S && prev_map_len == 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%dM%lluS", best.second, tail_clip);
				else if (add_S && (prev_map_len > 0 && tail_clip == 0))
					sprintf(clipped_jmp, "%lluS%dM", prev_map_len, best.second);
				else 
				{
					cerr <<"clip not correct"<<endl;
					cerr << this->tostring(0, 0)<<endl;
				}
			}

			ori_splice_way = clipped_jmp;
			 
			//unsigned short mis_match;

			canon_rate = 1;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			isexonic = true;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = true;

			canon_count = 0;			
			noncanon_count = 0;

			flankstrings.clear();

			return true;
		}
	}

	return false;
}

bool SamRec::clip_fusion_alignment(bool add_S/*JunctionSeed**/)
{
	const string& readstr = mapped_seq;

	size_t fusion_prefix_len, fusion_suffix_len;

	if (fusion_prefix_len > 50 || fusion_suffix_len > 50)
	{
		if (fusion_prefix_len >= fusion_suffix_len)
		{
			if (!add_S)
			{
				qual_str = qual_str.substr(0, fusion_prefix_len);

				mapped_seq = mapped_seq.substr(0, fusion_prefix_len);
			}

			mappedlen = fusion_prefix_len;

			start = spliceway_vec.front().first;

			end = spliceway_vec.back().first + spliceway_vec.back().second;

			char clipped_jmp[100];

			sprintf(clipped_jmp, "%lluS", fusion_suffix_len);

			isexonic = false;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = false;
			is_fusion = false;
			is_fusion_newfmt = false;

			if ((ori_splice_way.find("I") != string::npos))
				issmallins = true;
			if ((ori_splice_way.find("D") != string::npos))
				issmalldel = true;
			if (ori_splice_way.find("N") != string::npos)
				isspliced = true;

			if (add_S)
				ori_splice_way.append(clipped_jmp);

			//unsigned short mis_match;

			canon_rate = 1;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			canon_count = 0;

			noncanon_count = 0;

			ori_splice_way2.clear();

			spliceway_vec2.clear();

			splice_way2.clear();

			return true;
		}
		else
		{
			spliceway_vec = spliceway_vec2;

			ori_splice_way = ori_splice_way2;

			splice_way = splice_way2;

			strand1 = strand2;

			strand_t = strand_t2;

			mappedlen1 = mappedlen2;

			if (!add_S)
			{
				qual_str = qual_str.substr(fusion_prefix_len + 1, fusion_suffix_len);

				mapped_seq = mapped_seq.substr(fusion_prefix_len + 1, fusion_suffix_len);
			}

			mappedlen = fusion_suffix_len;

			start = spliceway_vec.front().first;

			end = spliceway_vec.back().first + spliceway_vec.back().second;

			char clipped_jmp[100];

			sprintf(clipped_jmp, "%lluS", fusion_prefix_len);

			isexonic = false;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = false;
			is_fusion = false;
			is_fusion_newfmt = false;

			if ((ori_splice_way.find("I") != string::npos))
				issmallins = true;
			if ((ori_splice_way.find("D") != string::npos))
				issmalldel = true;
			if (ori_splice_way.find("N") != string::npos)
				isspliced = true;

			if (add_S)
				ori_splice_way = clipped_jmp + ori_splice_way;

			//unsigned short mis_match;

			canon_rate = 1;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			canon_count = 0;

			noncanon_count = 0;

			ori_splice_way2.clear();

			spliceway_vec2.clear();

			splice_way2.clear();

			return true;
		}
	}

	return false;
}

bool SamRec::clip_by_end_mismatch(bool add_S, size_t max_insert_len/*JunctionSeed* filtered_junc*/)
{
	vector<pair<size_t, char> >::iterator mis_info_iter;

	for (mis_info_iter = mis_info.begin(); mis_info_iter != mis_info.end(); ++mis_info_iter)
	{
		if (mis_info_iter->first >= max_insert_len)
			break;
	}

	string clipped_str;

	if (mis_info_iter != mis_info.begin())
	{
		spliceway_vec.front().first += (mis_info_iter - 1)->first + 1;

		start += (mis_info_iter - 1)->first + 1;

		spliceway_vec.front().second -= (mis_info_iter - 1)->second;

		mis_match = mis_info_iter - mis_info.begin();

		char clipped_chr[1000];

		sprintf(clipped_chr, "%lluI", (mis_info_iter - 1)->first + 1);

		clipped_str.append(clipped_chr);
		//spliceway_vec.pu
	}

	const string& readstr = mapped_seq;

	if (max_anchor >= 50)
	{
		vector<pair<size_t, int> >::iterator vp_iter, max_vp_iter;

		int max_frame_len = 0;
	
		for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		{
			if (vp_iter->second > max_frame_len)
			{
				max_frame_len = vp_iter->second;
				max_vp_iter = vp_iter;
			}
		}

		pair<size_t, int> best = *max_vp_iter;

		if (max_frame_len >= 50)
		{
			size_t prev_map_len = 0;

			for (vp_iter =  spliceway_vec.begin(); vp_iter != max_vp_iter; ++vp_iter)
			{
				prev_map_len += abs(vp_iter->second);
			}

			if (!add_S)
			{
				qual_str = qual_str.substr(prev_map_len, max_vp_iter->second);

				mapped_seq = mapped_seq.substr(prev_map_len, max_vp_iter->second);
			}

			mappedlen = max_vp_iter->second;

			start = max_vp_iter->first;

			end = max_vp_iter->first + max_vp_iter->second - 1;

			min_anchor = max_vp_iter->second;

			max_anchor = max_vp_iter->second;

			spliceway_vec.clear();

			spliceway_vec.push_back(best);

			char clipped_jmp[100];

			if (true)
			{
				size_t tail_clip = mapped_seq.length() - prev_map_len - best.second;

				if (add_S && tail_clip > mapped_seq.length() )
					cerr <<" tail clip length"<<endl;

				if (!add_S)
					sprintf(clipped_jmp, "%dM", best.second);
				else if (add_S && prev_map_len > 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%lluS%dM%lluS", prev_map_len, best.second, tail_clip);
				else if (add_S && prev_map_len == 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%dM%lluS", best.second, tail_clip);
				else if (add_S && (prev_map_len > 0 && tail_clip == 0))
					sprintf(clipped_jmp, "%lluS%dM", prev_map_len, best.second);
				else 
				{
					cerr <<"clip not correct"<<endl;
					cerr << this->tostring(0, 0)<<endl;
				}
			}

			ori_splice_way = clipped_jmp;
			 
			//unsigned short mis_match;

			canon_rate = 1;

			ave_intron_len = 0;

			filter_score = 0;

			ave_junc_mis = 0;

			isexonic = true;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = true;

			canon_count = 0;			
			noncanon_count = 0;

			return true;
		}
	}

	return false;
}

void 
SamRec::setfusionbit()
{
	strand_t |= IS_PAIRED;

	strand_t |= IS_PAIRED_MAPPED;

	if (strand_t2 & IS_REVERSE)
		strand_t |= IS_MATE_REVERSE;

	strand_t |= IS_FIRST_END;

	//if (mate_match == "*")
	//{
	//	if (chrom_name2 == chrom_name)
	//		mate_match = "=";
	//	else
	//		mate_match = chrom_name2;

	//	mate_offset = start2;

	//	mate_diff = fusion_prefix_st - fusion_suffix_end;
	//}

	strand_t2 |= IS_PAIRED;

	strand_t2 |= IS_PAIRED_MAPPED;

	if (strand_t & IS_REVERSE)
		strand_t2 |= IS_MATE_REVERSE;

	strand_t2 |= IS_SECOND_END;

	//if (mate_match2 == "*")
	//{
	//	if (chrom_name2 == chrom_name)
	//		mate_match2 = "=";
	//	else
	//		mate_match2 = chrom_name2;

	//	mate_offset2 = start;

	//	mate_diff2 = fusion_suffix_end - fusion_prefix_st;
	//}

}

bool SamRec::modify_jumpcode_by_filtered_junc()
{
	vector<JunctionSeed*>::iterator junc_iter;

	bool filtered = false;

	if (issmallins)
		filtered = true;

	//for (junc_iter = corresponding_juncs.begin(); !filtered && junc_iter != corresponding_juncs.end(); ++junc_iter)
	//{
	//	//if ((*junc_iter)->m_filtered_type != NOT_FILTERED)
	//	//	filtered = true;

	//	if ((*junc_iter)->m_filtered_type == FILTERED_BY_SMALL_ANCHOR)
	//	{

	//	}

	//	//if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_ANCHOR)
	//	//{
	//	//	if ((*junc_iter)->m_max_prefix_len < m_min_junc_anchor)
	//	//		;
	//	//}

	//	//if ((*sam_rec_ptr_iter)->filter_type.empty())
	//	//	filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);
	//	//else if ((*sam_rec_ptr_iter)->filter_type.size() == 1)
	//	//{
	//	//	if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == NOT_FILTERED)
	//	//		filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);
	//	//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_ANCHOR)
	//	//		;//trim
	//	//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_DELETION)
	//	//		;//convert
	//	//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_LARGE_MIN_ANCHOR_DIFF)
	//	//		;//filter
	//	//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_LARGE_MULTIPLE_PAIRED)
	//	//		;//filter
	//	//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR)
	//	//		;
	//	//}
	//	//else
	//	//{
	//	//}
	//}

	return true;
}

bool SamRec::Swap()
{
	swap(splice_way2, splice_way);

	swap(ori_splice_way2, ori_splice_way);

	swap(chrom_name2, chrom_name);

	swap(splice_way2, splice_way);

	swap(start2, start);

	swap(end2, end);

	swap(strand_t2, strand_t);

	string mapseq1 = mapped_seq.substr(0, mappedlen1);

	string mapseq2 = mapped_seq.substr(mapped_seq.length() - mappedlen2, mappedlen2);

	mapped_seq = mapseq2 + mapseq1;

	string qualstr1 = qual_str.substr(0, mappedlen1);

	string qualstr2 = qual_str.substr(qual_str.length() - mappedlen2, mappedlen2);

	qual_str = qualstr2 + qualstr1;

	swap(mappedlen1, mappedlen2);

	swap(insertlen1, insertlen2);

	swap(fusion_prefix_len, fusion_suffix_len);

	swap(spliceway_vec2, spliceway_vec);

	swap(mappedlen1, mappedlen2);

	swap(strand1, strand2);

	reverse(flankstrings.begin(), flankstrings.end());

	if (xs_tag == '+')
		xs_tag = '-';
	else if (xs_tag == '-')
		xs_tag = '+';

	for (size_t i = 0; i < flankstrings.size();++i)
	{
		flankstrings[i] = revcomp(flankstrings[i]);
	}

	if (strand1 == '+')
		strand1 = '-';
	else
		strand1 = '+';

	if (strand2 == '+')
		strand2 = '-';
	else
		strand2 = '+';

	if (strand1 == '+')
		strand_t = 0;
	else
		strand_t = 16;

	if (strand2 == '+')
		strand_t2 = 0;
	else
		strand_t2 = 16;

	strand_t |= IS_PRIMARY;

	strand_t2 |= IS_PRIMARY;

	if (strand1 == '+')
	{
		fusion_prefix_st = start;

		fusion_prefix_end = end;
	}
	else
	{
		fusion_prefix_st = end;

		fusion_prefix_end = start;
	}

	if (strand2 == '+')
	{
		fusion_suffix_st = start2;

		fusion_suffix_end = end2;
	}
	else
	{
		fusion_suffix_st = end2;

		fusion_suffix_end = start2;
	}

	is_swapped = true;

	return true;
}

bool SamRec::is_in_region(vector<size_t>& sorted_region)
{
	vector<pair<size_t, int> >::iterator splice_iter;

	for (splice_iter = spliceway_vec.begin(); splice_iter != spliceway_vec.end(); ++splice_iter)
	{
		size_t cur_pos = splice_iter->first;

		vector<size_t>::iterator lower_bound_iter = lower_bound(sorted_region.begin(), sorted_region.end(), cur_pos);

		if ((lower_bound_iter == sorted_region.end()) || (*lower_bound_iter != cur_pos && ((lower_bound_iter - sorted_region.begin()) % 2) == 0))
			return false;

		if (splice_iter->second > 0)
		{
			cur_pos = splice_iter->first + splice_iter->second - 1;

			lower_bound_iter = lower_bound(sorted_region.begin(), sorted_region.end(), cur_pos);

			if ((lower_bound_iter == sorted_region.end()) || (*lower_bound_iter != cur_pos && ((lower_bound_iter - sorted_region.begin()) % 2) == 0))
				return false;
		}
	}

	return true;
}

bool SamRec::generate_endpoint_expressison(vector<size_t>& sorted_region, size_t offset, vector<vector<int> >& exon_express)
{
	#ifdef DEBUG
	cout << "generate exon expression"<<endl;
	#endif

	//vector<pair<size_t, int> >::iterator splice_iter;

	size_t cur_pos = offset;

	vector<size_t>::iterator lower_bound_iter = lower_bound(sorted_region.begin(), sorted_region.end(), cur_pos);

	#ifdef DEBUG

	for (size_t i = 0; i < sorted_region.size(); ++i)
	{
		cout << sorted_region[i] << '\t' ;
	}
 
	cout << endl;

	cout << "lower_bound_iter:"<<*lower_bound_iter<<endl;

	#endif

	if ((lower_bound_iter == sorted_region.end()) || (*lower_bound_iter != cur_pos && ((lower_bound_iter - sorted_region.begin()) % 2) == 0))
		return false;
	else
	{
		#ifdef DEBUG

		cout << "exon_start"<<endl;

		#endif

		size_t exon_start;
		if (*lower_bound_iter == cur_pos && (((lower_bound_iter - sorted_region.begin()) % 2) == 0))
			exon_start = *lower_bound_iter;
		else
			exon_start = *(lower_bound_iter - 1);

		#ifdef DEBUG

		cout << "exon_start:" << exon_start << endl;

		cout << "offset:" << offset << endl;

		cout << "offset - exon_start:" << offset - exon_start << endl;

		cout << "exon index:"<< (lower_bound_iter - sorted_region.begin()) / 2 << endl;

		cout << exon_express[(lower_bound_iter - sorted_region.begin()) / 2].size() << endl;

		#endif

		if (offset - exon_start >  exon_express[(lower_bound_iter - sorted_region.begin()) / 2].size())
		{
			cerr << "index exceed"<<endl;

			exit(1);
		}
		++(exon_express[(lower_bound_iter - sorted_region.begin()) / 2][offset - exon_start]);
	}

	cout << "generate exon expression finished"<<endl;

	return true;
}
	

string SamRec::toflankstring() const
{
	string comb_flankstring;

	for (size_t i = 0; i < flankstrings.size(); ++i)
		comb_flankstring += flankstrings[i] + ",";

	return comb_flankstring;
}