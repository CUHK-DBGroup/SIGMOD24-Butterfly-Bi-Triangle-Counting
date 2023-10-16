#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <time.h>
#include <cstdio>
#include <cassert>
#include <cstdio>
#include <stdio.h>
#include <random>
#include <numeric>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <map>
#include <sstream>
#include <fstream>
//#include <unistd.h>
//#include <boost\random\mersenne_twister.hpp>
//#include <boost\random\uniform_int_distribution.hpp>
//#include <boost\random\uniform_real_distribution.hpp>
//#include <boost\random\random_device.hpp>

using namespace std;
//using namespace boost::random;

#define SZ(x) ((int)x.size())
#define ll long long
#define ull unsigned long long
#define ld long double
#define eps 1e-11
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

string alias = "";

string input_address, output_address, bitrigt_address, bfcgt_address;

set < pair <int, int> > edges;
vector < pair <int, int> > list_of_edges;
unordered_map < int, int > vertices[2];
vector <int> index_map;
vector <int> index_map2;

vector <int> vertices_in_left;
vector <int> vertices_in_right;
vector < vector <int> > adj;
vector < vector < int > > sampled_adj_list;
vector <bool> visited;
vector <int> list_of_vertices;
vector <int> vertex_counter;
vector <int> tot_deg;
vector<int> weight_pair_id;

bool swap_flag = 0;
ll n_vertices;
ll n_edges;
//ll 
ll n_wedge_in_partition[2];
ld exact_n_bf;
ll largest_index_in_partition[2];

vector <int> clr;
vector <int> hashmap_C;
vector <ll> sum_wedges;
vector <ll> sum_deg_neighbors;
vector <int> aux_array_two_neighboorhood;

vector<int> random_label;
void clear_everything() {
	largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
	n_vertices = 0;
	n_edges = 0;
	edges.clear();
	vertices[0].clear(); vertices[1].clear();
	index_map.clear();
	vertices_in_left.clear();
	vertices_in_right.clear();
	adj.clear();
	sampled_adj_list.clear();
	visited.clear();
	list_of_edges.clear();
	vertex_counter.clear();
	clr.clear();
	hashmap_C.clear();
	sum_wedges.clear();
	sum_deg_neighbors.clear();
	aux_array_two_neighboorhood.clear();
}

vector<int> LP_label_last;
vector<int> LP_label_now;
void resize_all() {
	clr.resize(n_vertices);
	hashmap_C.resize(n_vertices);
	random_label.resize(n_vertices);
	aux_array_two_neighboorhood.resize(n_vertices);
	sum_wedges.resize(n_vertices);
	visited.resize(n_vertices);
	index_map.resize(n_vertices);
	index_map2.resize(n_vertices);
	sum_deg_neighbors.resize(n_vertices);
	tot_deg.resize(n_vertices);
	LP_label_last.resize(n_vertices + 5);
	LP_label_now.resize(n_vertices + 5);
	weight_pair_id.resize(n_edges + 1);
}

// ------------- Read the graph ---------------------
void add_vertex(int A, int side) {
	if (vertices[side].find(A) == vertices[side].end()) {
		if (side == 0) vertices_in_left.push_back(A);
		else vertices_in_right.push_back(A);
		vertices[side][A] = 1;
	}
}

void get_index(int& A, int side) {
	if (vertices[side].find(A) == vertices[side].end()) {
		vertices[side][A] = largest_index_in_partition[side] ++;
	}
	A = vertices[side][A];
}

void add_edge(int& A, int& B) {
	add_vertex(A, 0);
	add_vertex(B, 1);
	if (edges.find(make_pair(A, B)) == edges.end()) {
		edges.insert(make_pair(A, B));
		n_edges++;
	}
}

bool all_num(string& s) {
	for (int i = 0; i < SZ(s); i++) if ((s[i] >= '0' && s[i] <= '9') == false) return false;
	return true;
}

void clean_graph() {
	freopen(input_address.c_str(), "r", stdin);
	freopen(output_address.c_str(), "w", stdout);

	//output_address
	cerr << input_address << "\n";
	string s;
	cin.clear();
	while (getline(cin, s)) {
		stringstream ss; ss << s;
		vector <string> vec_str;
		for (string z; ss >> z; vec_str.push_back(z));
		if (SZ(vec_str) >= 2) {
			bool is_all_num = true;
			for (int i = 0; i < min(2, SZ(vec_str)); i++) is_all_num &= all_num(vec_str[i]);
			if (is_all_num) {
				int A, B;
				ss.clear(); ss << vec_str[0]; ss >> A;
				ss.clear(); ss << vec_str[1]; ss >> B;
				if (swap_flag) {
					swap(A, B);
				}
				add_edge(A, B);
			}
		}
	}
	vertices[0].clear();
	vertices[1].clear();
	int lv_num, rv_num;
	largest_index_in_partition[0] = 0;
	lv_num = largest_index_in_partition[1] = SZ(vertices_in_left);
	n_vertices = SZ(vertices_in_left) + SZ(vertices_in_right);
	adj.resize(n_vertices, vector <int>());

	for (auto edge : edges) {
		int A = edge.first;
		int B = edge.second;
		get_index(A, 0);
		get_index(B, 1);
		cout << A << " " << B - lv_num << "\n";
		adj[A].push_back(B);
		adj[B].push_back(A);
		list_of_edges.push_back(make_pair(A, B));
	}
	resize_all();

	cout.flush();
	fclose(stdout);
	n_wedge_in_partition[0] = 0;
	for (int i = 0; i < largest_index_in_partition[0]; i++) {
		n_wedge_in_partition[0] += (((ll)SZ(adj[i])) * (SZ(adj[i]) - 1)) >> 1;
	}
	n_wedge_in_partition[1] = 0;
	for (int i = largest_index_in_partition[0]; i < largest_index_in_partition[1]; i++) {
		n_wedge_in_partition[1] += ((ll)SZ(adj[i]) * (SZ(adj[i]) - 1)) >> 1;
	}
	for (int i = 0; i < n_vertices; i++) {
		sort(adj[i].begin(), adj[i].end());
		sum_deg_neighbors[i] = 0;
		for (auto neighbor : adj[i]) {
			sum_deg_neighbors[i] += SZ(adj[neighbor]);
		}
	}
	cerr << " for test # edges :: " << SZ(list_of_edges) << " left :: " << SZ(vertices_in_left) << " right :: " << SZ(vertices_in_right) << endl;
	sort(list_of_edges.begin(), list_of_edges.end());
	cerr << n_wedge_in_partition[0] << " " << n_wedge_in_partition[1] << "\n";

	fclose(stdin);
}

void get_graph() {
	freopen(input_address.c_str(), "r", stdin);
	//freopen(output_address.c_str(), "w", stdout);

	//output_address
	cerr << input_address << "\n";
	string s;
	cin.clear();
	while (getline(cin, s)) {
		stringstream ss; ss << s;
		vector <string> vec_str;
		for (string z; ss >> z; vec_str.push_back(z));
		if (SZ(vec_str) >= 2) {
			bool is_all_num = true;
			for (int i = 0; i < min(2, SZ(vec_str)); i++) is_all_num &= all_num(vec_str[i]);
			if (is_all_num) {
				int A, B;
				ss.clear(); ss << vec_str[0]; ss >> A;
				ss.clear(); ss << vec_str[1]; ss >> B;
				if (swap_flag) {
					swap(A, B);
				}
				add_edge(A, B);
			}
		}
	}
	vertices[0].clear();
	vertices[1].clear();
	int lv_num, rv_num;
	largest_index_in_partition[0] = 0;
	lv_num = largest_index_in_partition[1] = SZ(vertices_in_left);
	n_vertices = SZ(vertices_in_left) + SZ(vertices_in_right);
	adj.resize(n_vertices, vector <int>());

	for (auto edge : edges) {
		int A = edge.first;
		int B = edge.second;
		get_index(A, 0);
		get_index(B, 1);
		//cout << A << " " << B - lv_num << "\n";
		adj[A].push_back(B);
		adj[B].push_back(A);
		list_of_edges.push_back(make_pair(A, B));
	}
	resize_all();
	//cerr << "finish format" << endl;
	//cerr << edges.size() << "\n";
	//cout.flush();
	//fclose(stdout);
	n_wedge_in_partition[0] = 0;
	ld deg1000 = 0, deg5000 = 0;

	for (int i = 0; i < largest_index_in_partition[0]; i++) {
		if (SZ(adj[i]) >= 1000) deg1000 += SZ(adj[i]);
		if (SZ(adj[i]) >= 5000) deg5000 += SZ(adj[i]);

		n_wedge_in_partition[0] += (((ll)SZ(adj[i])) * (SZ(adj[i]) - 1)) >> 1;
	}
	n_wedge_in_partition[1] = 0;
	for (int i = largest_index_in_partition[0]; i < largest_index_in_partition[1]; i++) {
		n_wedge_in_partition[1] += ((ll)SZ(adj[i]) * (SZ(adj[i]) - 1)) >> 1;
	}
	for (int i = 0; i < n_vertices; i++) {
		sort(adj[i].begin(), adj[i].end());
		sum_deg_neighbors[i] = 0;
		for (auto neighbor : adj[i]) {
			sum_deg_neighbors[i] += SZ(adj[neighbor]);
		}
	}
	cerr << " for test # edges :: " << SZ(list_of_edges) << " left :: " << SZ(vertices_in_left) << " right :: " << SZ(vertices_in_right) << endl;
	sort(list_of_edges.begin(), list_of_edges.end());
	cerr << n_wedge_in_partition[0] << " " << n_wedge_in_partition[1] << "\n";
	//cerr << "propotion:" << deg1000 / SZ(list_of_edges) << " " << deg5000 / SZ(list_of_edges) << "\n";

	fclose(stdin);
}
random_device rd_edge;
mt19937_64 eng_ran_bfc_per_edge(rd_edge());

void read_the_graph(string alias) {
	clear_everything();
	cerr << " Insert the input (bipartite network) file location" << endl;
	input_address = "/home/usr/bgdata/" + alias + "/out." + alias;
	cerr << " Insert the output file" << endl;
	output_address = "/home/usr/bgdata/" + alias + "/out2." + alias;

	get_graph();
	cerr << " -------------------------------------------------------------------------- \n";
	cerr << " The graph is processed - there are " << n_vertices << " vertices and " << n_edges << " edges  \n";
	cerr << " -------------------------------------------------------------------------- \n";
}

double Modularity_graph(vector<int>& Label) {
	double it1 = 0;
	for (auto j : list_of_edges) {
		if (Label[j.first] == Label[j.second]) {
			it1 += 1;
		}
	}
	unordered_map<int, double> kt, dt;
	for (int i = 0; i < largest_index_in_partition[0]; i++) {
		if (kt.find(Label[i]) == kt.end())
			kt[Label[i]] = 0;
		kt[Label[i]] += SZ(adj[i]);
	}
	for (int i = largest_index_in_partition[0]; i < largest_index_in_partition[1]; i++) {
		if (dt.find(Label[i]) == dt.end())
			dt[Label[i]] = 0;
		dt[Label[i]] += SZ(adj[i]);
	}
	double it2 = 0;
	for (auto i : kt) {
		it2 += kt[i.first] * dt[i.first] / n_edges;
	}
	return (it1 - it2) / n_edges;
}

int cmp_large(pair<int, int> a, pair<int, int> b) {
	return a.first > b.first;
}

void LP_algorithm() {
	for (int i = 0; i < n_vertices; i++) {
		LP_label_last[i] = i;
	}
	random_device rdw;
	mt19937 genedg(rdw());
	int rnd = 0;
	double mod_last = 0, mod_new = 0;
	while (rnd <= 10000) {
		rnd++;
		bool flag = 0;
		unordered_map<int, int> i_label;
		vector<int> i_max_label;
		int mx_label_num = 0;
		for (int i = 0; i < n_vertices; i++) {
			for (auto j : adj[i]) {
				if (i_label.find(LP_label_last[j]) == i_label.end()) {
					i_label[LP_label_last[j]] = 0;
				}
				i_label[LP_label_last[j]] += 1;
				mx_label_num = max(mx_label_num, i_label[LP_label_last[j]]);
			}
			for (auto j : i_label) {
				if (j.second == mx_label_num)
					i_max_label.push_back(j.first);
			}
			uniform_int_distribution<int> dis(0, i_max_label.size() - 1);
			int new_label = i_max_label[dis(genedg)];
			if (new_label != LP_label_last[i]) {
				flag = 1;
				LP_label_last[i] = new_label;
			}
			i_max_label.clear();
			mx_label_num = 0;
			i_label.clear();
		}
		mod_new = Modularity_graph(LP_label_last);
		if (abs(mod_last - mod_new) <= 0.001 || mod_last - mod_new > 0.001) break;
		mod_last = mod_new;
		if (!flag) break;
	}

	unordered_map<int, int> tot_label;
	int tot_cluster = 0;
	for (int i = 0; i < n_vertices; i++) {
		if (tot_label.find(LP_label_last[i]) == tot_label.end()) {
			tot_label[LP_label_last[i]] = tot_cluster;
			tot_cluster++;
		}
	}

	vector<pair<int, int>> cluster_size;
	cluster_size.resize(tot_cluster + 5);
	for (int i = 0; i < tot_cluster; i++) {
		cluster_size[i] = make_pair(0, i);
	}
	for (int i = 0; i < n_vertices; i++) {
		cluster_size[tot_label[LP_label_last[i]]].first++;
	}
	sort(cluster_size.begin(), cluster_size.begin() + tot_cluster, cmp_large);
	int rnk = 0;
	for (auto i : cluster_size) {
		tot_label[i.second] = rnk;
	}

	for (int i = 0; i < n_vertices; i++) {
		LP_label_last[i] = tot_label[LP_label_last[i]];
	}

	int PRINT_CLUSTER = 10;
	unordered_map<int, int> big_cluster_map;
	for (int i = tot_cluster - 1; i >= tot_cluster - 1 - PRINT_CLUSTER; i--) {
		big_cluster_map[cluster_size[i].second] = 1;
	}
	vector<vector<pair<int, int>>> cls_edges;
	cls_edges.resize(tot_cluster + 5);
	for (auto j : list_of_edges) {
		if (LP_label_last[j.first] == LP_label_last[j.second]) {
			cls_edges[LP_label_last[j.first]].push_back(j);
		}
	}
	return;
}

void init_deg(vector<int>& deg_, vector<int>& Label, int side, int cluster_size) {
	for (int i = 0; i < cluster_size; i++) deg_[i] = 0;
	for (int i = (side == 0 ? 0 : largest_index_in_partition[0]); i <= largest_index_in_partition[side]; i++) {
		//cerr << i << " " << Label[i] << " " << SZ(adj[i]) << "\n";
		deg_[Label[i]] += SZ(adj[i]);
	}
}

void print_cluster(vector<int>& Label, int cluster_size, string cluster_name) {
	vector<vector<pair<int, int>>> cls_edges;
	vector<int> cls_size;
	cls_edges.resize(cluster_size + 5);
	for (auto j : list_of_edges) {
		if (Label[j.first] == Label[j.second]) {
			cls_edges[Label[j.first]].push_back(j);
		}
	}
	for (int i = 0; i < cluster_size; i++) {
		string cluster_graph = "/home/usr/bgdata/" + alias + "/" + alias + "_" + cluster_name + "/out." + to_string(i) + alias;
		std::ofstream cls_out(cluster_graph);
		for (auto j : cls_edges[i]) {
			cls_out << j.first << " " << j.second << "\n";
		}
		cls_out << endl;
		cls_out.close();
	}
}

void BRIM(int cluster_size, int init_cluster, int cluster_method) {
	//int cluster_size;
	vector<int> Label;
	Label.resize(n_vertices);
	if (init_cluster == 0) {
		random_device rdw;
		mt19937 genedg(rdw());
		uniform_int_distribution<int> dis(0, cluster_size - 1);
		for (int i = 0; i < n_vertices; i++) {
			Label[i] = i % cluster_size;
		}
	}
	else {
		for (int i = 0; i < n_vertices; i++) {
			Label[i] = LP_label_last[i] >= cluster_size ? (LP_label_last[i] % cluster_size) : LP_label_last[i];
		}
	}
	vector<int> deg_R, deg_L;
	deg_L.resize(cluster_size);
	deg_R.resize(cluster_size);
	init_deg(deg_L, Label, 0, cluster_size);
	init_deg(deg_R, Label, 1, cluster_size);
	double mod_new = 0;
	double mod_last = 0;
	mod_last = Modularity_graph(Label);

	while (true) {
		for (int i = 0; i < largest_index_in_partition[0]; i++) {
			vector<int> intra_deg;
			intra_deg.resize(cluster_size);
			for (auto j : adj[i]) {
				intra_deg[Label[j]] += 1;
			}
			double mx = -1000000;
			int new_c = 0;
			for (int c = 0; c < cluster_size; c++) {
				double mod_i = intra_deg[c] - double(SZ(adj[i])) * deg_R[c] / n_edges;

				if (mx < mod_i) {
					mx = mod_i;
					new_c = c;
				}
			}
			Label[i] = new_c;
			intra_deg.clear();
		}
		init_deg(deg_L, Label, 0, cluster_size);

		for (int i = largest_index_in_partition[0]; i < largest_index_in_partition[1]; i++) {
			vector<int> intra_deg;
			intra_deg.resize(cluster_size);
			for (auto j : adj[i]) {
				intra_deg[Label[j]] += 1;
			}
			double mx = -1000000;
			int new_c = 0;
			for (int c = 0; c < cluster_size; c++) {
				double mod_i = intra_deg[c] - double(SZ(adj[i])) * deg_L[c] / n_edges;
				if (mx < mod_i) {
					mx = mod_i;
					new_c = c;
				}
			}
			Label[i] = new_c;
			intra_deg.clear();
		}
		init_deg(deg_R, Label, 1, cluster_size);
		mod_new = Modularity_graph(Label);
		if (abs(mod_last - mod_new) <= 0.0001) break;
		mod_last = mod_new;
	}

	vector<vector<pair<int, int>>> cls_edges;
	vector<int> cls_size;
	cls_edges.resize(cluster_size + 5);
	cls_size.resize(cluster_size + 5);
	for (int i = 0; i < n_vertices; i++) {
		cls_size[Label[i]]++;
	}
	print_cluster(Label, cluster_size, "BRIM"+to_string(cluster_method));
}

int main(int argc, char* args[]) {
	std::ios::sync_with_stdio(false);
	int file = 7;
	int cnt = 0;
	vector<ld> total_er;
	swap_flag = false;

	while (cnt < argc) {
		if (strcmp(args[cnt], "-file") == 0) {
			file = atoi(args[++cnt]);
		}
		else if (strcmp(args[cnt], "-swap") == 0) {
			swap_flag = 1;
		}
		cnt++;
		//if (strcmp(ags[cnt], ""))
	}
	// set graph file name
	string graph_in = "/home/usr/bgdata/" + alias + "/out." + alias;
	read_the_graph(alias);

	clock_t start, end;
	start = clock(); 

	LP_algorithm();
	end = clock();  
	cout << "LP time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  
	BRIM(10,0, 3);

	end = clock();
	cout << "BLP time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cerr << " Take a look at the output file ..." << endl;
	return 0;
}