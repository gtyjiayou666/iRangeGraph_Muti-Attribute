#pragma once

#include <vector>
#include "utils.h"
#include "searcher.hpp"
#include "memory.hpp"
#include <bitset>
#include <stack>
#include <cmath>
#include <climits>

namespace iRangeGraph
{
    typedef std::pair<float, std::pair<int, int>> PFII;
    template <typename dist_t>
    class iRangeGraph_Search_Muti_Interval
    {
    public:
        DataLoaderMultiInterval *storage;
        SegmentTree *tree;
        size_t max_elements_{0};
        size_t dim_{0};
        size_t M_out{0};
        size_t ef_construction{0};

        size_t size_data_per_element_{0};
        size_t size_links_per_element_{0};
        size_t data_size_{0};

        size_t size_links_per_layer_{0};
        size_t offsetData_{0};

        char *data_memory_{nullptr};

        hnswlib::L2Space *space;
        hnswlib::DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_{nullptr};

        size_t metric_distance_computations{0};
        size_t metric_hops{0};

        size_t metric_distance_computations1{0};
        size_t metric_hops1{0};
        int prefetch_lines{0};

        iRangeGraph_Search_Muti_Interval(std::string vectorfilename, std::string edgefilename, DataLoaderMultiInterval *store, int M) : storage(store)
        {
            std::ifstream vectorfile(vectorfilename, std::ios::in | std::ios::binary);
            if (!vectorfile.is_open())
                throw Exception("cannot open " + vectorfilename);
            std::ifstream edgefile(edgefilename, std::ios::in | std::ios::binary);
            if (!edgefile.is_open())
                throw Exception("cannot open " + edgefilename);

            vectorfile.read((char *)&max_elements_, sizeof(int));
            vectorfile.read((char *)&dim_, sizeof(int));

            tree = new SegmentTree(max_elements_);
            tree->BuildTree(tree->root);

            space = new hnswlib::L2Space(dim_);
            fstdistfunc_ = space->get_dist_func();
            dist_func_param_ = space->get_dist_func_param();
            M_out = M;

            data_size_ = (dim_ + 7) / 8 * 8 * sizeof(float);
            size_links_per_layer_ = M_out * sizeof(tableint) + sizeof(linklistsizeint);
            size_links_per_element_ = (size_links_per_layer_ * (tree->max_depth + 1) + 31) / 32 * 32;
            size_data_per_element_ = size_links_per_element_ + data_size_;
            offsetData_ = size_links_per_element_;
            prefetch_lines = data_size_ >> 4;

            data_memory_ = (char *)memory::align_mm<1 << 21>(max_elements_ * size_data_per_element_);
            if (data_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            for (int pid = 0; pid < max_elements_; pid++)
            {
                for (int layer = 0; layer <= tree->max_depth; layer++)
                {
                    linklistsizeint *data = get_linklist(pid, layer);
                    edgefile.read((char *)data, sizeof(tableint));
                    int size = getListCount(data);
                    if (size > M_out)
                        throw Exception("real linklist size is bigger than defined M_out");
                    for (int i = 0; i < size; i++)
                    {
                        char *current_neighbor_ = (char *)(data + 1 + i);
                        edgefile.read(current_neighbor_, sizeof(tableint));
                    }
                    // std::sort(data + 1, data + 1 + size);
                }

                char *data = getDataByInternalId(pid);
                // vectorfile.read(data, data_size_);
                vectorfile.read(data, dim_ * sizeof(float));
            }

            edgefile.close();
            vectorfile.close();
            std::cout << "load index finished ..." << std::endl;
        }

        ~iRangeGraph_Search_Muti_Interval()
        {
            free(data_memory_);
            data_memory_ = nullptr;
        }

        inline char *getDataByInternalId(tableint internal_id) const
        {
            return (data_memory_ + internal_id * size_data_per_element_ + offsetData_);
        }

        linklistsizeint *get_linklist(tableint internal_id, int layer) const
        {
            return (linklistsizeint *)(data_memory_ + internal_id * size_data_per_element_ + layer * size_links_per_layer_);
        }

        int getListCount(linklistsizeint *ptr) const
        {
            return *((int *)ptr);
        }

        int GetOverLap(int l, int r, int ql, int qr)
        {
            int L = std::max(l, ql);
            int R = std::min(r, qr);
            return R - L + 1;
        }

        std::vector<tableint> SelectEdge(int pid, int ql, int qr, int edge_limit, searcher::Bitset<uint64_t> &visited_set)
        {
            TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<tableint> selected_edges;
            selected_edges.reserve(edge_limit);
            do
            {
                cur_node = nxt_node;
                bool contain = false;
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, ql, qr) == GetOverLap(nxt_node->lbound, nxt_node->rbound, ql, qr))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);

                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);

                for (size_t j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    // if (visitedpool[neighborId] == visited_tag)
                    //     continue;
                    if (visited_set.get(neighborId))
                        continue;
                    selected_edges.emplace_back(neighborId);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }

            } while (cur_node->lbound < ql || cur_node->rbound > qr);

            return selected_edges;
        }

        int MaxStep{20};
        std::vector<double> probability;
        void setprob()
        {
            probability.resize(MaxStep);

            for (int x = 0; x < MaxStep; x++)
            {
                probability[x] = 1 / (1 + exp(x));
            }
        }

        // purepost = True -> p=1   purepost =  False -> 0<=p<=1
        bool purepost{false};

        int ProbFunc(int x)
        {
            if (purepost)
                return 1;
            if (x >= MaxStep)
                return 0;
            double randNum = (double)rand() / RAND_MAX;
            return randNum < probability[x] ? 1 : 0;
        }
        std::vector<bool> flags = std::vector<bool>(10000000, false);
        std::vector<std::pair<tableint, bool>> SelectEdgeF(int pid, int ql, int qr, int edge_limit, int current_step)
        {
            iRangeGraph::TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<std::pair<tableint, bool>> selected_edges;
            do
            {
                cur_node = nxt_node;
                bool contain = false;
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); i++)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, ql, qr) == GetOverLap(nxt_node->lbound, nxt_node->rbound, ql, qr))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);

                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);

                for (size_t j = 1; j <= size; j++)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    int prob = 1;
                    int next_step = current_step + 1;
                    bool inrange = flags[neighborId];
                    if (!inrange)
                        prob = ProbFunc(next_step);
                    if (!prob)
                        continue;

                    selected_edges.emplace_back(neighborId, inrange);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }
            } while (cur_node->lbound < ql || cur_node->rbound > qr);
            return selected_edges;
        }

        std::priority_queue<PFI> TopDown_nodeentries_search(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFI, std::vector<PFI>, std::greater<PFI>> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(dis, pid);
                top_candidates.emplace(dis, pid);
            }

            float lowerBound = top_candidates.top().first;

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second;
                auto selected_edges = SelectEdge(current_pid, QL, QR, edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];

                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        std::priority_queue<PFI> MIDG(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<int> &QL, std::vector<int> &QR, int edge_limit, int id)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);
            std::priority_queue<PFIMI, std::vector<PFIMI>, Compare> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                    continue;
                visited_set.set(pid);
                ++metric_hops;
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                ++metric_distance_computations;
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid, id);
                    top_candidates.emplace(dis, pid);
                    lowerBound = top_candidates.top().first;
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid, id);
                    top_candidates.emplace(dis, pid);
                    top_candidates.pop();
                    lowerBound = top_candidates.top().first;
                }
            }
            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.dis > lowerBound)
                    break;
                candidate_set.pop();
                int current_pid = current_point_pair.id;
                int range_id = current_point_pair.range_id;
                auto selected_edges = SelectEdge(current_pid, QL[range_id], QR[range_id], edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, current_point_pair.range_id);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, current_point_pair.range_id);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
                std::vector<std::pair<tableint, int>> edges;
                edges.reserve(edge_limit / 2);
                JINS(QL, QR, current_pid, edge_limit / 2, 2, edges);
                num_edges = edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(edges[i].first), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = edges[i].first;
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;
                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, edges[i].second);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, edges[i].second);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        std::priority_queue<PFI> MIDG_P(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<int> &QL, std::vector<int> &QR, int edge_limit, int id)
        {
            int interval_num = QL.size();
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);
            std::priority_queue<PFIMI, std::vector<PFIMI>, Compare> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                    continue;
                visited_set.set(pid);
                ++metric_hops;
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                ++metric_distance_computations;
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid, id);
                    top_candidates.emplace(dis, pid);
                    lowerBound = top_candidates.top().first;
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid, id);
                    top_candidates.emplace(dis, pid);
                    top_candidates.pop();
                    lowerBound = top_candidates.top().first;
                }
            }
            bool erfen = interval_num < 16;
            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.dis > lowerBound)
                    break;
                candidate_set.pop();
                int current_pid = current_point_pair.id;
                int range_id = current_point_pair.range_id;
                if (range_id == -1)
                {
                    if (erfen)
                    {
                        for (int i = 0; i < interval_num; i++)
                        {
                            if (current_pid >= QL[i])
                            {
                                if (current_pid <= QR[i])
                                {
                                    range_id = i;
                                    break;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        auto it = std::lower_bound(QL.begin(), QL.end(), current_pid);
                        if (it != QL.end())
                        {
                            int i = static_cast<int>(std::distance(QL.begin(), it));
                            if (*it == current_pid)
                            {
                                range_id = i;
                            }
                            else
                            {
                                range_id = i - 1;
                            }
                        }
                        else
                        {
                            range_id = interval_num - 1;
                        }
                    }
                }
                auto selected_edges = SelectEdge(current_pid, QL[range_id], QR[range_id], edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
                std::vector<tableint> edges;
                edges.reserve(edge_limit / 2);
                JINS_P(QL, QR, current_pid, edge_limit / 2, 2, edges);
                num_edges = edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;
                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        struct GapInfo
        {
            int index;
            int len;
        };
        struct CompareGap
        {
            bool operator()(const GapInfo &a, const GapInfo &b)
            {
                return a.len < b.len;
            }
        };
        void merge_small_gaps(const std::vector<int> &Q_L, const std::vector<int> &Q_R,
                              std::vector<int> &Q_L1, std::vector<int> &Q_R1,
                              std::vector<int> &mapping, double p)
        {
            int n = Q_L.size();
            mapping.assign(n, 0);

            double total_data = Q_R[n - 1] - Q_L[0] + 1;
            double total_gaps = 0;

            std::priority_queue<GapInfo, std::vector<GapInfo>, CompareGap> max_heap_gaps;

            for (int i = 0; i < n - 1; ++i)
            {
                int g_l = Q_R[i] + 1;
                int g_r = Q_L[i + 1] - 1;
                max_heap_gaps.push({i, g_r - g_l + 1});
                total_gaps = total_gaps + g_r - g_l + 1;
            }
            int total_nogaps = total_data - total_gaps;

            std::vector<bool> is_keeper(n - 1, false);

            int s = total_data - (p * total_nogaps);
            if (s <= 0)
            {
                Q_L1.push_back(Q_L[0]);
                Q_R1.push_back(Q_R[n - 1]);
                return;
            }
            int current_gaps = 0;
            while (!max_heap_gaps.empty())
            {
                GapInfo gap = max_heap_gaps.top();
                current_gaps = current_gaps + gap.len;
                if (current_gaps >= s)
                {
                    is_keeper[gap.index] = true;
                    break;
                }
                max_heap_gaps.pop();
                is_keeper[gap.index] = true;
            }
            int current_l = Q_L[0];
            int current_r = Q_R[0];

            int new_idx = 0;
            mapping[0] = new_idx;
            for (int i = 0; i < n - 1; ++i)
            {
                if (is_keeper[i])
                {
                    Q_L1.push_back(current_l);
                    Q_R1.push_back(current_r);
                    new_idx++;
                    current_l = Q_L[i + 1];
                    current_r = Q_R[i + 1];
                    mapping[i + 1] = new_idx;
                }
                else
                {
                    current_r = Q_R[i + 1];
                    mapping[i + 1] = new_idx;
                }
            }
            Q_L1.push_back(current_l);
            Q_R1.push_back(current_r);
        }

        std::priority_queue<PFI> MIDG_G(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<int> &Q_L, std::vector<int> &Q_R, int edge_limit, int id)
        {
            std::vector<int> Q_L1, Q_R1, mapping;
            merge_small_gaps(Q_L, Q_R, Q_L1, Q_R1, mapping, 2);
            int interval_num = Q_L1.size();

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFIMI, std::vector<PFIMI>, Compare> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                {
                    continue;
                }
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    top_candidates.emplace(dis, pid);
                    lowerBound = top_candidates.top().first;
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    top_candidates.emplace(dis, pid);
                    top_candidates.pop();
                    lowerBound = top_candidates.top().first;
                }
            }
            bool erfen = interval_num < 16;
            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.dis > lowerBound)
                    break;
                candidate_set.pop();
                int current_pid = current_point_pair.id;
                int range_id = current_point_pair.range_id;
                if (range_id == -1)
                {
                    if (erfen)
                    {
                        for (int i = 0; i < interval_num; i++)
                        {
                            if (current_pid >= Q_L1[i])
                            {
                                if (current_pid <= Q_R1[i])
                                {
                                    range_id = i;
                                    break;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        auto it = std::lower_bound(Q_L1.begin(), Q_L1.end(), current_pid);
                        if (it != Q_L1.end())
                        {
                            int i = static_cast<int>(std::distance(Q_L1.begin(), it));
                            if (*it == current_pid)
                            {
                                range_id = i;
                            }
                            else
                            {
                                range_id = i - 1;
                            }
                        }
                        else
                        {
                            range_id = interval_num - 1;
                        }
                    }
                }
                auto selected_edges = SelectEdge(current_pid, Q_L1[range_id], Q_R1[range_id], edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }

                if (interval_num == 1)
                {
                    continue;
                }
                std::vector<tableint> edges;
                edges.reserve(edge_limit / 2);
                // JINS_P(Q_L, Q_R, current_pid, edge_limit / 2, 2, edges);
                JINS_P(Q_L1, Q_R1, current_pid, edge_limit / 2, 2, edges);
                int num_edges1 = edges.size();
                for (int i = 0; i < std::min(num_edges1, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges1; ++i)
                {
                    int neighbor_id = edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;
                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }


        void JINS(std::vector<int> &QL, std::vector<int> &QR, int pid, int edge_limit, int deep, std::vector<std::pair<tableint, int>> &edges)
        {
            int Qsize = QL.size();
            if (deep == 0 || Qsize <= 1)
            {
                return;
            }
            int start_idx = 0;
            int end_idx = Qsize - 1;
            int left = QL[0];
            int right = QR[end_idx];
            TreeNode *cur_node = nullptr;
            TreeNode *nxt_node = tree->root;
            cur_node = nxt_node;
            bool contain = false;
            while (cur_node)
            {
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, left, right) == GetOverLap(nxt_node->lbound, nxt_node->rbound, left, right))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);
                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);
                for (int j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < left || neighborId > right)
                        continue;
                    // count--;
                    int inRange = -1;
                    for (int i = 0; i < Qsize; i++)
                    {
                        if (neighborId >= QL[i])
                        {
                            if (neighborId <= QR[i])
                            {
                                inRange = i;
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    if (inRange != -1)
                    {
                        edges.emplace_back(neighborId, inRange);
                    }
                    else
                    {
                        JINS(QL, QR, neighborId, edge_limit, deep - 1, edges);
                    }
                    if (edges.size() >= edge_limit)
                    {
                        return;
                    }
                }
                if (nxt_node == nullptr)
                {
                    return;
                }
                cur_node = nxt_node;
            }
        }

        void JINS_P(std::vector<int> &QL, std::vector<int> &QR, int pid, int edge_limit, int deep, std::vector<tableint> &edges)
        {
            int Qsize = QL.size();
            if (deep == 0 || Qsize <= 1)
            {
                return;
            }
            int start_idx = 0;
            int end_idx = Qsize - 1;
            int left = QL[0];
            int right = QR[end_idx];
            TreeNode *cur_node = nullptr;
            TreeNode *nxt_node = tree->root;
            cur_node = nxt_node;
            bool contain = false;
            while (cur_node)
            {
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, left, right) == GetOverLap(nxt_node->lbound, nxt_node->rbound, left, right))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);
                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);
                for (int j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < left || neighborId > right)
                        continue;
                    // count--;
                    if (flags[neighborId])
                    {
                        edges.emplace_back(neighborId);
                    }
                    else
                    {
                        JINS_P(QL, QR, neighborId, edge_limit, deep - 1, edges);
                    }
                    if (edges.size() >= edge_limit)
                    {
                        return;
                    }
                }
                if (nxt_node == nullptr)
                {
                    return;
                }
                cur_node = nxt_node;
                int low = start_idx;
                int high = end_idx;
                start_idx = end_idx + 1;

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1);
                    if (QR[mid] >= cur_node->lbound)
                    {
                        start_idx = mid;
                        high = mid - 1;
                    }
                    else
                    {
                        low = mid + 1;
                    }
                }

                if (start_idx > end_idx)
                {
                    return;
                }

                low = start_idx;
                high = end_idx;
                end_idx = start_idx;

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1);
                    if (QL[mid] > cur_node->rbound)
                    {
                        high = mid - 1;
                    }
                    else
                    {
                        end_idx = mid;
                        low = mid + 1;
                    }
                }

                left = QL[start_idx];
                right = QR[end_idx];
            }
        }

        std::priority_queue<PFI> Postfilter(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<int> &QL, std::vector<int> &QR, int edge_limit)
        {
            int interval_num = QL.size() - 1;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFI, std::vector<PFI>, std::greater<PFI>> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);
            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid);
                    if (flags[pid])
                    {
                        top_candidates.emplace(dis, pid);
                        lowerBound = top_candidates.top().first;
                    }
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid);
                    if (flags[pid])
                    {
                        top_candidates.emplace(dis, pid);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second;
                auto selected_edges = SelectEdge(current_pid, QL[0], QR[interval_num], edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];

                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id);
                        if (flags[neighbor_id])
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        std::priority_queue<PFI> TopDown_search(const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit, std::vector<iRangeGraph::TreeNode *> &filterednodes)
        {
            searcher::Bitset<uint64_t> visited_set(max_elements_);
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFII, std::vector<PFII>, std::greater<PFII>> candidate_set;
            std::priority_queue<PFI> top_candidates;

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                    continue;
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(std::make_pair(dis, std::make_pair(pid, -1)));
                if (flags[pid])
                    top_candidates.emplace(dis, pid);
            }

            float lowerBound = std::numeric_limits<float>::max();

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                metric_hops++;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second.first;
                int current_step = current_point_pair.second.second;
                auto selected_edges = SelectEdgeF(current_pid, QL, QR, edge_limit, current_step);

                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i].first), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i].first;
                    bool inrange = selected_edges[i].second;

                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    metric_distance_computations++;

                    if (top_candidates.size() < ef || dis < lowerBound)
                    {
                        int next_step = current_step + 1;
                        if (inrange)
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            next_step = -1;
                        }
                        candidate_set.emplace(std::make_pair(dis, std::make_pair(neighbor_id, next_step)));

                        if (top_candidates.size() > ef)
                            top_candidates.pop();
                        if (top_candidates.size())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }

            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        void search_Prefilter(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;
                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "_Prefilter.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        std::vector<int> ql, qr;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<std::priority_queue<PFI>> res_list;
                        res_list.reserve(interval_num);
                        for (int x = 0; x < interval_num; ++x)
                        {
                            std::vector<TreeNode *> filterednodes = tree->range_filter(tree->root, ql[x], qr[x]);
                            std::priority_queue<PFI> res = TopDown_nodeentries_search(filterednodes, storage->query_points[i].data(), ef, k, ql[x], qr[x], edge_limit);
                            res_list.push_back(res);
                        }
                        std::priority_queue<PFI> res;
                        float lowerBound = std::numeric_limits<float>::max();
                        for (int x = 0; x < interval_num; x++)
                        {
                            while (res_list[x].size())
                            {
                                if (res.size() >= k)
                                {
                                    if (res_list[x].top().first < lowerBound)
                                    {
                                        res.push(res_list[x].top());
                                        res.pop();
                                    }
                                    res_list[x].pop();
                                }
                                else
                                {
                                    res.push(res_list[x].top());
                                    res_list[x].pop();
                                }
                                lowerBound = res.top().first;
                            }
                        }
                        while (res.size() > k)
                        {
                            res.pop();
                        }
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                            {
                                continue;
                            }
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;
                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }

                outfile.close();
            }
        }

        void search_MIDG_P(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;

                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "MIDG-P.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {

                        std::fill(flags.begin(), flags.end(), false);
                        std::vector<int> ql, qr;
                        int size = 0;
                        int id = 0;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                            if (storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first > size)
                            {
                                size = storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first;
                                id = x;
                            };
                            int start = storage->query_range[x][suffix][i].first;
                            int end = storage->query_range[x][suffix][i].second;
                            std::fill(flags.begin() + start, flags.begin() + end + 1, true);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<TreeNode *> node = tree->range_filter(tree->root, ql[id], qr[id]);
                        std::priority_queue<PFI> res = MIDG_P(node, storage->query_points[i].data(), ef, k, ql, qr, edge_limit, id);
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;
                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }
                outfile.close();
            }
        }

        void search_MIDG(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;
                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "_MIDG.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {

                        std::fill(flags.begin(), flags.end(), false);
                        std::vector<int> ql, qr;
                        int size = 0;
                        int id = 0;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                            if (storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first > size)
                            {
                                size = storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first;
                                id = x;
                            };
                            int start = storage->query_range[x][suffix][i].first;
                            int end = storage->query_range[x][suffix][i].second;
                            std::fill(flags.begin() + start, flags.begin() + end + 1, true);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<TreeNode *> node = tree->range_filter(tree->root, ql[id], qr[id]);
                        std::priority_queue<PFI> res = MIDG(node, storage->query_points[i].data(), ef, k, ql, qr, edge_limit, id);
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;
                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }
                outfile.close();
            }
        }

        void
        searchor_MIDG_G(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;

                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "_MIDG_G.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {

                        std::fill(flags.begin(), flags.end(), false);
                        std::vector<int> ql, qr;
                        int size = 0;
                        int id = 0;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                            if (storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first > size)
                            {
                                size = storage->query_range[x][suffix][i].second - storage->query_range[x][suffix][i].first;
                                id = x;
                            };
                            int start = storage->query_range[x][suffix][i].first;
                            int end = storage->query_range[x][suffix][i].second;
                            std::fill(flags.begin() + start, flags.begin() + end + 1, true);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<TreeNode *> node = tree->range_filter(tree->root, ql[id], qr[id]);
                        std::priority_queue<PFI> res = MIDG_G(node, storage->query_points[i].data(), ef, k, ql, qr, edge_limit, id);
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;
                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }
                outfile.close();
            }
        }


        void search_Postfilter(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;
                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "_postfilter.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::vector<int> HOP;
                std::vector<int> DCO;
                std::vector<float> QPS;
                std::vector<float> RECALL;

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;
                    metric_hops1 = 0;
                    metric_distance_computations1 = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        std::fill(flags.begin(), flags.end(), false);

                        std::vector<int> ql, qr;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                            int start = storage->query_range[x][suffix][i].first;
                            int end = storage->query_range[x][suffix][i].second;
                            std::fill(flags.begin() + start, flags.begin() + end + 1, true);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<TreeNode *> filterednodes = tree->range_filter(tree->root, ql[0], qr[interval_num - 1]);
                        std::priority_queue<PFI> res = Postfilter(filterednodes, storage->query_points[i].data(), ef, k, ql, qr, edge_limit);
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    HOP.emplace_back(hop);
                    DCO.emplace_back(dco);
                    QPS.emplace_back(qps);
                    RECALL.emplace_back(recall);

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }

                outfile.close();
            }
        }

        void search_Postfilter_P(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int interval_num = storage->query_range.size();
            for (auto range : storage->query_range[0])
            {
                int suffix = range.first;
                std::vector<std::vector<int>> &gt = storage->groundtruth[suffix];
                std::string savepath = saveprefix + std::to_string(suffix) + "_postfilterP.csv";
                CheckPath(savepath);
                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);

                std::vector<int> HOP;
                std::vector<int> DCO;
                std::vector<float> QPS;
                std::vector<float> RECALL;

                std::cout << "suffix = " << suffix << std::endl;
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;
                    metric_hops1 = 0;
                    metric_distance_computations1 = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        std::fill(flags.begin(), flags.end(), false);
                        std::vector<int> ql, qr;
                        for (int x = 0; x < interval_num; x++)
                        {
                            ql.push_back(storage->query_range[x][suffix][i].first);
                            qr.push_back(storage->query_range[x][suffix][i].second);
                            int start = storage->query_range[x][suffix][i].first;
                            int end = storage->query_range[x][suffix][i].second;
                            std::fill(flags.begin() + start, flags.begin() + end + 1, true);
                        }
                        int k = storage->query_K;
                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        std::vector<TreeNode *> filterednodes = tree->range_filter(tree->root, ql[0], qr[interval_num - 1]);
                        std::priority_queue<PFI> res = TopDown_search(storage->query_points[i].data(), ef, k, ql[0], qr[interval_num - 1], edge_limit, filterednodes);
                        gettimeofday(&t2, NULL);
                        auto duration = GetTime(t1, t2);
                        searchtime += duration;
                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    HOP.emplace_back(hop);
                    DCO.emplace_back(dco);
                    QPS.emplace_back(qps);
                    RECALL.emplace_back(recall);

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                    if (recall < 0.9)
                        break;
                }

                outfile.close();
            }
        }
    };
}