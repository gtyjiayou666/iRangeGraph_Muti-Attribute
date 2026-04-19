#pragma once

#include "space_l2.h"
#include <filesystem>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <map>
#include <thread>
class Exception : public std::runtime_error
{
public:
    Exception(const std::string &msg) : std::runtime_error(msg) {}
};

void CheckPath(std::string filename)
{
    std::filesystem::path pathObj(filename);
    std::filesystem::path dirPath = pathObj.parent_path();
    if (!std::filesystem::exists(dirPath))
    {
        try
        {
            if (std::filesystem::create_directories(dirPath))
            {
                std::cout << "Directory created: " << dirPath << std::endl;
            }
            else
            {
                std::cerr << "Failed to create directory: " << dirPath << std::endl;
            }
        }
        catch (std::filesystem::filesystem_error &e)
        {
            throw Exception(e.what());
        }
    }
}

float GetTime(timeval &begin, timeval &end)
{
    return end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
}

namespace iRangeGraph
{
    typedef std::pair<float, int> PFI;
    typedef unsigned int tableint;
    typedef unsigned int linklistsizeint;




    struct PFIMI
    {
        float dis;
        int id;
        int range_id;
        PFIMI(float d, int i, int range_id)
            : dis(d), id(i), range_id(range_id) {}
    };

    
    struct Compare
    {
        bool operator()(const PFIMI &a, const PFIMI &b)
        {
            return a.dis > b.dis;
        }
    };


    class DataLoaderMultiInterval
    {
    public:
        int Dim, query_nb, query_K;
        std::vector<std::vector<float>> query_points;
        int data_nb;
        std::vector<std::vector<float>> data_points;
        std::vector<std::unordered_map<int, std::vector<std::pair<int, int>>>> query_range;
        std::unordered_map<int, std::vector<std::vector<int>>> groundtruth;
        int interval_num;

        DataLoaderMultiInterval(int interval_num) : interval_num(interval_num)
        {
            query_range.resize(interval_num);
        }
        ~DataLoaderMultiInterval() {}

        // query vector filename format: 4 bytes: query number; 4 bytes: dimension; query_nb*Dim vectors
        void LoadQuery(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&query_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            std::cout << "dom=" << Dim << std::endl;
            query_points.resize(query_nb);
            for (int i = 0; i < query_nb; i++)
            {
                query_points[i].resize(Dim);
                infile.read((char *)query_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        // Used only when computing groundtruth and constructing index. Do not use this to load data for search process
        void LoadData(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&data_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            data_points.resize(data_nb);
            for (int i = 0; i < data_nb; i++)
            {
                data_points[i].resize(Dim);
                infile.read((char *)data_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        void LoadQueryRange(std::string fileprefix)
        {
            std::vector<int> s;
            for (int i = 0; i < 15; i++)
                s.emplace_back(i);
            for (auto suffix : s)
            {
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                int interval_num;
                infile.read(reinterpret_cast<char *>(&interval_num), sizeof(interval_num));
                for (int i = 0; i < query_nb; i++)
                {

                    for (int j = 0; j < interval_num; ++j)
                    {
                        int ql, qr;
                        infile.read(reinterpret_cast<char *>(&ql), sizeof(ql));
                        infile.read(reinterpret_cast<char *>(&qr), sizeof(qr));

                        if (infile.fail())
                        {
                            throw std::runtime_error("File " + filename + " corrupted at query " + std::to_string(i));
                        }
                        query_range[j][suffix].emplace_back(ql, qr);
                    }
                }
                infile.close();
            }
        }

        void LoadQueryRangeLargeGaps(std::string fileprefix)
        {
            std::vector<int> s;
            for (int i = 0; i < 15; i++)
                s.emplace_back(i);
            for (auto suffix : s)
            {
                std::string filename = fileprefix + std::to_string(suffix) + "largegaps.bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                int interval_num;
                infile.read(reinterpret_cast<char *>(&interval_num), sizeof(interval_num));
                for (int i = 0; i < query_nb; i++)
                {

                    for (int j = 0; j < interval_num; ++j)
                    {
                        int ql, qr;
                        infile.read(reinterpret_cast<char *>(&ql), sizeof(ql));
                        infile.read(reinterpret_cast<char *>(&qr), sizeof(qr));
                        if (infile.fail())
                        {
                            throw std::runtime_error("File " + filename + " corrupted at query " + std::to_string(i));
                        }
                        query_range[j][suffix].emplace_back(ql, qr);
                    }
                }
                infile.close();
            }
        }
        
        void LoadGroundtruth(std::string fileprefix)
        {
            for (auto t : query_range[0])
            {
                int suffix = t.first;
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                groundtruth[suffix].resize(query_nb);
                for (int i = 0; i < query_nb; i++)
                {
                    groundtruth[suffix][i].resize(query_K);
                    infile.read((char *)groundtruth[suffix][i].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }

        void LoadGroundtruthpro(std::string fileprefix)
        {
            for (auto t : query_range[0])
            {
                int suffix = t.first;
                std::string filename = fileprefix + std::to_string(suffix) + "pro.bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                groundtruth[suffix].resize(query_nb);
                for (int i = 0; i < query_nb; i++)
                {
                    groundtruth[suffix][i].resize(query_K);
                    infile.read((char *)groundtruth[suffix][i].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }
    };

    class QueryGeneratorMultiInterval
    {
    public:
        int data_nb, query_nb, interval_num;
        hnswlib::L2Space *space;

        QueryGeneratorMultiInterval(int data_num, int query_num, int interval_num) : data_nb(data_num), query_nb(query_num), interval_num(interval_num) {}
        ~QueryGeneratorMultiInterval() {}

        std::vector<int> split_positive(int total, int count, std::default_random_engine &e)
        {
            if (total < count)
            {
                throw Exception("Cannot split " + std::to_string(total) + " into " + std::to_string(count) + " positive integers.");
            }
            std::vector<int> cuts(count - 1);
            std::uniform_int_distribution<int> dist(1, total - 1);
            for (int i = 0; i < count - 1; ++i)
            {
                cuts[i] = dist(e);
            }
            std::sort(cuts.begin(), cuts.end());
            cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

            while (static_cast<int>(cuts.size()) < count - 1)
            {
                int x = dist(e);
                if (std::find(cuts.begin(), cuts.end(), x) == cuts.end())
                {
                    cuts.push_back(x);
                    std::sort(cuts.begin(), cuts.end());
                }
            }

            std::vector<int> parts;
            int last = 0;
            for (int cut : cuts)
            {
                parts.push_back(cut - last);
                last = cut;
            }
            parts.push_back(total - last);
            return parts;
        }

        std::vector<int> split_nonnegative(int total, int count, std::default_random_engine &e)
        {
            if (total < 0)
            {
                throw Exception("Total gap cannot be negative!");
            }
            if (count == 1)
            {
                return {total};
            }
            std::vector<int> cuts(count - 1);
            std::uniform_int_distribution<int> dist(0, total + count - 2);
            for (int i = 0; i < count - 1; ++i)
            {
                cuts[i] = dist(e);
            }
            std::sort(cuts.begin(), cuts.end());
            cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

            while (static_cast<int>(cuts.size()) < count - 1)
            {
                int x = dist(e);
                if (std::find(cuts.begin(), cuts.end(), x) == cuts.end())
                {
                    cuts.push_back(x);
                    std::sort(cuts.begin(), cuts.end());
                }
            }

            std::vector<int> parts;
            int last = -1;
            for (int cut : cuts)
            {
                parts.push_back(cut - last - 1);
                last = cut;
            }
            parts.push_back(total + count - 1 - last - 1);
            return parts;
        }
        std::vector<int> split_fixed_big_gaps(int total, int count, std::default_random_engine &e,
                                              int min_big = 3, int max_big = 10)
        {
            if (total < 0)
            {
                throw Exception("Total gap cannot be negative!");
            }
            if (count <= 0)
            {
                return {};
            }
            if (count == 1)
            {
                return {total};
            }
            if (total == 0)
            {
                return std::vector<int>(count, 0);
            }

            int actual_max_big = std::min(max_big, count);
            int actual_min_big = std::min(min_big, actual_max_big);

            std::uniform_int_distribution<int> dist_k(actual_min_big, actual_max_big);
            int K = dist_k(e);

            std::uniform_real_distribution<double> dist_big(100.0, 1000.0);
            std::uniform_real_distribution<double> dist_small(0.0, 1.0);

            std::vector<double> weights(count, 0.0);
            std::vector<int> indices(count);
            std::iota(indices.begin(), indices.end(), 0);

            std::shuffle(indices.begin(), indices.end(), e);

            double weight_sum = 0.0;
            for (int i = 0; i < count; ++i)
            {
                double w;
                if (i < K)
                {
                    w = dist_big(e);
                }
                else
                {
                    w = dist_small(e);
                }
                weights[indices[i]] = w;
                weight_sum += w;
            }

            std::vector<int> parts(count);
            int allocated_sum = 0;

            for (int i = 0; i < count; ++i)
            {
                double ratio = weights[i] / weight_sum;
                int val = static_cast<int>(total * ratio);
                parts[i] = val;
                allocated_sum += val;
            }

            int remainder = total - allocated_sum;
            if (remainder > 0)
            {
                std::vector<int> big_indices;
                big_indices.reserve(K);
                for (int i = 0; i < count; ++i)
                {
                    if (i < K)
                        big_indices.push_back(indices[i]);
                }

                if (big_indices.empty())
                {
                    std::shuffle(indices.begin(), indices.end(), e);
                    for (int i = 0; i < remainder; ++i)
                        parts[indices[i % count]]++;
                }
                else
                {
                    std::shuffle(big_indices.begin(), big_indices.end(), e);
                    for (int i = 0; i < remainder; ++i)
                    {
                        parts[big_indices[i % K]]++;
                    }
                }
            }

            return parts;
        }

        void GenerateMultiInterval(std::string saveprefix, int I_min_num)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::vector<float> s =
                {0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01};
            std::vector<std::pair<int, int>> rs;
            for (int i = 0; i < 15; i++)
            {
                rs.emplace_back(I_min_num * s[i], i);
            }
            for (auto t : rs)
            {
                int total_len = t.first;
                int suffix = t.second;
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::cout << "save query range to " << savepath << std::endl;
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                outfile.write((char *)(&interval_num), sizeof(int));

                if (data_nb < total_len + 1)
                {
                    throw Exception("data_nb too small to fit two non-overlapping ranges with total length " + std::to_string(total_len));
                }

                for (int i = 0; i < query_nb; i++)
                {
                    std::vector<int> lengths = split_positive(total_len, interval_num, e);

                    int gap_total = (I_min_num - total_len);
                    std::vector<int> gaps = split_nonnegative(gap_total, interval_num - 1, e);
                    gaps.push_back(0);
                    std::vector<std::pair<int, int>> ranges;
                    std::uniform_int_distribution<int> dist(0, data_nb - I_min_num);

                    int pos = dist(e);
                    for (int j = 0; j < interval_num; ++j)
                    {
                        int start = pos;
                        int end = start + lengths[j] - 1;
                        ranges.emplace_back(start, end);
                        pos = end + 1 + gaps[j];
                    }

                    for (const auto &r : ranges)
                    {
                        outfile.write(reinterpret_cast<const char *>(&r.first), sizeof(int));
                        outfile.write(reinterpret_cast<const char *>(&r.second), sizeof(int));
                    }
                }

                outfile.close();
            }
        }

        void GenerateMultiIntervalLargeGaps(std::string saveprefix, int I_min_num)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::vector<float> s =
                {0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01};
            std::vector<std::pair<int, int>> rs;
            for (int i = 0; i < 15; i++)
            {
                rs.emplace_back(I_min_num * s[i], i);
            }
            for (auto t : rs)
            {
                int total_len = t.first; 
                int suffix = t.second;
                std::string savepath = saveprefix + std::to_string(suffix) + "largegaps.bin";
                CheckPath(savepath);
                std::cout << "save query range to " << savepath << std::endl;
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                outfile.write((char *)(&interval_num), sizeof(int));

                if (data_nb < total_len + 1)
                {
                    throw Exception("data_nb too small to fit two non-overlapping ranges with total length " + std::to_string(total_len));
                }

                for (int i = 0; i < query_nb; i++)
                {
                    std::vector<int> lengths = split_positive(total_len, interval_num, e);
                    int gap_total = (I_min_num - total_len);
                    std::vector<int> gaps = split_fixed_big_gaps(gap_total, interval_num - 1, e);

                    gaps.push_back(0);
                    std::vector<std::pair<int, int>> ranges;
                    std::uniform_int_distribution<int> dist(0, data_nb - I_min_num);

                    int pos = dist(e);
                    for (int j = 0; j < interval_num; ++j)
                    {
                        int start = pos;
                        int end = start + lengths[j] - 1;
                        ranges.emplace_back(start, end);
                        pos = end + 1 + gaps[j]; 
                    }

                    for (const auto &r : ranges)
                    {
                        outfile.write(reinterpret_cast<const char *>(&r.first), sizeof(int));
                        outfile.write(reinterpret_cast<const char *>(&r.second), sizeof(int));
                    }
                }

                outfile.close();
            }
        }
        float dis_compute(std::vector<float> &v1, std::vector<float> &v2)
        {
            hnswlib::DISTFUNC<float> fstdistfunc_ = space->get_dist_func();
            float dis = fstdistfunc_((char *)v1.data(), (char *)v2.data(), space->get_dist_func_param());
            return dis;
        }

        void GenerateGroundtruthThread(std::string saveprefix, DataLoaderMultiInterval &storage)
        {
            space = new hnswlib::L2Space(storage.Dim);
            const int num_suffixes = 15;
            std::vector<std::thread> threads;
            threads.reserve(num_suffixes);

            for (int suffix = 0; suffix < num_suffixes; ++suffix)
            {
                threads.emplace_back([&, suffix]()
                                     {
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open()) {
                    std::cerr << "Error: cannot open " << savepath << std::endl;
                    return;
                }

                std::cout << "Generating ground truth for suffix " << suffix << std::endl;

                for (int i = 0; i < storage.query_nb; ++i)
                {
                    std::priority_queue<std::pair<float, int>> ans; // max-heap

                    for (int k = 0; k < interval_num; ++k)
                    {
                        auto [ql, qr] = storage.query_range[k].at(suffix)[i];
                        for (int j = ql; j <= qr; ++j)
                        {
                            float dis = dis_compute(storage.query_points[i], storage.data_points[j]);
                            ans.emplace(dis, j);
                            if ((int)ans.size() > storage.query_K)
                            {
                                ans.pop();
                            }
                        }
                    }

                    std::vector<int> results;
                    while (!ans.empty())
                    {
                        results.push_back(ans.top().second);
                        ans.pop();
                    }
                    std::reverse(results.begin(), results.end());

                    for (int id : results)
                    {
                        outfile.write(reinterpret_cast<const char*>(&id), sizeof(int));
                    }
                }

                outfile.close();
                std::cout << "Finished suffix " << suffix << std::endl; });
            }

            for (auto &t : threads)
            {
                if (t.joinable())
                    t.join();
            }
        }

        void GenerateGroundtruthThreadLargeGaps(std::string saveprefix, DataLoaderMultiInterval &storage)
        {
            space = new hnswlib::L2Space(storage.Dim);
            const int num_suffixes = 15;
            std::vector<std::thread> threads;
            threads.reserve(num_suffixes);

            for (int suffix = 0; suffix < num_suffixes; ++suffix)
            {
                threads.emplace_back([&, suffix]()
                                     {
                std::string savepath = saveprefix + std::to_string(suffix) + "largegaps.bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open()) {
                    std::cerr << "Error: cannot open " << savepath << std::endl;
                    return;
                }

                std::cout << "Generating ground truth for suffix " << suffix << std::endl;

                for (int i = 0; i < storage.query_nb; ++i)
                {
                    std::priority_queue<std::pair<float, int>> ans; // max-heap

                    for (int k = 0; k < interval_num; ++k)
                    {
                        auto [ql, qr] = storage.query_range[k].at(suffix)[i];
                        for (int j = ql; j <= qr; ++j)
                        {
                            float dis = dis_compute(storage.query_points[i], storage.data_points[j]);
                            ans.emplace(dis, j);
                            if ((int)ans.size() > storage.query_K)
                            {
                                ans.pop();
                            }
                        }
                    }

                    std::vector<int> results;
                    while (!ans.empty())
                    {
                        results.push_back(ans.top().second);
                        ans.pop();
                    }
                    std::reverse(results.begin(), results.end());
                    for (int id : results)
                    {
                        outfile.write(reinterpret_cast<const char*>(&id), sizeof(int));
                    }
                }

                outfile.close();
                std::cout << "Finished suffix " << suffix << std::endl; });
            }

            for (auto &t : threads)
            {
                if (t.joinable())
                    t.join();
            }
        }
    };


















    class DataLoader
    {
    public:
        int Dim, query_nb, query_K;
        std::vector<std::vector<float>> query_points;
        int data_nb;
        std::vector<std::vector<float>> data_points;
        std::unordered_map<int, std::vector<std::pair<int, int>>> query_range;
        std::unordered_map<int, std::vector<std::vector<int>>> groundtruth;

        DataLoader() {}
        ~DataLoader() {}

        // query vector filename format: 4 bytes: query number; 4 bytes: dimension; query_nb*Dim vectors
        void LoadQuery(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&query_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            query_points.resize(query_nb);
            for (int i = 0; i < query_nb; i++)
            {
                query_points[i].resize(Dim);
                infile.read((char *)query_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        // Used only when computing groundtruth and constructing index. Do not use this to load data for search process
        void LoadData(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&data_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            data_points.resize(data_nb);
            for (int i = 0; i < data_nb; i++)
            {
                data_points[i].resize(Dim);
                infile.read((char *)data_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        // By default generation, 0.bin~9.bin denotes 2^0~2^-9 range fractions, 17.bin denotes mixed range fraction.
        // Before reading the query ranges, make sure query vectors have been read.
        void LoadQueryRange(std::string fileprefix)
        {
            std::vector<int> s;
            for (int i = 0; i < 10; i++)
                s.emplace_back(i);
            s.emplace_back(17);
            for (auto suffix : s)
            {
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                for (int i = 0; i < query_nb; i++)
                {
                    int ql, qr;
                    infile.read((char *)&ql, sizeof(int));
                    infile.read((char *)&qr, sizeof(int));
                    query_range[suffix].emplace_back(ql, qr);
                }
                infile.close();
            }
        }

        // 0.bin~9.bin correspond to groundtruth for 2^0~2^-9 range fractions, 17.bin for mixed fraction
        void LoadGroundtruth(std::string fileprefix)
        {
            for (auto t : query_range)
            {
                int suffix = t.first;
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                groundtruth[suffix].resize(query_nb);
                for (int i = 0; i < query_nb; i++)
                {
                    groundtruth[suffix][i].resize(query_K);
                    infile.read((char *)groundtruth[suffix][i].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }
    };

    class QueryGenerator
    {
    public:
        int data_nb, query_nb;
        hnswlib::L2Space *space;

        QueryGenerator(int data_num, int query_num) : data_nb(data_num), query_nb(query_num) {}
        ~QueryGenerator() {}

        void GenerateRange(std::string saveprefix)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::vector<std::pair<int, int>> rs;
            int current_len = data_nb;
            for (int i = 0; i < 10; i++)
            {
                if (current_len < 10)
                    throw Exception("dataset size is too small, increase the amount of data objects!");
                rs.emplace_back(current_len, i);
                current_len /= 2;
            }
            for (auto t : rs)
            {
                int len = t.first, suffix = t.second;
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::cout << "save query range to" << savepath << std::endl;
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                std::uniform_int_distribution<int> u_start(0, data_nb - len);
                for (int i = 0; i < query_nb; i++)
                {
                    int ql = u_start(e);
                    int qr = ql + len - 1;
                    if (ql >= data_nb || qr >= data_nb)
                        throw Exception("Query range out of bound");
                    outfile.write((char *)&ql, sizeof(int));
                    outfile.write((char *)&qr, sizeof(int));
                }
                outfile.close();
            }

            rs.clear();
            current_len = data_nb;
            for (int i = 0; i < 10; i++)
            {
                rs.emplace_back(current_len, i);
                current_len /= 2;
            }
            std::string savepath = saveprefix + "17.bin";
            CheckPath(savepath);
            std::cout << "save query range to" << savepath << std::endl;
            std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
            if (!outfile.is_open())
                throw Exception("cannot open " + savepath);

            for (auto t : rs)
            {
                int len = t.first;
                std::uniform_int_distribution<int> u_start(0, data_nb - len);

                for (int i = 0; i < query_nb / 10; i++)
                {
                    int ql = u_start(e);
                    int qr = ql + len - 1;
                    if (ql >= data_nb || qr >= data_nb)
                        throw Exception("Query range out of bound");
                    outfile.write((char *)&ql, sizeof(int));
                    outfile.write((char *)&qr, sizeof(int));
                }
            }
            outfile.close();
        }

        float dis_compute(std::vector<float> &v1, std::vector<float> &v2)
        {
            hnswlib::DISTFUNC<float> fstdistfunc_ = space->get_dist_func();
            float dis = fstdistfunc_((char *)v1.data(), (char *)v2.data(), space->get_dist_func_param());
            return dis;
        }

        void GenerateGroundtruth(std::string saveprefix, DataLoader &storage)
        {
            space = new hnswlib::L2Space(storage.Dim);
            for (auto t : storage.query_range)
            {
                int suffix = t.first;
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                std::cout << "generating for " << t.first << std::endl;
                for (int i = 0; i < query_nb; i++)
                {
                    auto rp = t.second[i];
                    int ql = rp.first, qr = rp.second;
                    std::priority_queue<std::pair<float, int>> ans;
                    for (int j = ql; j <= qr; j++)
                    {
                        float dis = dis_compute(storage.query_points[i], storage.data_points[j]);
                        ans.emplace(dis, j);
                        if (ans.size() > storage.query_K)
                            ans.pop();
                    }
                    while (ans.size())
                    {
                        auto id = ans.top().second;
                        ans.pop();
                        outfile.write((char *)&id, sizeof(int));
                    }
                }
                outfile.close();
            }
        }
    };

    class TreeNode
    {
    public:
        int node_id;
        int lbound, rbound;
        int depth;
        std::vector<TreeNode *> childs;
        TreeNode(int l, int r, int d) : lbound(l), rbound(r), depth(d) {}
    };

    class SegmentTree
    {
    public:
        int ways_ = 2;
        TreeNode *root{nullptr};
        int max_depth{-1};
        std::vector<TreeNode *> treenodes;

        SegmentTree(int data_nb)
        {
            root = new TreeNode(0, data_nb - 1, 0);
        }

        void BuildTree(TreeNode *u)
        {
            if (u == nullptr)
                throw Exception("Tree node is a nullptr");
            treenodes.emplace_back(u);
            max_depth = std::max(max_depth, u->depth);
            int L = u->lbound, R = u->rbound;
            size_t Len = R - L + 1;
            if (L == R)
                return;
            int gap = (R - L + 1) / ways_;
            int res = (R - L + 1) % ways_;

            for (int l = L; l <= R;)
            {
                int r = l + gap - 1;
                if (res > 0)
                {
                    r++;
                    res--;
                }
                r = std::min(r, R);
                TreeNode *childnode = new TreeNode(l, r, u->depth + 1);
                u->childs.emplace_back(childnode);
                BuildTree(childnode);
                l = r + 1;
            }
        }

        std::vector<TreeNode *> range_filter(TreeNode *u, int ql, int qr)
        {
            if (u->lbound >= ql && u->rbound <= qr)
                return {u};
            std::vector<TreeNode *> res;
            if (u->lbound > qr)
                return res;
            if (u->rbound < ql)
                return res;
            for (auto child : u->childs)
            {
                auto t = range_filter(child, ql, qr);
                while (t.size())
                {
                    res.emplace_back(t.back());
                    t.pop_back();
                }
            }
            return res;
        }
    };
}