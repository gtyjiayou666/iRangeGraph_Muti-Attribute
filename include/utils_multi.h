#pragma once

#include "utils.h"
#include "z_order.h"

typedef struct Attr_Constraint
{
    std::vector<std::pair<int, int>> attr_constraints;
} Attr_Constraint;
typedef std::pair<float, int> PFI;
namespace iRangeGraph_multi
{

    struct TwoRangeQuery
    {
        int l1, r1, l2, r2;
    };

    class DataLoader
    {
    public:
        int Dim, query_nb, query_K;
        std::vector<std::vector<float>> query_points;
        int data_nb;
        std::vector<std::vector<float>> data_points;
        std::vector<int> original_id;

        int attr_nb{0};
        std::vector<std::vector<int>> attributes;
        std::vector<int> keys;
        std::vector<int> vids;
        std::vector<std::vector<int>> maxmin;
        std::vector<int> zorder;
        std::vector<int> lexorder;

        hnswlib::L2Space *space;

        std::unordered_map<std::string, std::vector<Attr_Constraint>> query_range_or;
        std::unordered_map<std::string, std::vector<std::pair<int, int>>> query_range;
        std::unordered_map<std::string, std::vector<std::vector<int>>> ground_truth;

        std::unordered_map<std::string, std::vector<std::pair<int, int>>> mapped_queryrange;

        DataLoader() {}
        ~DataLoader() {}

        void CalZorder()
        {
            int d = attributes[0].size();
            if (d <= 1)
            {
                return;
            }
            int count = attributes.size();
            zorder.reserve(count);
            for (size_t i = 0; i < count; i++)
            {
                zorder.push_back(z_order_encode(attributes[i][0], attributes[i][1]));
            }
        }

        float dis_compute(std::vector<float> &v1, std::vector<float> &v2)
        {
            hnswlib::DISTFUNC<float> fstdistfunc_ = space->get_dist_func();
            float dis = fstdistfunc_((char *)v1.data(), (char *)v2.data(), space->get_dist_func_param());
            return dis;
        }

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
            space = new hnswlib::L2Space(Dim);
            infile.close();
        }

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
            attributes.resize(data_nb);
            infile.close();
        }
        void LoadAttribute(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
            {
                throw Exception("cannot open " + filename);
            }
            int val1;
            infile.read((char *)&val1, sizeof(int));
            for (int i = 0; i < data_nb; i++)
            {
                int val;
                infile.read((char *)&val, sizeof(int));
                attributes[i].emplace_back(val);
            }
            infile.close();
            attr_nb++;
        }
        void LoadKey(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
            {
                throw Exception("cannot open " + filename);
            }
            int val1;
            infile.read((char *)&val1, sizeof(int));
            std::cout<< data_nb << std::endl;
            for (int i = 0; i < data_nb; i++)
            {
                int val;
                infile.read((char *)&val, sizeof(int));
                keys.emplace_back(val);
            }
            infile.close();
        }

        void LoadVid(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
            {
                throw Exception("cannot open " + filename);
            }
            int val1;
            infile.read((char *)&val1, sizeof(int));
            for (int i = 0; i < data_nb; i++)
            {
                int val;
                infile.read((char *)&val, sizeof(int));
                vids.emplace_back(val);
            }
            infile.close();
        }

        bool check_amount(std::map<std::pair<std::string, std::string>, std::vector<TwoRangeQuery>> &mp)
        {
            for (auto t : mp)
            {
                if (t.second.size() < query_nb)
                    return false;
            }
            return true;
        }


        void LoadRangesOr(std::string file_attr0, std::string file_attr1, const std::string &query_key = "or")
        {
            std::ifstream infile0(file_attr0);
            std::ifstream infile1(file_attr1);

            if (!infile0.is_open() || !infile1.is_open())
            {
                throw std::runtime_error("Cannot open range files: " + file_attr0 + " or " + file_attr1);
            }

            std::string line0, line1;
            int count = 0;

            std::cout << "Loading OR ranges from:\n  - " << file_attr0 << "\n  - " << file_attr1 << std::endl;

            while (std::getline(infile0, line0) && std::getline(infile1, line1))
            {
                if (count == 0)
                {
                    count++;
                    continue;
                }

                std::vector<int> vals0 = parse_csv_line(line0); 
                if (vals0.size() < 2)
                    continue;
                int l0 = vals0[0];
                int r0 = vals0[1];

                std::vector<int> vals1 = parse_csv_line(line1);
                if (vals1.size() < 2)
                    continue;
                int l1 = vals1[0];
                int r1 = vals1[1];

                Attr_Constraint constraint;

                constraint.attr_constraints.emplace_back(l0, r0); 
                constraint.attr_constraints.emplace_back(l1, r1);

                query_range_or[query_key].push_back(constraint);

                count++;
            }

            infile0.close();
            infile1.close();

            std::cout << "Loaded " << (count - 1) << " OR constraints into map key: " << query_key << std::endl;
        }

        std::vector<int> parse_csv_line(std::string line)
        {
            std::stringstream ss(line);
            std::string cell;
            std::vector<int> values;

            while (std::getline(ss, cell, ','))
            {
                cell.erase(std::remove(cell.begin(), cell.end(), '\r'), cell.end());
                try
                {
                    values.push_back(std::stoi(cell));
                }
                catch (...)
                {
                }
            }
            return values;
        }
        void LoadRanges(std::string saveprefix)
        {
            std::string savepath = saveprefix;

            std::ifstream infile(savepath);
            if (!infile.is_open())
            {
                throw std::runtime_error("cannot open " + savepath);
            }

            std::string line;
            int count = 0;

            while (std::getline(infile, line))
            {
                if (count == 0)
                {
                    count++;
                    continue;
                }

                std::stringstream ss(line);
                std::string cell;
                std::vector<double> values;
                while (std::getline(ss, cell, ','))
                {
                    cell.erase(std::remove(cell.begin(), cell.end(), '\r'), cell.end());
                    values.push_back(std::stod(cell));
                }

                if (values.size() >= 2)
                {
                    double val_a = values[0];
                    double val_b = values[1];

                    int l1 = static_cast<int>(val_a);
                    int r1 = static_cast<int>(val_b);
                    query_range["or"].emplace_back(l1, r1);
                }
                count++;
            }
            infile.close();

            std::cout << "Loaded " << (count - 1) << " ranges from CSV." << std::endl;
        }

        void Generate_GroundtruthOr(std::string saveprefix)
        {
            for (const auto &entry : query_range_or)
            {
                std::string domain = entry.first;
                const std::vector<Attr_Constraint> &constraint_list = entry.second;

                std::string savepath = saveprefix + domain + ".bin";
                CheckPath(savepath);

                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                {
                    throw std::runtime_error("cannot open " + savepath); 
                }

                for (int i = 0; i < query_nb; i++)
                {
                    const Attr_Constraint &current_constraint = constraint_list[i];

                    std::priority_queue<std::pair<float, int>> ans;

                    for (int pid = 0; pid < data_nb; pid++)
                    {
                        bool is_match = false; 

                        for (const auto &range : current_constraint.attr_constraints)
                        {
                            int attr_idx = &range - &current_constraint.attr_constraints[0]; 
                            int val = attributes[pid][attr_idx];     
                            int ql = range.first; 
                            int qr = range.second;
                            if (val >= ql && val <= qr)
                            {
                                is_match = true;
                                break;
                            }
                        }

                        if (!is_match)
                            continue;

                        float dis = dis_compute(query_points[i], data_points[pid]);
                        ans.emplace(dis, vids[pid]);

                        if (ans.size() > query_K)
                        {
                            ans.pop();
                        }
                    }

                    if (ans.size() < query_K)
                    {
                        std::cout << "Warning: Query " << i << " only found " << ans.size() << " results." << std::endl;
                    }

                    while (ans.size())
                    {
                        auto id = ans.top().second;
                        ans.pop();
                        outfile.write((char *)&id, sizeof(int));
                    }

                }
                outfile.close();
                std::cout << "Finished domain: " << domain << std::endl;
            }
        }

        void Generate_GroundtruthOrzorder(std::string saveprefix)
        {
            for (const auto &entry : query_range_or)
            {
                std::string domain = entry.first;
                const std::vector<Attr_Constraint> &constraint_list = entry.second;

                std::string savepath = saveprefix + domain + ".bin";
                CheckPath(savepath);

                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                {
                    throw std::runtime_error("cannot open " + savepath); 
                }

                for (int i = 0; i < query_nb; i++)
                {
                    const Attr_Constraint &current_constraint = constraint_list[i];

                    std::priority_queue<std::pair<float, int>> ans;

                    for (int pid = 0; pid < data_nb; pid++)
                    {
                        bool is_match = false; 

                        for (const auto &range : current_constraint.attr_constraints)
                        {
                            int attr_idx = &range - &current_constraint.attr_constraints[0];
                            int val = attributes[pid][attr_idx];                   
                            int ql = range.first;
                            int qr = range.second;
                            if (val >= ql && val <= qr)
                            {
                                is_match = true;
                                break; 
                            }
                        }

                        if (!is_match)
                            continue;

                        float dis = dis_compute(query_points[i], data_points[pid]);
                        ans.emplace(dis, pid);

                        if (ans.size() > query_K)
                        {
                            ans.pop();
                        }
                    }

                    if (ans.size() < query_K)
                    {
                        std::cout << "Warning: Query " << i << " only found " << ans.size() << " results." << std::endl;
                    }

                    while (ans.size())
                    {
                        auto id = ans.top().second;
                        ans.pop();
                        outfile.write((char *)&id, sizeof(int));
                    }
                }
                outfile.close();
                std::cout << "Finished domain: " << domain << std::endl;
            }
        }
        void Generate_Groundtruth(std::string saveprefix)
        {
            for (auto t : query_range)
            {
                std::string domain = t.first;
                std::string savepath = saveprefix + domain + ".bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                for (int i = 0; i < query_nb; i++)
                {
                    auto constraints = t.second[i];
                    std::priority_queue<std::pair<float, int>> ans;
                    for (int pid = 0; pid < data_nb; pid++)
                    {
                        bool flag = true;
                        int ql = constraints.first, qr = constraints.second;
                        if (attributes[pid][1] < ql || attributes[pid][1] > qr)
                            flag = false;

                        if (!flag)
                            continue;
                        float dis = dis_compute(query_points[i], data_points[pid]);
                        ans.emplace(dis, pid);
                        if (ans.size() > query_K)
                            ans.pop();
                    }
                    if (ans.size() < query_K)
                    {
                        std::cout << "error :" << ans.size() << std::endl;
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


        void LoadGroundtruthOr(std::string saveprefix)
        {
            for (auto t : query_range_or)
            {
                std::string domain = t.first;
                std::string savepath = saveprefix + domain + ".bin";
                std::ifstream infile(savepath, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                ground_truth[domain].resize(query_nb);
                for (int qid = 0; qid < query_nb; qid++)
                {
                    ground_truth[domain][qid].resize(query_K);
                    infile.read((char *)ground_truth[domain][qid].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }
        void LoadGroundtruth(std::string saveprefix)
        {
            for (auto t : query_range)
            {
                std::string domain = t.first;
                std::string savepath = saveprefix + domain + ".bin";
                std::ifstream infile(savepath, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                ground_truth[domain].resize(query_nb);
                for (int qid = 0; qid < query_nb; qid++)
                {
                    ground_truth[domain][qid].resize(query_K);
                    infile.read((char *)ground_truth[domain][qid].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }
    };
}
