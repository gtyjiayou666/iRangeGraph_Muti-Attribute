#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <climits>
#include <sstream>
#include <tuple>

std::vector<int> read_attr_csv(
    const std::string &filename,
    int col_index = 0,
    bool has_header = true,
    char delimiter = ','
)
{
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open CSV: " + filename);

    std::vector<int> data;
    std::string line;

    if (has_header) {
        if (!std::getline(file, line)) throw std::runtime_error("CSV file is empty");
    }

    while (std::getline(file, line)) {
        if (!line.empty() && line[line.length() - 1] == '\r') line.erase(line.length() - 1);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        int current_col = 0;
        bool found = false;

        while (std::getline(iss, token, delimiter)) {
            if (current_col == col_index) {
                try {
                    data.push_back(std::stoi(token));
                    found = true;
                } catch (...) {
                    throw std::runtime_error("Invalid integer format");
                }
                break;
            }
            current_col++;
        }
        if (!found) throw std::runtime_error("Column index out of range");
    }
    return data;
}

std::vector<std::vector<float>> read_vectors_bin(const std::string &filename, uint32_t n, uint32_t &dim)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open: " + filename);
    
    uint32_t file_n;
    file.read(reinterpret_cast<char *>(&file_n), sizeof(uint32_t));
    file.read(reinterpret_cast<char *>(&dim), sizeof(uint32_t));

    uint32_t count = std::min(n, file_n);
    
    std::vector<std::vector<float>> vectors(count, std::vector<float>(dim));
    for (uint32_t i = 0; i < count; ++i) {
        file.read(reinterpret_cast<char *>(vectors[i].data()), dim * sizeof(float));
    }
    file.close();
    return vectors;
}

void write_attr_bin(const std::string &filename, const std::vector<int> &data)
{
    std::ofstream file(filename, std::ios::binary);
    uint32_t n = static_cast<uint32_t>(data.size());
    file.write(reinterpret_cast<const char *>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char *>(data.data()), n * sizeof(int));
    file.close();
}

void write_vectors_bin(const std::string &filename, const std::vector<std::vector<float>> &vectors)
{
    std::ofstream file(filename, std::ios::binary);
    uint32_t n = static_cast<uint32_t>(vectors.size());
    uint32_t dim = static_cast<uint32_t>(vectors[0].size());
    file.write(reinterpret_cast<const char *>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char *>(&dim), sizeof(uint32_t));
    for (const auto &vec : vectors) {
        file.write(reinterpret_cast<const char *>(vec.data()), dim * sizeof(float));
    }
    file.close();
}

void write_id_map_bin(const std::string &filename, const std::vector<uint32_t> &ids)
{
    std::ofstream file(filename, std::ios::binary);
    uint32_t n = static_cast<uint32_t>(ids.size());
    file.write(reinterpret_cast<const char *>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char *>(ids.data()), n * sizeof(uint32_t));
    file.close();
    std::cout << "   [Saved ID Map] " << filename << "\n";
}

void sort_and_save(
    const std::vector<int> &attr0,
    const std::vector<int> &attr1,
    const std::vector<std::vector<float>> &vectors,
    const std::string &output_prefix,
    const std::string &suffix
)
{
    uint32_t n = attr0.size();
    std::cout << "\n[PROCESS] Sorting by " << (suffix == "_x" ? "Attr0 -> Attr1" : "Attr1 -> Attr0") << "...\n";

    std::vector<std::tuple<int, int, uint32_t>> sort_keys;
    sort_keys.reserve(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        if (suffix == "_x") {
            sort_keys.emplace_back(attr0[i], attr1[i], i);
        } else {
            sort_keys.emplace_back(attr1[i], attr0[i], i);
        }
    }

    std::sort(sort_keys.begin(), sort_keys.end());

    std::vector<int> sorted_attr0(n);
    std::vector<int> sorted_attr1(n);
    std::vector<std::vector<float>> sorted_vectors(n);
    std::vector<uint32_t> original_ids(n);

    for (uint32_t new_idx = 0; new_idx < n; ++new_idx)
    {
        uint32_t old_idx = std::get<2>(sort_keys[new_idx]);
        
        sorted_attr0[new_idx] = attr0[old_idx];
        sorted_attr1[new_idx] = attr1[old_idx];
        sorted_vectors[new_idx] = vectors[old_idx];
        original_ids[new_idx] = old_idx;
    }

    write_attr_bin(output_prefix + suffix + ".attr_0.bin", sorted_attr0);
    write_attr_bin(output_prefix + suffix + ".attr_1.bin", sorted_attr1);
    write_vectors_bin(output_prefix + suffix + ".bin", sorted_vectors);
    write_id_map_bin(output_prefix + suffix + ".id_map.bin", original_ids);
}

int main()
{
    std::string csv_path = "/data/wit/wit_embeddings_metadata.csv";
    std::string bin_path = "/data/wit/wit_embeddings.bin";
    std::string output_prefix = "/data/wit/wit.lexorder";

    uint32_t dim = 0;

    try {
        std::cout << "[INFO] Reading CSV attributes...\n";
        std::vector<int> attributes_0 = read_attr_csv(csv_path, 0, true);
        std::vector<int> attributes_1 = read_attr_csv(csv_path, 1, true);

        uint32_t n = 1000000;
        attributes_0.resize(n);
        attributes_1.resize(n);

        std::cout << "[INFO] Reading Vectors...\n";
        auto vectors = read_vectors_bin(bin_path, n, dim);
        
        sort_and_save(attributes_0, attributes_1, vectors, output_prefix, "_x");

        sort_and_save(attributes_0, attributes_1, vectors, output_prefix, "_y");

        std::cout << "\n[SUCCESS] All tasks completed.\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}