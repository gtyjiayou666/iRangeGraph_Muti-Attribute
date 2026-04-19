#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <climits>
#include <sstream>   
#include <z_order.h> 

std::vector<int> read_attr_csv(
    const std::string &filename,
    uint32_t n,
    int col_index = 0,    
    bool has_header = true,
    char delimiter = ',' 
)
{
    std::ifstream file(filename);
    if (!file)
    {
        throw std::runtime_error("Cannot open CSV: " + filename);
    }

    std::vector<int> data;
    std::string line;

    if (has_header)
    {
        if (!std::getline(file, line))
        {
            throw std::runtime_error("CSV file is empty");
        }
    }

    while (std::getline(file, line))
    {
        if (!line.empty() && line[line.length() - 1] == '\r')
        {
            line.erase(line.length() - 1);
        }

        if (line.empty())
            continue; 

        std::istringstream iss(line);
        std::string token;
        int current_col = 0;
        bool found = false;

        while (std::getline(iss, token, delimiter))
        {
            if (current_col == col_index)
            {
                try
                {
                    data.push_back(std::stoi(token));
                    found = true;
                }
                catch (const std::invalid_argument &e)
                {
                    throw std::runtime_error("Invalid integer format in column " + std::to_string(col_index) + ": " + token);
                }
                break;
            }
            current_col++;
        }

        if (!found)
        {
            throw std::runtime_error("Column index " + std::to_string(col_index) + " out of range in line: " + line);
        }
    }

    n = static_cast<uint32_t>(data.size());
    file.close();
    return data;
}

std::vector<std::vector<float>> read_vectors_bin(const std::string &filename, uint32_t n, uint32_t &dim)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        throw std::runtime_error("Cannot open: " + filename);
    }
    int x;
    file.read(reinterpret_cast<char *>(&x), sizeof(uint32_t));
    file.read(reinterpret_cast<char *>(&dim), sizeof(uint32_t));

    std::vector<std::vector<float>> vectors(n, std::vector<float>(dim));
    for (uint32_t i = 0; i < n; ++i)
    {
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
    std::cout << "[INFO] Saved: " << filename << "\n";
}

void write_vectors_bin(const std::string &filename, const std::vector<std::vector<float>> &vectors)
{
    std::ofstream file(filename, std::ios::binary);
    uint32_t n = static_cast<uint32_t>(vectors.size());
    uint32_t dim = static_cast<uint32_t>(vectors[0].size());
    file.write(reinterpret_cast<const char *>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char *>(&dim), sizeof(uint32_t));
    for (const auto &vec : vectors)
    {
        file.write(reinterpret_cast<const char *>(vec.data()), dim * sizeof(float));
    }
    file.close();
    std::cout << "[INFO] Saved: " << filename << "\n";
}



int main()
{
    const int d_attr = 2;
    uint32_t n = 1000000, dim = 0;

    std::vector<std::vector<int>> attributes(d_attr);
    std::vector<std::string> attr_files = {
        "/data/wit/wit_embeddings_metadata.csv",
        "/data/wit/wit_embeddings_metadata.csv"};
    for (int i = 0; i < d_attr; ++i)
    {
        uint32_t local_n;

        attributes[i] = read_attr_csv(attr_files[i], n, i + 3, true);
        std::cout << attributes[i][0] << std::endl;
        attributes[i].resize(1000000);
    }

    std::cout << "[INFO] Loaded " << n << " attributes from CSV.\n";

    auto vectors = read_vectors_bin("/data/wit/wit_embeddings.bin", n, dim);

    std::vector<std::pair<uint64_t, uint32_t>> z_indices;
    z_indices.reserve(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t x = (attributes[0][i] < 0) ? 0 : static_cast<uint32_t>(attributes[0][i]);
        uint32_t y = (attributes[1][i] < 0) ? 0 : static_cast<uint32_t>(attributes[1][i]);

        uint64_t z = z_order_encode(x, y);
        z_indices.emplace_back(z, i);
    }

    std::sort(z_indices.begin(), z_indices.end());

    std::vector<int> sorted_attr0(n), sorted_attr1(n);
    std::vector<std::vector<float>> sorted_vectors(n);

    for (uint32_t new_idx = 0; new_idx < n; ++new_idx)
    {
        uint32_t old_idx = z_indices[new_idx].second;
        sorted_attr0[new_idx] = attributes[0][old_idx];
        sorted_attr1[new_idx] = attributes[1][old_idx];
        sorted_vectors[new_idx] = vectors[old_idx];
    }

    write_attr_bin("/data/wit/wit.zorder.attr_0.bin", sorted_attr0);
    write_attr_bin("/data/wit/wit.zorder.attr_1.bin", sorted_attr1);
    write_vectors_bin("/data/wit/wit_zorder.bin", sorted_vectors);
    std::cout << "[SUCCESS] Data sorted by Z-order and saved.\n";
    return 0;
}