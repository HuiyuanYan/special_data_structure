#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "src/KDTree.hpp"
#include <assert.h>
using namespace std;
#define GENERATE_DATA
#define TIME_KEEP
void generateData(int K, int input_data_num, int test_data_num)
{
    // 生成K和数据
    int all_data_num = input_data_num + test_data_num;
    int **data = new int *[all_data_num];
    for (int i = 0; i < all_data_num; i++)
    {
        data[i] = new int[K];
    }
    for (int i = 0; i < all_data_num; i++)
    {
        for (int j = 0; j < K; j++)
        {
            data[i][j] = rand() % 101; // 生成不超过100的随机整数
        }
    }
    // 写入数据到文本文件
    ofstream outfile("data.txt");
    if (outfile.is_open())
    {
        outfile << K << endl; // 先写入K
        for (int i = 0; i < all_data_num; i++)
        {
            for (int j = 0; j < K; j++)
            {
                outfile << data[i][j] << " "; // 写入数据，用空格分隔
            }
            outfile << endl; // 每组数据写完后换行
        }
        outfile.close();
        cout << "Data has been wrriten to data.txt." << endl;
    }
    else
    {
        cout << "Can not open file data.txt." << endl;
    }
}

void readDataToVec(int K_read, int input_data_num, int test_data_num, vector<vector<int>> &input_data, vector<vector<int>> &test_data)
{
    // 生成K和数据
    int all_data_num = input_data_num + test_data_num;
    int **data_read = new int *[all_data_num];
    for (int i = 0; i < all_data_num; i++)
    {
        data_read[i] = new int[K_read];
    }
    // 从文本文件中读取数据
    ifstream infile("data.txt");
    if (infile.is_open())
    {
        infile >> K_read; // 先读取K
        for (int i = 0; i < all_data_num; i++)
        {
            vector<int> vec;
            for (int j = 0; j < K_read; j++)
            {
                int num;
                infile >> num;
                vec.push_back(num);
            }

            if (i < input_data_num)
                input_data.push_back(vec);
            else
                test_data.push_back(vec);
        }
        infile.close();
    }
    else
    {
        cout << "Can not open data.txt." << endl;
    }
}
void KDTree_search_nearest(int K, vector<vector<int>> &data1, vector<vector<int>> &data2, vector<vector<int>> &res)
{
#ifdef TIME_KEEP
    clock_t start, end1, end2;
    start = clock();
#endif
    KDTree<int> kdTree = KDTree<int>(K, data1);
#ifdef TIME_KEEP
    end1 = clock();
    cout << "The time to build KDTree is " << double(end1 - start) << " ms" << endl;
    cout << "Size of KDTree is " << kdTree.size() << endl;
#endif
    for (int i = 0; i < data2.size(); i++)
    {
        vector<int> nearest = kdTree.findNearestPoint(data2[i]);
        /*
        for (int num : nearest)
            cout << num << ' ';
        cout << endl;
        */
        res.push_back(nearest);
    }
#ifdef TIME_KEEP
    end2 = clock();
    cout << "The time to find the nearest neighbor points using KDTree-Search is " << double(end2 - end1) << " ms" << endl;
#endif
}
void Linear_search_nearest(int K, vector<vector<int>> &data1, vector<vector<int>> &data2, vector<vector<int>> &res)
{
#ifdef TIME_KEEP
    clock_t start, end;
    start = clock();
#endif
    for (int i = 0; i < data2.size(); i++)
    {
        double minDist = std::numeric_limits<double>::infinity();
        vector<int> *nearest;
        for (int j = 0; j < data1.size(); j++)
        {
            double dist = 0;
            for (int k = 0; k < K; k++)
            {
                dist += pow(data1[j][k] - data2[i][k], 2);
            }
            if (dist < minDist)
            {
                minDist = dist;
                nearest = &(data1[j]);
            }
        }
        res.push_back(*nearest);
        /*
        for (int num : *nearest)
            cout << num << ' ';
        cout << endl;
        */
    }
#ifdef TIME_KEEP
    end = clock();
    cout << "The time to find the nearest neighbor points using Linear-Search is " << double(end - start) << " ms" << endl;
#endif
}
double calDist(vector<int> &v1, vector<int> &v2, int K)
{
    double dist = 0;
    for (int i = 0; i < K; i++)
        dist += pow(v1[i] - v2[i], 2);
    return dist;
}
int main()
{
    srand(time(0)); // 设置随机数种子
    int K = 5, input_data_num = 500, test_data_num = 5000;
#ifdef GENERATE_DATA
    generateData(K, input_data_num, test_data_num);
#endif
    vector<vector<int>> input, test, res1, res2;
    readDataToVec(K, input_data_num, test_data_num, input, test);

    KDTree_search_nearest(K, input, test, res1);
    Linear_search_nearest(K, input, test, res2);

    assert(res1.size() == res2.size());
    for (int i = 0; i < res1.size(); i++)
    {
        assert(calDist(res1[i], test[i], K) == calDist(res2[i], test[i], K));
        /*
        bool error = false;
        for (int j = 0; j < K; j++)
        {
            if (res1[i][j] != res2[i][j])
                error = true;
        }
        if (error == true)
        {
            cout << "Error Occurred:" << endl;
            cout << "Target: ";
            for (int k = 0; k < K; k++)
                cout << test[i][k] << " ";
            cout << endl;
            cout << "Res1: ";
            for (int k = 0; k < K; k++)
                cout << res1[i][k] << " ";
            cout << "Dist: " << calDist(test[i], res1[i], K);
            cout << endl;
            cout << "Res2: ";
            for (int k = 0; k < K; k++)
                cout << res2[i][k] << " ";
            cout << "Dist: " << calDist(test[i], res2[i], K);
            cout << endl
                 << endl;
        }
        */
    }
}