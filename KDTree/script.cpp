#include <iostream>
#include <vector>
using namespace std;
int main()
{
    vector<int> v;
    v.push_back(1);
    v.erase(v.begin(), v.end());
    cout << v.size() << endl;
}