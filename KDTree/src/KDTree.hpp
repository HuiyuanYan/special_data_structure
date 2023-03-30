#include <vector>
#include <cmath>
#include <algorithm>
template <typename T>
struct KDNode
{
    std::vector<T> point;
    KDNode *left;
    KDNode *right;
    KDNode *father;
    bool exist;
    long long axis;
    long long size;
    KDNode(std::vector<T> point) : point(point), left(nullptr), right(nullptr), father(nullptr), axis(-1), exist(true), size(1){};
    void setAxis(long long axis) { this->axis = axis; }
    ~KDNode()
    {
        point.clear();
        left = right = nullptr;
        father = nullptr;
    }
};

template <typename T>
class KDTree
{
public:
    KDTree(long long K) : m_K(K), m_size(0), m_root(nullptr), m_rebuild(nullptr), m_if_rebuild(false) {}
    KDTree(long long K, std::vector<std::vector<T>> &points);
    ~KDTree();
    void printTree()
    {
        this->_recursivePrintTree(this->m_root, 0);
    }
    std::vector<T> findNearestPoint(std::vector<T> &target);
    std::vector<std::vector<T>> findKNearestPoints(std::vector<T> &target, long long &k);
    std::vector<std::vector<T>> findKDistPoints(std::vector<T> &target, double &dist_k);
    bool insertPoint(std::vector<T> &point);
    bool removePoint(std::vector<T> &point);
    bool exist(std::vector<T> &point);
    void clear();
    long long size() { return this->m_size; }
    long long K() { return this->m_K; }

private:
    void _recursivePrintTree(KDNode<T> *node, int depth);
    double _distance(std::vector<T> &data1, std::vector<T> &data2);
    void _setAxis(std::vector<KDNode<T> *> &nodes, long long l, long long r);
    long long _findMaxDistIdx(std::vector<double> dist);
    bool _is_cut(std::vector<T> &point1, std::vector<T> &point2, double r, long long axis);
    bool _dataCompare(std::vector<T> &data1, std::vector<T> &data2, long long axis);
    static bool _nodePtrCompare(KDNode<T> *node1, KDNode<T> *node2);
    bool _nodeCompare(KDNode<T> &node1, KDNode<T> &node2);

    void _get_nodes(KDNode<T> *root, std::vector<KDNode<T> *> &nodes);
    bool _recursiveInsertPoint(KDNode<T> *node, std::vector<T> &point);
    KDNode<T> *_findPoint(std::vector<T> &point);
    void _recursiveFindKDistNodes(KDNode<T> *node, std::vector<T> &target, double &dist_k, std::vector<KDNode<T> *> &k_dist_nodes);
    void _recursiveFindNearstPoint(KDNode<T> *node, std::vector<T> &target, double &minDist, KDNode<T> *&nearest);
    void _recursiveFindKNearstPoints(KDNode<T> *node, std::vector<T> &target, long long &k, std::vector<double> &k_nearest_dist, std::vector<std::vector<T>> &k_nearest);
    KDNode<T> *_recursiveBuildTree(std::vector<KDNode<T> *> &nodes, long long l, long long r);
    void _recursiveDestroyTree(KDNode<T> *node);

private:
    long long m_K;
    KDNode<T> *m_root;
    long long m_size;
    static float s_alpha;
    bool m_if_rebuild;
    KDNode<T> *m_rebuild;
};

template <typename T>
float KDTree<T>::s_alpha = 0.75;

template <typename T>
inline KDTree<T>::KDTree(long long K, std::vector<std::vector<T>> &points)
{
    this->m_K = K;
    std::vector<KDNode<T> *> nodes;
    for (std::vector<T> point : points)
    {
        KDNode<T> *node = new KDNode<T>(point);
        nodes.push_back(node);
    }
    this->m_root = this->_recursiveBuildTree(nodes, 0, nodes.size() - 1);
    this->m_rebuild = nullptr;
    this->m_if_rebuild = false;
    this->m_size = points.size();
}

template <typename T>
inline KDTree<T>::~KDTree()
{
    this->_recursiveDestroyTree(this->m_root);
    this->m_root = nullptr;
}

template <typename T>
inline std::vector<T> KDTree<T>::findNearestPoint(std::vector<T> &target)
{
    KDNode<T> *nearest = nullptr;
    double minDist = std::numeric_limits<double>::max();
    if (this->m_root == nullptr)
        return std::vector<T>();
    this->_recursiveFindNearstPoint(this->m_root, target, minDist, nearest);
    return nearest->point;
}

template <typename T>
inline std::vector<std::vector<T>> KDTree<T>::findKNearestPoints(std::vector<T> &target, long long &k)
{
    std::vector<std::vector<T>> k_nearest;
    std::vector<double> k_nearest_dist;
    this->_recursiveFindKNearstPoints(this->m_root, target, k, k_nearest_dist, k_nearest);
    return k_nearest;
}

template <typename T>
inline std::vector<std::vector<T>> KDTree<T>::findKDistPoints(std::vector<T> &target, double &dist_k)
{
    std::vector<KDNode<T> *> k_dist_nodes;
    this->_recursiveFindKDistNodes(this->m_root, dist_k, k_dist_nodes);
    std::vector<std::vector<T>> k_dist_points;
    for (KDNode<T> *node : k_dist_nodes)
    {
        k_dist_nodes.push_back(node->point);
    }
    return k_dist_points;
}

template <typename T>
inline bool KDTree<T>::insertPoint(std::vector<T> &point)
{
    KDNode<T> *node = this->_findPoint(point);
    if (node != nullptr)
        if (node->exist)
            return false;
        else
        {
            this->m_size += 1;
            return true;
        }

    assert(point.size() >= this->m_K);
    bool ret = false;
    if (this->m_root == nullptr)
    {
        this->m_root = KDNode<T>(point);
        m_root->axis = 0;
        ret = true;
    }
    else
    {
        ret = this->_recursiveInsertPoint(this->m_root, point);
    }
    if (this->m_if_rebuild == true)
    {
        std::vector<KDNode<T> *> nodes;
        this->_get_nodes(this->m_rebuild, nodes);
        if (this->m_rebuild->father == nullptr) // root
        {
            this->m_root = this->_recursiveBuildTree(nodes, 0, nodes.size() - 1);
        }
        else
        {
            KDNode<T> *father = this->m_rebuild->father;
            KDNode<T> *new_child = this->_recursiveBuildTree(nodes, 0, nodes.size() - 1);
            if (father->left == this->m_rebuild)
                father->left = new_child;
            if (father->right == this->m_rebuild)
                father->right = new_child;
        }
        this->m_if_rebuild = false;
        this->m_rebuild = nullptr;
    }
    return ret;
}

template <typename T>
inline bool KDTree<T>::removePoint(std::vector<T> &point)
{
    KDNode<T> *node = this->_findPoint(point);
    if (node != nullptr)
    {
        if (node->exist)
        {
            node->exist = false;
            this->m_size -= 1;
            return true;
        }
    }
    return false;
}

template <typename T>
inline bool KDTree<T>::exist(std::vector<T> &point)
{
    if (this->_findPoint(point) != nullptr)
        return true;
    return false;
}

template <typename T>
inline void KDTree<T>::clear()
{
    this->_recursiveDestroyTree(this->m_root);
    this->m_root = nullptr;
    this->m_size = 0;
}

template <typename T>
inline void KDTree<T>::_recursivePrintTree(KDNode<T> *node, int depth)
{
    if (node == nullptr)
        return;
    for (int i = 0; i < depth; i++)
        printf("\t");

    for (int i = 0; i < node->point.size(); i++)
        printf("%d ", node->point[i]);
    printf("\n");
    this->_recursivePrintTree(node->left, depth + 1);
    this->_recursivePrintTree(node->right, depth + 1);
}

template <typename T>
inline double KDTree<T>::_distance(std::vector<T> &data1, std::vector<T> &data2)
{
    double dist = 0;
    for (int i = 0; i < this->m_K; i++)
        dist += std::pow(data1[i] - data2[i], 2);
    return dist;
}

template <typename T>
inline void KDTree<T>::_setAxis(std::vector<KDNode<T> *> &nodes, long long l, long long r)
{
    if (l >= r)
        return;

    double max_variance = 0;
    long long axis_chosen = 0;
    for (int i = 0; i < this->m_K; i++)
    {
        double mean_value = 0;
        double variance = 0;
        for (int j = l; j <= r; j++)
            mean_value += nodes[j]->point[i];
        mean_value = mean_value / nodes.size();
        for (int j = l; j <= r; j++)
            variance += std::pow(nodes[j]->point[i] - mean_value, 2);
        if (variance > max_variance)
        {
            max_variance = variance;
            axis_chosen = i;
        }
    }
    // set axis
    for (int i = l; i <= r; i++)
        nodes[i]->axis = axis_chosen;
    return;
}

template <typename T>
inline long long KDTree<T>::_findMaxDistIdx(std::vector<double> dist)
{
    long long max_dist_idx = 0;
    for (int i = 1; i < dist.size(); i++)
        if (dist[i] > dist[max_dist_idx])
            max_dist_idx = i;
    return max_dist_idx;
}

template <typename T>
inline bool KDTree<T>::_is_cut(std::vector<T> &point1, std::vector<T> &point2, double r, long long axis)
{
    return std::fabs(point1[axis] - point2[axis]) < r;
}

template <typename T>
inline bool KDTree<T>::_dataCompare(std::vector<T> &data1, std::vector<T> &data2, long long axis)
{
    return data1[axis] < data2[axis];
}

template <typename T>
inline bool KDTree<T>::_nodePtrCompare(KDNode<T> *node1, KDNode<T> *node2)
{
    return node1->point[node1->axis] < node2->point[node1->axis];
}

template <typename T>
inline bool KDTree<T>::_nodeCompare(KDNode<T> &node1, KDNode<T> &node2)
{
    return node1.point[node1.axis] < node2.point[node1.axis];
}

template <typename T>
inline void KDTree<T>::_get_nodes(KDNode<T> *root, std::vector<KDNode<T> *> &nodes)
{
    nodes.clear();
    if (root == nullptr)
        return;
    std::vector<KDNode<T> *> st;

    if (root->exist)
        st.push_back(root);
    KDNode<T> *top;
    while (!st.empty)
    {
        top = st.back();
        st.pop_back();

        if (top->right != nullptr && top->right->exist)
            st.push_back(top->right);
        if (top->left != nullptr && top->left->exist)
            st.push_back(top->left);
        nodes.push_back(top);
    }
    return;
}

template <typename T>
inline bool KDTree<T>::_recursiveInsertPoint(KDNode<T> *node, std::vector<T> &point)
{
    bool ret = false;
    if (this->_dataCompare(node->point, point, node->axis) == true)
    {
        if (node->left == nullptr)
        {
            KDNode<T> *new_node = new KDNode<T>(point);
            new_node->axis = (node->axis + 1) % this->m_K;
            new_node->father = node;
            node->left = new_node;
            this->m_size += 1;
            ret |= true;
        }
        else
        {
            ret |= this->_recursiveInsertPoint(node->left, point);
        }
    }
    else
    {
        if (node->right == nullptr)
        {
            KDNode<T> *new_node = new KDNode<T>(point);
            new_node->axis = (node->axis + 1) % this->m_K;
            new_node->father = node;
            node->right = new_node;
            this->m_size += 1;
            ret |= true;
        }
        else
        {
            ret |= this->_recursiveInsertPoint(node->right, point);
        }
    }
    if ((node->left->size * this->s_alpha >= node->size) | (node->right->size * this->s_alpha >= node->size))
    {
        this->m_if_rebuild |= true;
        this->m_rebuild = node;
    }
    if (ret == true)
        node->size += 1;
    return ret;
}

/*Just find the point, don't care if it exists.*/
template <typename T>
inline KDNode<T> *KDTree<T>::_findPoint(std::vector<T> &point)
{
    std::vector<KDNode<T> *> st;
    if (this->m_root != nullptr && this->m_root->exist)
        st.push_back(m_root);
    KDNode<T> *top, ret;
    while (!st.empty())
    {
        top = st.back();
        st.pop_back();
        bool is_same = false;

        for (long long i = 0; i < this->m_K; i++)
            if (top->point[i] != point[i])
            {
                is_same = true;
                break;
            }
        if (is_same == true)
            return top;
        if (this->_dataCompare(top->point, point, top->axis) && top->left != nullptr)
            st.push_back(top->left);
        else if (this->_dataCompare(top->point, point, top->axis) && top->left != nullptr)
            st.push_back(top->right);
    }
    return nullptr;
}

template <typename T>
inline void KDTree<T>::_recursiveFindKDistNodes(KDNode<T> *node, std::vector<T> &target, double &dist_k, std::vector<KDNode<T> *> &k_dist_nodes)
{
    if (node == nullptr)
    {
        return;
    }
    if (this->_dataCompare(target, node->point, node->axis) == true)
        this->_recursiveFindKDistNodes(node->left, target, dist_k, k_dist_nodes);
    else
        this->_recursiveFindKDistNodes(node->right, target, dist_k, k_dist_nodes);

    if (node->exist)
    {
        double dist = this->_distance(target, node->point);
        if (dist <= dist_k)
        {
            k_dist_nodes.push_back(node);
        }
    }
    if (this->_is_cut(target, node->point, dist_k, node->axis))
    {
        if (this->_dataCompare(target, node->point, node->axis) == true)
            this->_recursiveFindKDistNodes(node->right, target, dist_k, k_dist_nodes);
        else
            this->_recursiveFindKDistNodes(node->left, target, dist_k, k_dist_nodes);
    }
}

template <typename T>
inline void KDTree<T>::_recursiveFindNearstPoint(KDNode<T> *node, std::vector<T> &target, double &minDist, KDNode<T> *&nearest)
{
    if (node == nullptr)
    {
        return;
    }
    if (this->_dataCompare(target, node->point, node->axis) == true)
        this->_recursiveFindNearstPoint(node->left, target, minDist, nearest);
    else
        this->_recursiveFindNearstPoint(node->right, target, minDist, nearest);

    if (node->exist)
    {
        double dist = this->_distance(target, node->point);
        if (dist < minDist)
        {
            minDist = dist;
            nearest = node;
        }
    }

    if (this->_is_cut(target, node->point, minDist, node->axis))
    {
        if (this->_dataCompare(target, node->point, node->axis) == true)
            this->_recursiveFindNearstPoint(node->right, target, minDist, nearest);
        else
            this->_recursiveFindNearstPoint(node->left, target, minDist, nearest);
    }
}

template <typename T>
inline void KDTree<T>::_recursiveFindKNearstPoints(KDNode<T> *node, std::vector<T> &target, long long &k, std::vector<double> &k_nearest_dist, std::vector<std::vector<T>> &k_nearest)
{
    if (node == nullptr)
    {
        return;
    }
    if (this->_dataCompare(target, node->point, node->axis) == true)
        this->_recursiveFindKNearstPoints(node->left, target, k, k_nearest_dist, k_nearest);
    else
        this->_recursiveFindKNearstPoints(node->right, target, k, k_nearest_dist, k_nearest);

    if (node->exist)
    {
        double dist = this->_distance(target, node->point);
        if (k_nearest_dist.size() < k)
        {
            k_nearest_dist.push_back(dist);
            k_nearest.push_back(node->point);
        }
        else
        {
            long long max_dist_idx = this->_findMaxDistIdx(k_nearest_dist);
            if (dist < k_nearest_dist[max_dist_idx])
            {
                k_nearest_dist[max_dist_idx] = dist;
                k_nearest[max_dist_idx] = node->point;
            }
        }
    }
    long long max_dist_idx = this->_findMaxDistIdx(k_nearest_dist);
    if (this->_is_cut(target, node->point, k_nearest_dist[max_dist_idx], node->axis))
    {
        if (this->_dataCompare(target, node->point, node->axis) == true)
            this->_recursiveFindKNearstPoints(node->right, target, k, k_nearest_dist, k_nearest);
        else
            this->_recursiveFindKNearstPoints(node->right, target, k, k_nearest_dist, k_nearest);
    }
}

template <typename T>
inline KDNode<T> *KDTree<T>::_recursiveBuildTree(std::vector<KDNode<T> *> &nodes, long long l, long long r)
{
    if (l > r)
        return nullptr;
    long long mid = (l + r) / 2;
    this->_setAxis(nodes, l, r);
    std::nth_element(nodes.begin() + l, nodes.begin() + mid, nodes.begin() + r + 1, _nodePtrCompare);
    KDNode<T> *midNode = nodes[mid];
    midNode->left = this->_recursiveBuildTree(nodes, l, mid - 1);
    midNode->right = this->_recursiveBuildTree(nodes, mid + 1, r);
    if (midNode->left != nullptr)
    {
        midNode->left->father = midNode;
        midNode->size += midNode->left->size;
    }
    if (midNode->right != nullptr)
    {
        midNode->right->father = midNode;
        midNode->size += midNode->right->size;
    }
    return midNode;
}

template <typename T>
inline void KDTree<T>::_recursiveDestroyTree(KDNode<T> *node)
{
    if (node == nullptr)
        return;
    this->_recursiveDestroyTree(node->left);
    this->_recursiveDestroyTree(node->right);
    delete node;
}
