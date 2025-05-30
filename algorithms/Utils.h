#pragma once
#include <vector>
#include "SortAlgo.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <random>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;
const int inf = 1e8;
typedef function<bool(int, int)> Comp;

extern vector<vector<int>> buckets;
extern vector<double> median_estimates;


void new_pred();

void generate_statistical_preds();

void index_to_rank();

void output_to_file(vector<string> names, vector<vector<vector<ll>>> results, string filename);

double get_time();

typedef long long ll;
pair<vector<ll>, vector<ll>> parsePopulationData(const std::string &filename, int old);

vector<vector<bool>> parseTennisRelation(const string &filename, int size);

struct Node {
    int value, priority, st, ed, depth, size;
    Node *left, *right, *parent;

    Node(int value) : value(value), priority(rand()), st(-inf), ed(inf), left(nullptr), right(nullptr), parent(nullptr), depth(1), size(1) {}
};
extern vector<Node *> ai_to_node;


class Treap {
public:
    Treap() : root(nullptr), rightmost(nullptr) {}

    Node* insert(int value, const Comp& dirty, const Comp& clean, int target = -1) {
        Node *newNode;
        tie(root, newNode) = insert(root, value, -inf, inf, dirty, clean, false, target);
        return newNode;
    }

    Node *root, *rightmost;

    int getDepth(Node *node) {
        if (node == nullptr)
            return 0;
        return node->depth;
    }
    int getPriority(Node *node) {
        if (node == nullptr)
            return -1e9;
        return node->priority;
    }

    Node *rotateRight(Node *y) {
        if (y == nullptr)
            return nullptr;
        Node *x = y->left;
        y->left = x->right;
        if (x->right != nullptr) {
            x->right->parent = y;
        }
        x->right = y;
        x->parent = y->parent;
        y->parent = x;

        x->ed = y->ed;
        y->st = x->value;

        y->depth = max(getDepth(y->left), getDepth(y->right)) + 1;
        x->depth = max(getDepth(x->left), getDepth(x->right)) + 1;
        return x;
    }

    Node *rotateLeft(Node *x) {
        if (x == nullptr)
            return nullptr;
        Node *y = x->right;
        x->right = y->left;
        if (y->left != nullptr) {
            y->left->parent = x;
        }
        y->left = x;
        y->parent = x->parent;
        x->parent = y;

        y->st = x->st;
        x->ed = y->value;

        x->depth = max(getDepth(x->left), getDepth(x->right)) + 1;
        y->depth = max(getDepth(y->left), getDepth(y->right)) + 1;
        
        return y;
    }

    pair<Node*, Node*> insert(Node *node, int value, int st, int ed, const Comp& dirty, const Comp& clean, bool isClean, int target) {
        if (node == nullptr) {
            Node *newNode = new Node(value);
            newNode->st = st;
            newNode->ed = ed;
            return make_pair(newNode, newNode);
        }
        // cerr << "inserting " << value << " " << node->value << " " << isClean << endl;

        if (node->value == target)
            isClean = true;

        Node *newNode;
        Comp comp = isClean ? clean : dirty;
        if (node->value == target)
            //comp is the inverse of dirty
            comp = [&](int a, int b) {return !dirty(a, b);};

        if (comp(value, node->value))
        {
            tie(node->left, newNode) = insert(node->left, value, st, node->value, dirty, clean, isClean, target);
            node->left->parent = node;

        } 
        else {
            tie(node->right, newNode) = insert(node->right, value, node->value, ed, dirty, clean, isClean, target);
            node->right->parent = node;
        }

        if (getPriority(node->left)> node->priority) {
            node = rotateRight(node);
        }
        if (getPriority(node->right) > node->priority) {
            node = rotateLeft(node);
        }
        // Update st and ed values
        if (node->left)     node->st = node->left->st;
        if (node->right)    node->ed = node->right->ed;

        assert(node!=nullptr);
        node->depth = max(getDepth(node->left), getDepth(node->right)) + 1;
        return make_pair(node, newNode);
    }

    

    void LinearOutput(vector<int> *ranking) const {LinearOutput(root, ranking);}
    static void LinearOutput(const Node *node, vector<int> *ranking) {
        if (node == nullptr) return;
        LinearOutput(node->left, ranking);
        ranking->push_back(node->value);
        LinearOutput(node->right, ranking);
    }

    static void LinearPrint(const Node *node) {
        if (node == nullptr) return;
        LinearPrint(node->left);
        cerr << node->value << " " << node->st << " " << node->ed << " par=" << (node->parent == nullptr ? -1 : node->parent->value) << endl;
        LinearPrint(node->right);
    }

    void deleteNode(Node *node) {
        if (!node)
            return ;
        Node *parent = node;
        deleteNodeHelper(node);
        while (parent) {
            // cerr << "???" << endl;
            parent->depth = max(getDepth(parent->left), getDepth(parent->right)) + 1;
            parent = parent->parent;
        }
    }

    void deleteNodeHelper(Node *node_to_be_deleted) {
        if (!node_to_be_deleted)
            return ;

        Node *parent = node_to_be_deleted->parent;

        if (!node_to_be_deleted->left || !node_to_be_deleted->right) {
            Node *temp = node_to_be_deleted->left ? node_to_be_deleted->left : node_to_be_deleted->right;

            if (parent) {
                if (parent->left == node_to_be_deleted)
                    parent->left = temp;
                else
                    parent->right = temp;

                if (temp)
                    temp->parent = parent;
                parent->depth = max(getDepth(parent->left), getDepth(parent->right)) + 1;
            } else {
                root = temp;
                if (temp)
                    temp->parent = nullptr;
            }
            delete node_to_be_deleted;
        }
        else {
            if (node_to_be_deleted->left->priority < node_to_be_deleted->right->priority) {
                node_to_be_deleted = rotateRight(node_to_be_deleted);
                deleteNodeHelper(node_to_be_deleted->right);

            } else {
                node_to_be_deleted = rotateLeft(node_to_be_deleted);
                deleteNodeHelper(node_to_be_deleted->left);
            }
            if (node_to_be_deleted)
                node_to_be_deleted->depth = max(getDepth(node_to_be_deleted->left), getDepth(node_to_be_deleted->right)) + 1;
        }

        if (parent) {
            if (parent->left)  parent->st = parent->left->st;
            if (parent->right) parent->ed = parent->right->ed;
        }
        return ;
    }
};


extern vector<Node*> nodes;
class ScapegoatTree {
public:
    ScapegoatTree() : root(nullptr), rightmost(nullptr) {}

    Node* insert(int value, const function<bool(int, int)>& dirty, const function<bool(int, int)>& clean, int target = -1) {
        Node *newNode;
        tie(root, newNode) = insert(root, value, -inf, inf, dirty, clean, false, target, false);
        return newNode;
    }
    Node* freeze_insert(int value, const function<bool(int, int)>& dirty, const function<bool(int, int)>& clean, int target = -1) {
        Node *newNode;
        tie(root, newNode) = insert(root, value, -inf, inf, dirty, clean, false, target, true);
        return newNode;
    }
    

    Node *root, *rightmost;

    int getDepth(Node *node) {
        if (node == nullptr)
            return 0;
        return node->depth;
    }

    int getSize(Node *node) {
        if (node == nullptr)
            return 0;
        return node->size;
    }

    Node *rotateRight(Node *y) {
        if (y == nullptr)
            return nullptr;
        Node *x = y->left;
        y->left = x->right;
        if (x->right != nullptr) {
            x->right->parent = y;
        }
        x->right = y;
        x->parent = y->parent;
        y->parent = x;

        x->ed = y->ed;
        y->st = x->value;

        y->depth = max(getDepth(y->left), getDepth(y->right)) + 1;
        y->size = getSize(y->left) + getSize(y->right) + 1;
        x->depth = max(getDepth(x->left), getDepth(x->right)) + 1;
        x->size = getSize(x->left) + getSize(x->right) + 1;
        return x;
    }

    Node *rotateLeft(Node *x) {
        if (x == nullptr)
            return nullptr;
        Node *y = x->right;
        x->right = y->left;
        if (y->left != nullptr) {
            y->left->parent = x;
        }
        y->left = x;
        y->parent = x->parent;
        x->parent = y;

        y->st = x->st;
        x->ed = y->value;

        x->depth = max(getDepth(x->left), getDepth(x->right)) + 1;
        x->size = getSize(x->left) + getSize(x->right) + 1;
        y->depth = max(getDepth(y->left), getDepth(y->right)) + 1;
        y->size = getSize(y->left) + getSize(y->right) + 1;
        
        return y;
    }

    const double ALPHA = 1.2;
    bool isScapegoat(Node *node) {
        if (!node) return false;
        double large = max(getDepth(node->left), getDepth(node->right));
        double small = min(getDepth(node->left), getDepth(node->right));
        return large >= ALPHA * (small + 1);
    
        // double largeSize = max(getSize(node->left), getSize(node->right));
        // double smallSize = min(getSize(node->left), getSize(node->right));
        // return largeSize > ALPHA * (smallSize + 1);
    }


    pair<Node*, Node*> insert(Node *node, int value, int st, int ed, const function<bool(int, int)>& dirty, const function<bool(int, int)>& clean, bool isClean, int target, bool isFreeze) {
        // cerr << "in insert" << endl;
        if (node == nullptr) {
            Node *newNode = new Node(value);
            newNode->st = st;
            newNode->ed = ed;
            // cerr << "in nullptr, completed insert" << endl;
            return make_pair(newNode, newNode);
        }

        if (node->value == target)
            isClean = true;

        function<bool(int, int)> comp = isClean ? clean : dirty;
        if (node->value == target)
            comp = [&](int a, int b) {return !dirty(a, b);};

        Node *newNode;
        if (comp(value, node->value)) {
            tie(node->left, newNode) = insert(node->left, value, st, node->value, dirty, clean, isClean, target, isFreeze);
            node->left->parent = node;
        } else {
            tie(node->right, newNode) = insert(node->right, value, node->value, ed, dirty, clean, isClean, target, isFreeze);
            node->right->parent = node;
        }

        node->depth = max(getDepth(node->left), getDepth(node->right)) + 1;
        node->size = getSize(node->left) + getSize(node->right) + 1;

        if (!isFreeze && isScapegoat(node)) {
            node = rebuild(node);
            updateAncestors(node);
        }

        if (node->left) node->st = node->left->st;
        if (node->right) node->ed = node->right->ed;

        assert(node != nullptr);
        return make_pair(node, newNode);
    }

    Node* rebuild(Node *node) {
        if (!node) return nullptr;
        // cerr << "rebuilding " << node->value << endl;
        nodes.resize(0);
        
        flatten(node, nodes);
        Node *par = node->parent;
        int st = node->st, ed = node->ed;
        Node *new_root = buildBalanced(nodes, 0, nodes.size() - 1, st, ed);
        new_root->parent = par;
        
        return new_root;
    }

    void flatten(Node *node, vector<Node*>& nodes) {
        if (!node) return;
        flatten(node->left, nodes);
        nodes.push_back(node);
        flatten(node->right, nodes);
        node->left = node->right = nullptr;
        node->depth = node->size = 1;
    }

    Node* buildBalanced(vector<Node*>& nodes, int l, int r, int st, int ed) {
        if (l > r) return nullptr;
        int mid = (l + r) / 2;
        Node *root = nodes[mid];
        root->st = st;
        root->ed = ed;
        root->left = buildBalanced(nodes, l, mid - 1, st, root->value);
        if (root->left) root->left->parent = root;
        root->right = buildBalanced(nodes, mid + 1, r, root->value, ed);
        if (root->right) root->right->parent = root;

        root->depth = max(getDepth(root->left), getDepth(root->right)) + 1;
        root->size = getSize(root->left) + getSize(root->right) + 1;
        return root;
    }

    void updateAncestors(Node *node) {
        if (!node)
            return;
        if (node->parent) {
            node->parent->depth = max(getDepth(node->parent->left), getDepth(node->parent->right)) + 1;
            node->parent->size = getSize(node->parent->left) + getSize(node->parent->right) + 1;
            updateAncestors(node->parent);
        }
    }

    void LinearOutput(vector<int> *ranking) const {LinearOutput(root, ranking);}
    static void LinearOutput(const Node *node, vector<int> *ranking) {
        if (node == nullptr) return;
        LinearOutput(node->left, ranking);
        ranking->push_back(node->value);
        LinearOutput(node->right, ranking);
    }

    static void LinearPrint(const Node *node) {
        if (node == nullptr) return;
        LinearPrint(node->left);
        cerr << node->value << " " << node->st << " " << node->ed << " par=" << (node->parent == nullptr ? -1 : node->parent->value) << endl;
        LinearPrint(node->right);
    }

    Node *find_kth_smallest(Node *node, int k) {
        if (node == nullptr) return nullptr;
        int left_size = getSize(node->left);
        if (left_size == k - 1) return node;
        if (left_size >= k) return find_kth_smallest(node->left, k);
        return find_kth_smallest(node->right, k - left_size - 1);
    }

    Node *find_kth_largest(Node *node, int k) {
        if (node == nullptr) return nullptr;
        int right_size = getSize(node->right);
        if (right_size == k - 1) return node;
        if (right_size >= k) return find_kth_largest(node->right, k);
        return find_kth_largest(node->left, k - right_size - 1);
    }

    Node* find_largest_small(Node *node, int value, const Comp &comp) {
        if (node == nullptr) return nullptr;
        if (comp(node->value, value)) {
            Node *res = find_largest_small(node->right, value, comp);
            if (res) return res;
            return node;
        }
        return find_largest_small(node->left, value, comp);
    }

    Node* find_smallest_large(Node *node, int value, const Comp &comp) {
        if (node == nullptr) return nullptr;
        if (comp(value, node->value)) {
            Node *res = find_smallest_large(node->left, value, comp);
            if (res) return res;
            return node;
        }
        return find_smallest_large(node->right, value, comp);
    }

    Node *merge(Node *left, Node *right) {
        if (left == nullptr) return right;
        if (right == nullptr) return left;
        if (left->depth > right->depth) {
            left->right = merge(left->right, right);
            left->right->parent = left;
            left->depth = max(getDepth(left->left), getDepth(left->right)) + 1;
            left->size = getSize(left->left) + getSize(left->right) + 1;
            return left;
        } else {
            right->left = merge(left, right->left);
            right->left->parent = right;
            right->depth = max(getDepth(right->left), getDepth(right->right)) + 1;
            right->size = getSize(right->left) + getSize(right->right) + 1;
            return right;
        }
    }

    void del(Node * r) {
        if (r == root) {
            root = merge(root->left, root->right);
            return ;
        }
        Node *p = r->parent;
        if (p->left == r) {
            p->left = merge(r->left, r->right);
            if (p->left) p->left->parent = p;
        } else {
            p->right = merge(r->left, r->right);
            if (p->right) p->right->parent = p;
        }
    }

    Node *getLCA(Node *x, Node *y) {
        if (x == y) return x;
        if (x->size < y->size) swap(x, y);
        return getLCA(x, y->parent);
    }

    bool compare(Node *x, Node *y) {
        /*
            return true if x is on the left of lca(x, y)
            if x is y's ancestor, return if y is on the right of lca(x, y)
        */
        // assert(x != y);
        if (x == y)
            return false;
        if (x == nullptr || y == nullptr)
            return false;
        Node *lca = getLCA(x, y);
        if (lca == x) {
            while (y != lca->left && y != lca->right)
                y = y->parent;
            if (y == lca->left)
                return false;
            else 
                return true;
        }
        while (x != lca->left && x != lca->right)
            x = x->parent;
        if (x == lca->left)
            return true;
        else 
            return false;
    }

    // --- Added helper methods for rank-based operations ---
    // Split the tree by rank. The left subtree contains the first k nodes (in-order).
    // Returns {leftSubtreeRoot, rightSubtreeRoot}.
    pair<Node*, Node*> split(Node *node, int k) {
        if (!node) return {nullptr, nullptr};
        int leftSize = getSize(node->left);
        if (k <= leftSize) {
            auto pr = split(node->left, k);
            node->left = pr.second;
            if (node->left) node->left->parent = node;
            if (pr.first)  pr.first->parent = nullptr; // detach
            // update metadata
            node->depth = max(getDepth(node->left), getDepth(node->right)) + 1;
            node->size  = getSize(node->left) + getSize(node->right) + 1;
            if (node->left) node->st = node->left->st;
            if (node->right) node->ed = node->right->ed;
            return {pr.first, node};
        } else {
            auto pr = split(node->right, k - leftSize - 1);
            node->right = pr.first;
            if (node->right) node->right->parent = node;
            if (pr.second) pr.second->parent = nullptr; // detach
            node->depth = max(getDepth(node->left), getDepth(node->right)) + 1;
            node->size  = getSize(node->left) + getSize(node->right) + 1;
            if (node->left) node->st = node->left->st;
            if (node->right) node->ed = node->right->ed;
            return {node, pr.second};
        }
    }

    // Convenience wrapper: split starting from current root
    pair<Node*, Node*> split_root(int k) { return split(root, k); }

    // Join left + mid + right where all keys in left < mid < right.
    Node* join(Node *left, Node *mid, Node *right) {
        if (mid) {
            mid->left = mid->right = nullptr;
            mid->parent = nullptr;
            mid->depth = mid->size = 1;
            mid->st = mid->ed = mid->value;
        }
        Node *tmp = merge(left, mid);
        Node *res = merge(tmp, right);
        if (res) res->parent = nullptr; // ensure root parent is null
        return res;
    }
};