// Randomly generate trees with n vertices and weighted edges (that represent Euclidean distance in a unit object)
// Calculate the average weight of the MST tree using Kruskal's algorithm

#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <time.h>

// random number generator
std::random_device rand_graph;
std::mt19937 gen(rand_graph());
std::uniform_real_distribution<float> randNum(0.0, 1.0);

// distance formula helper function
float distance(std::pair<std::pair<float, float>, std::pair<float, float>> p1,
               std::pair<std::pair<float, float>, std::pair<float, float>> p2)
{
    float dx = p2.first.first - p1.first.first;
    float dy = p2.first.second - p1.first.second;
    float dz = p2.second.first - p1.second.first;
    float da = p2.second.second - p1.second.second;

    return std::sqrt(dx * dx + dy * dy + dz * dz + da * da);
}


// maximum number of vertices
const int MAXN = 262144; 

//initializes parent and rank arrays
int parent[MAXN];
int rank[MAXN];

// kruskal's helper functions
void makeSet(int x)
{
    parent[x] = x;
    rank[x] = 0;
}

int find(int x)
{
    if (parent[x] != x)
    {
        parent[x] = find(parent[x]);
    }
    return parent[x];
}

int link(int x, int y)
{
    if (x == y)
    {
        return y;
    }
    if (rank[x] > rank[y])
    {
        std::swap(x, y);
    }
    if (rank[x] == rank[y])
    {
        rank[y]++;
    }
    parent[x] = y;
    return y;
}

void unionFind(int x, int y)
{
    link(find(x), find(y));
}

// compare helper function, helps us sort our list of edges in increasing order
bool compare(std::pair<float, std::pair<float, float>> a,
             std::pair<float, std::pair<float, float>> b)
{
    return a.first < b.first;
}

// implementation of kruskal's algorithm
std::pair<float, float> kruskals(std::vector<std::pair<float, std::pair<float, float>>> edges, int n)
{

    // sort the edges by weight in increasing order
    std::sort(edges.begin(), edges.end(), compare);

    // populate the parent & rank arrays
    for (int i = 0; i < n; i++)
    {
        makeSet(i);
    }

    // create the mst
    std::vector<std::pair<float, std::pair<int, int>>> mst;

    int numEdges = 0;
    for (int i = 0; i < edges.size(); i++)
    {
        int u = edges[i].second.first;
        int v = edges[i].second.second;
        float weight = edges[i].first;
        if (find(u) != find(v))
        {
            unionFind(u, v);
            mst.push_back(std::make_pair(weight, std::make_pair(u, v)));
            numEdges++;
            if (numEdges == n - 1)
            {
                break;
            }
        }
    }

    // calculates the weight of the whole MST and the weight of the max edge
    float mstWeight = 0.0;
    float maxEdge = 0.0;
    for (int i = 0; i < mst.size(); i++)
    {
        if (mst[i].first > maxEdge)
        {
            maxEdge = mst[i].first;
        }
        mstWeight += mst[i].first;
    }

    return {mstWeight, maxEdge};
}

// function for generating random graphs
std::vector<std::pair<float, std::pair<float, float>>> getGraph(int n, int dim)
{
    // creates vector of pairs which represents vertices
    std::vector<std::pair<std::pair<float, float>, std::pair<float, float>>> vertices(n);
    
    //checks dimension inputted
    for (int i = 0; i < n; i++)
    {
        float x = 0.0;
        float y = 0.0;
        float z = 0.0;
        float a = 0.0;

        if (dim == 0)
        {
            break;
        }
        if (dim >= 2)
        {
            x = randNum(gen);
            y = randNum(gen);
        }
        if (dim >= 3)
        {
            z = randNum(gen);
        }
        if (dim >= 4)
        {
            a = randNum(gen);
        }

        vertices[i] = {{x, y}, {z, a}};
    }

    /* optEdges stores our list of optimized edges in a vector of pairs 
    storing the weight and then another pair with the endpoints */
    std::vector<std::pair<float, std::pair<float, float>>> optEdges;

    float weight;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (dim == 0)
            {
                weight = randNum(gen);
            }
            else
            {
                weight = distance(vertices[i], vertices[j]);
            }
            if ((dim == 0 && weight > 0.005 && n >= 8192) ||
                (dim == 2 && weight > 0.02 && n >= 8192) ||
                (dim == 3 && weight > 0.15 && n >= 8192) ||
                (dim == 4 && weight > 0.2 && n >= 8192) ||
                (dim == 0 && weight > 0.001 && n >= 131072) ||
                (dim == 2 && weight > 0.005 && n >= 131072) ||
                (dim == 3 && weight > 0.035 && n >= 131072) ||
                (dim == 4 && weight > 0.07 && n >= 131072))
            {
                continue;
            }
            else
            {
                optEdges.push_back(std::make_pair(weight, std::make_pair(float(i), float(j))));
            }
        }
    }
    return optEdges;
}

int main(int argc, char *argv[])
{

    if (argc != 5)
    {
        printf("Usage: ./randmst flag numpoints numtrials dimension\n");
        return 1;
    }
    int flag = atoi(argv[1]);
    int numpoints = atoi(argv[2]);
    int numtrials = atoi(argv[3]);
    int dimension = atoi(argv[4]);

    if (dimension < 0 || dimension == 1 || dimension > 4)
    {
        printf("Valid Dimensions: 0, 2, 3, 4 \n");
        return 1;
    }

    // creates parent and rank arrays
    for (int i = 0; i < MAXN; i++)
    {
        parent[i] = i;
        rank[i] = 0;
    }

    srand((unsigned)time(NULL));
    clock_t begin = clock();

    float total = 0.0;
    float maxEdgeTrials = 0.0;
    
    // generates graphs based on numtrials and calculates 
    // Avg MST weight and weight of the max edge 
    std::pair<float, float> output;
    for (int t = 0; t < numtrials; t++)
    {
        auto graph = getGraph(numpoints, dimension);
        output = kruskals(graph, numpoints);
        total += output.first;
        if (output.second > maxEdgeTrials)
        {
            maxEdgeTrials = output.second;
        }
    }

    clock_t end = clock();
    printf("Avg MST weight for %i trials: %f\n", numtrials, total / (float)numtrials);
    printf("Max MST edge weight for %i trials: %f\n", numtrials, maxEdgeTrials);
    printf("Runtime: %f\n", float(end - begin) / CLOCKS_PER_SEC);

    return 0;
}