#include "../../inc/Algorithms/IsBiconnected.h"

IsBiconnected::IsBiconnected(Graph *g) : Algorithm(g) {

    intialiseAdjList();
}

void IsBiconnected::apply()
{
    bool isBiConnected;

    isBiConnected = isBC();

    if(isBiConnected == true){
        cout << "Graph is biconnected" << endl;
    }
    else
    {
       cout << "Graph is not biconnected" << endl;
    }

}

void IsBiconnected::intialiseAdjList(){
    for (int i = 0; i < graph->numVertices; ++i){
        //if (graph->edgeList[i][0] == graph->edgeList[i+1][0] )
        int v = graph->edgeList[i][0];
        int u = graph->edgeList[i][1];

                adj[v].push_back(u);
                adj[u].push_back(v);
    }
}

bool IsBiconnected::isBCUtil(int u, bool visited[], int disc[],int low[],int parent[])
{

    static int time = 0;

    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[u] = true;

    disc[u] = low[u] = ++time;

    // Go through all vertices aadjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of u

        // If v is not visited yet, then make it a child of u
        if (!visited[v])
        {
            children++;
            parent[v] = u;

            if (isBCUtil(v, visited, disc, low, parent))
               return true;

            low[u]  = min(low[u], low[v]);

            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == NULL && children > 1)
               return true;

            // (2) If u is not root and low value of one of its child is
            // more than discovery value of u.
            if (parent[u] != NULL && low[v] >= disc[u])
               return true;
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u]  = min(low[u], disc[v]);
    }
    return false;
}

// The main function that returns true if graph is Biconnected,
// otherwise false. It uses recursive function isBCUtil()
bool IsBiconnected::isBC()
{
    int V = graph->numVertices;
    //an array of booleans which denotes whether a vertex is visited or not
    bool *visited = new bool[V];
    //an array of V elements which stores the discovery time of every vertex
    int *disc = new int[V];
    //an array which stores the discovery time of the earliest discovered vertex
    int *low = new int[V];
    //an array which stores the parent of each vertex
    int *parent = new int[V];

    // Initialise parent and visited
    for (int i = 0; i < V; i++)
    {
        parent[i] = NULL;
        visited[i] = false;
    }

    // Call the recursive helper function to find if there is an articulation
    // point in given graph. We do DFS traversal starting from vertex 0
    if (isBCUtil(0, visited, disc, low, parent) == true)
        return false;

    // Now check whether the given graph is connected or not. An undirected
    // graph is connected if all vertices are reachable from any starting
    // point (0 as starting point)
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            return false;

    return true;
}
