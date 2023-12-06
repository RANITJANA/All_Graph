#include <bits/stdc++.h>
using namespace std;

// Graph Traversals

void BFS(int V, vector<int> adj[])
{
    int vis[V] = {0};
    queue<int> q;
    q.push(0);
    vis[0]++;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        cout << node << " ";
        for (auto it : adj[node])
        {
            if (!vis[it])
            {
                vis[it]++;
                q.push(it);
            }
        }
    }
    cout << endl;
}
void dfs(int start, vector<int> adj[], int vis[])
{
    vis[start] = 1;
    cout << start << " ";
    for (auto it : adj[start])
    {
        if (!vis[it])
            dfs(it, adj, vis);
    }
}
void DFS(int V, vector<int> adj[])
{
    int vis[V] = {0};
    for (int i = 0; i < V; i++)
    {
        if (!vis[i])
        {
            dfs(i, adj, vis);
        }
    }
    cout << endl;
}

//cucle detect

//assuming there is only one component
bool isCyclePresent_bfs(int V,vector<int>adj[]){
	int vis[V]; //here we will store parents
	for(int i=0;i<V;i++) vis[i]=-2;
	queue<int> q;
	q.push(0);
	vis[0]=-1;
	while(!q.empty()){
		int node=q.front();
		q.pop();
		for(auto it:adj[node]){
			if(vis[it]==-2){
				vis[it]=node;
				q.push(it);
			}
			else if(it!=vis[node]) return true;
		}
	}
	return false;
}


bool find_dfs(int node,int parent,vector<int>adj[],int vis[]){
	vis[node]=1;
	for(auto it:adj[node]){
		if(it==parent) continue;
		if(vis[it]) return true;
		if(find_dfs(it,node,adj,vis)) return true;
	}
	return false;
}

bool isCyclePresent_dfs(int V,vector<int>adj[]){
	int vis[V]={0};
	for(int i=0;i<V;i++){
		if(!vis[i] && find_dfs(i,-1,adj,vis)) return true;
	}
	return false;
}

// special type of graph

bool solveColor(int start, vector<int> adj[], int color[])
{
    color[start] = 0;
    queue<int> q;
    q.push(start);
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        for (auto it : adj[node])
        {
            if (color[it] == -1)
            {
                color[it] = !color[node];
            }
            else if (color[node] == color[it])
                return false;
        }
    }
    return true;
}
bool isBipartite(int V, vector<int> adj[]) // when 2 adjNodes are having two different color .. it is bipartite
{
    int color[V];
    for (int i = 0; i < V; i++)
        color[i] = -1;
    for (int i = 0; i < V; i++)
    {
        if (color[i] == -1)
        {
            if (!solveColor(i, adj, color))
                return false;
        }
    }
    return true;
}

// linear graph representaion

void topo(int start, vector<int> adj[], int vis[], stack<int> &s)
{
    vis[start] = 1;
    for (auto it : adj[start])
    {
        if (!vis[it])
            topo(it, adj, vis, s);
    }
    s.push(start);
}
void TopoLogicalSort_recursion(int V, vector<int> adj[])
{
    int vis[V] = {0};
    stack<int> s;
    for (int i = 0; i < V; i++)
    {
        if (!vis[i])
            topo(i, adj, vis, s);
    }
    while (!s.empty())
    {
        cout << s.top() << " ";
        s.pop();
    }
}
void TopoLogicalSort_KhansAlgo(int V, vector<int> adj[])
{
    int indegree[V] = {0};
    for (int i = 0; i < V; i++)
    {
        for (auto it : adj[i])
        {
            indegree[it]++;
        }
    }
    queue<int> q;
    for (int i = 0; i < V; i++)
    {
        if (indegree[i] == 0)
            q.push(i);
    }
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        cout << node << " ";
        for (auto it : adj[node])
        {
            indegree[it]--;
            if (indegree[it] == 0)
            {
                q.push(it);
            }
        }
    }
    cout << endl;
}

// minimum path cost from a given source algorithm

vector<int> DijkstraAlog(int V, vector<vector<int>> adj[], int S)
{
    vector<int> dist(V, 1e9);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, S}); // distance,node
    dist[S] = 0;
    while (!pq.empty())
    {
        int dis = pq.top().first;
        int node = pq.top().second;
        pq.pop();
        for (auto it : adj[node])
        {
            int adjNode = it[0];
            int edgeWeight = it[1];
            if (edgeWeight + dis < dist[adjNode])
            {
                dist[adjNode] = dis + edgeWeight;
                pq.push({dist[adjNode], adjNode});
            }
        }
    }
    return dist;
}
vector<int> BellmanFordAlgo(int V, vector<vector<int>> &edges, int S)
{
    vector<int> dist(V, 1e9);
    dist[S] = 0;
    for (int i = 0; i < V - 1; i++) // for worst case we have to traverse n-1 times.
    {
        for (auto it : edges)
        {
            int u = it[0];
            int v = it[1];
            int wt = it[2];
            if (dist[u] != 1e9 && dist[u] + wt < dist[v])
                dist[v] = dist[u] + wt;
        }
    }
    for (auto it : edges)
    {
        int u = it[0];
        int v = it[1];
        int wt = it[2];
        if (dist[u] != 1e9 && dist[u] + wt < dist[v])
        {
            cout << "Negative edge cycle is Present.\n";
            return {-1};
        }
    }
    return dist;
}

// multisource sortest algorithm

void FloydWarShallAlgo(vector<vector<int>> &matrix)
{
    // in place algo.
    //  if unreahcable .. the value is -1.
    int n = matrix.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (matrix[i][j] == -1)
                matrix[i][j] = 1e8;
        }
    }
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                matrix[i][j] = min(matrix[i][j], matrix[i][k] + matrix[k][j]);
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (matrix[i][j] == 1e8)
                matrix[i][j] = -1;
        }
    }
}

// minimum spanning tree algoritm

pair<vector<vector<int>>, int> prim_s_Alog(int V, vector<vector<int>> adj[])
{
    vector<vector<int>> mst;
    int sum = 0;
    int vis[V] = {0};
    priority_queue<pair<int, pair<int, int>>, vector<pair<int, pair<int, int>>>, greater<pair<int, pair<int, int>>>> pq;
    pq.push({0, {0, -1}}); // edgeweight,node,parent
    while (!pq.empty())
    {
        int wt = pq.top().first;
        int node = pq.top().second.first;
        int parent = pq.top().second.second;
        pq.pop();
        if (vis[node])
            continue;
        if (parent != -1)
            mst.push_back({parent, node});
        sum += wt;
        vis[node] = 1;
        for (auto it : adj[node])
        {
            int adjNode = it[0];
            int edgeWeight = it[1];
            if (!vis[adjNode])
                pq.push({edgeWeight, {adjNode, node}});
        }
    }
    return {mst, sum};
}

// another data structure to be know for kruskal's algo.

// union by rank
class DisjointSet_ByRank
{
    vector<int> rank, parent;

public:
    DisjointSet_byRank(int n)
    {
        rank.resize(n + 1, 0);
        parent.resize(n + 1);
        for (int i = 0; i <= n; i++)
            parent[i] = i;
    }
    int findUpar(int node)
    {
        if (parent[node] == node)
            return node;
        return parent[node] = findUpar(parent[node]);
    }
    void unionByRank(int u, int v)
    {
        int ulp_u = findUpar(u);
        int ulp_v = findUpar(v);
        if (ulp_u == ulp_v)
            return;
        if (rank[ulp_u] < rank[ulp_v])
        {
            parent[ulp_u] = ulp_v;
        }
        else if (rank[ulp_v] < rank[ulp_u])
        {
            parent[ulp_v] = ulp_u;
        }
        else
        {
            parent[ulp_v] = ulp_u;
            rank[ulp_u]++;
        }
    }
};
class DisjointSet_BySize
{
    vector<int> size, parent;

public:
    DisjointSet_BySize(int n)
    {
        size.resize(n + 1, 1);
        parent.resize(n + 1);
        for (int i = 0; i <= n; i++)
            parent[i] = i;
    }
    int findUpar(int node)
    {
        if (node == parent[node])
            return node;
        return parent[node] = findUpar(parent[node]);
    }
    void unionBySize(int u, int v)
    {
        int ulp_u = findUpar(u);
        int ulp_v = findUpar(v);
        if (ulp_u == ulp_v)
            return;
        if (size[ulp_u] < size[ulp_v])
        {
            parent[ulp_u] = ulp_v;
            size[ulp_v] += size[ulp_u];
        }
        else
        {
            parent[ulp_v] = ulp_u;
            size[ulp_u] += size[ulp_v];
        }
    }
};
bool cmp(vector<int> &a, vector<int> &b)
{
    return a[0] <= b[0];
}
int kruskalAlgo(int V, vector<vector<int>> &edges) // edge format is : weight - u - v
{
    sort(edges.begin(), edges.end(), cmp);
    DisjointSet_BySize ds(V);
    int sum = 0;
    for (auto it : edges)
    {
        int wt = it[0];
        int u = it[1];
        int v = it[2];
        if (ds.findUpar(u) != ds.findUpar(v))
        {
            ds.unionBySize(u, v);
            sum += wt;
        }
    }
    return sum;
}

// Trie Data Structure
struct node
{
    node *links[26];
    bool flag = false;

    // instant member function
    bool isLetterPresent(char ch)
    {
        return links[ch - 'a'] != NULL;
    }
    void setLetter(char ch)
    {
        node *newnode = new node();
        links[ch - 'a'] = newnode;
    }
    node *getNext(char ch)
    {
        return links[ch - 'a'];
    }
    void setFlag()
    {
        flag = true;
    }
};
class Trie
{
    node *root;

public:
    Trie()
    {
        root = new node();
    }
    void insert(string word)
    {
        int n = word.size();
        node *temp = root;
        for (int i = 0; i < n; i++)
        {
            if (!temp->isLetterPresent(word[i]))
            {
                temp->setLetter(word[i]);
            }
            temp = temp->getNext(word[i]);
        }
        temp->setFlag();
    }
    bool search(string word)
    {
        node *temp = root;
        for (int i = 0; i < word.size(); ++i)
        {
            if (!temp->isLetterPresent(word[i]))
                return false;
            temp = temp->getNext(word[i]);
        }
        return temp->flag;
    }
    bool startWith(string prefix)
    {
        node *temp = root;
        for (int i = 0; i < prefix.size(); ++i)
        {
            if (!temp->isLetterPresent(prefix[i]))
                return false;
            temp = temp->getNext(prefix[i]);
        }
        return true;
    }
};
int main()
{
    /* Trie t;
     t.insert("ranit");
     t.insert("raitp");
     cout << t.search("ranit") << endl
          << t.search("ran") << endl
          << t.startWith("rai") << endl;
     */

    int m, n;
    cin >> m >> n;
    vector<vector<int>> edges;
    for (int i = 0; i < n; i++)
    {
        int wt, u, v;
        cin >> wt >> u >> v;
        edges.push_back({wt, u, v});
    }
    cout << kruskalAlgo(m, edges) << endl;

    return 0;
}
