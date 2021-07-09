#define _CRT_SECURE_NO_DEPRECATE
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;


const int EdgeNumber = 116;

const double EPSILON = 1e-6;


struct point {
    double x, y;
    point(const double _x = 0, const double _y = 0) :x(_x), y(_y) {};
    point(const point& _point) :x(_point.x), y(_point.y) {};
    const point operator+(const point& p)const {
        return point(x + p.x, y + p.y);
    };
    const point operator-(const point& p)const {
        return point(x - p.x, y - p.y);
    };
    const point operator*(const double time)const {
        return point(x * time, y * time);
    }
    const point operator/(const double time)const {
        return point(x / time, y / time);
    }
    const double dot(const point& p)const {
        return x * p.x + y * p.y;
    }
    point& operator+=(const point& p) {
        return (*this) = (*this) + p;
    }
    const double length()const {
        return sqrt(x * x + y * y);
    }
    const double operator*(const point& rhs)const;
};
typedef point vec;
const double vec::operator*(const point& rhs)const {
    return x * rhs.x + y * rhs.y;
}
const point center(const point _s, const point _e) {
    return (_s + _e) / 2.0;
}

struct Edge {
    point _start, _end;
    vector<point>_subdivs;
    vector<int>_compatibleEdges;
    Edge(const point& start, const point& end) :_start(start), _end(end) {}
    Edge(const Edge& rhs) :_start(rhs._start), _end(rhs._end),_subdivs(rhs._subdivs),_compatibleEdges(rhs._compatibleEdges){};
    const point MidPoint()const;
    void add_subdivisions();
    void add_spring_forces(vector<point>& forces_, double K_);
    void add_electrostatic_forces(vector<point>& forces_, Edge edge_, double eps);
    void update(vector<point>& forces_, double S_);
    const point vec()const;
};

const point Edge::vec()const {
    return _end - _start;
}

const point Edge::MidPoint()const {
    return center(_start, _end);
}

void Edge::add_subdivisions() {
    int oldSubdivsNum = (int)_subdivs.size();

    //cerr << "oldSubdivsNum " << oldSubdivsNum << endl;

    //还未开始细分
    if (oldSubdivsNum == 0)
        _subdivs.assign(1, center(_start, _end));
    else
    {
        //已经开始细分
        int newSubdivsNum = 2 * oldSubdivsNum, subdivIndex = 0, v1Index = -1, v2Index = 0;

        //当存在oldSubdivsNum个细分点的时候，整个初始线段其实是被分为oldSubdivsNum+1段小的线段，
        //segmentLength的含义是新的线段是旧的小的线段的倍数长度
        double segmentLength = double(oldSubdivsNum + 1) / double(newSubdivsNum + 1);
        vector<point> subdivisions(newSubdivsNum);


        point v1 = _start, v2 = _subdivs[0];

        //r为segmentLength，意义同上，后续r的意义为跳转到第几个的倍数

        //这是重点：求得只是细分点，而不是点与点之间的线段关系，一些点在运行的过程中有已经被逐渐的丢掉了
        double r = segmentLength;
        while (subdivIndex < newSubdivsNum)
        {
            //cerr << "subdivIndex " << subdivIndex << endl;
            subdivisions[subdivIndex] = v1 + (v2 - v1) * r;
            subdivIndex++;

            //不好解释
            //cerr << "r + segmentLength:" << r + segmentLength << " " << (r + segmentLength - 1 > 0) << " " << (r + segmentLength - 1) << endl;


            //double用于比较大于小于奇怪的bug，不能用 if(r + segmentLength > 1.0)进行判断，奇奇怪怪的
            if (r + segmentLength - 1.0 > EPSILON)
            {
                r = segmentLength - (1.0 - r);
                v1Index++;
                v2Index++;

               // cerr << "v1Index " << v1Index << endl;
               // cerr << "v2Index " << v2Index << endl;

                if (v1Index >= 0)
                    v1 = _subdivs[v1Index];
                if (v2Index < oldSubdivsNum)
                    v2 = _subdivs[v2Index];
                else
                    v2 = _end;
            }
            else
                r += segmentLength; //如果先前的小的线段的长度还够一个新的小的线段的大小
        }
        _subdivs = subdivisions;
    }
}


void Edge::add_spring_forces(vector<point>& forces_, double K_)
{
    int len = (int)_subdivs.size();

    //与论文的公式貌似有一点点的出入
    double kP = K_ / ((_end - _start).length() * double(len + 1));

    if (len == 1)
        forces_[0] += (_start + _end - _subdivs[0] * 2.0) * kP;
    else
    {
        // first division point
        forces_[0] += (_start + _subdivs[1] - _subdivs[0] * 2.0) * kP;
        // inner division points
        for (int i = 1; i < len - 1; i++)
            forces_[i] += (_subdivs[i - 1] + _subdivs[i + 1] - _subdivs[i] * 2.0) * kP;
        // last division point
        forces_[len - 1] += (_subdivs[len - 2] + _end - _subdivs[len - 1] * 2.0) * kP;
    }
}


void Edge::add_electrostatic_forces(vector<point>& forces_, Edge edge_, double epsilon_)
{
    int len = (int)_subdivs.size();
    double dlen;

    //cerr << "add_electrostatic_forces start "<<len << endl;
    for (int i = 0; i < len; i++)
    {
        auto dist = (edge_._subdivs[i] - _subdivs[i]);
        dlen = dist.length();
        if (dlen > epsilon_)
            forces_[i] += dist / dlen;
    }
    //cerr << "add_electrostatic_forces end " << endl;
}


void Edge::update(vector<point>& forces_, double S_)
{
    int len = (int)_subdivs.size();
    double flen = 0.0;
    for (int i = 0; i < len; i++)
    {
        flen = forces_[i].length();
        
        /*******************************************************************/
        //与合力的大小没有关系，只是在力的方向上移动固定的距离
        if (flen > EPSILON)
            _subdivs[i] += forces_[i] * S_ / flen;
    }
}

const point project(const point& point_, const point& lineStart_, const point& lineEnd_)
{
    double L = (lineStart_ - lineEnd_).length();
    double r = ((lineStart_.y - point_.y) * (lineStart_.y - lineEnd_.y)
        + (lineStart_.x - point_.x) * (lineStart_.x - lineEnd_.x)) / (L * L);
    return lineStart_ + (lineEnd_ - lineStart_) * r;
}

double edge_visibility(const Edge& edge1_, const Edge& edge2_)
{
    const point I0 = project(edge1_._start, edge2_._start, edge2_._end);
    const point I1 = project(edge1_._end, edge2_._start, edge2_._end);
    const point midI = center(I0, I1);
    const point midP = center(edge2_._start, edge2_._end);
    return std::max(0.0, 1.0 - 2.0 * (midP - midI).length() / (I0 - I1).length());
}

double angle_compatilibity(const Edge& edge1_, const Edge& edge2_)
{
    const vec v1 = edge1_.vec();
    const vec v2 = edge2_.vec();
    return fabs(v1 * v2 / (v1.length() * v2.length()));
}

double scale_compatibility(const Edge& edge1_, const Edge& edge2_)
{
    double l1 = edge1_.vec().length();
    double l2 = edge2_.vec().length();
    double lavg = (l1 + l2) / 2.0;

    //防止出现除以0的情况
    if (lavg > EPSILON)
        return 2.0 / (lavg / min(l1, l2) + max(l1, l2) / lavg);
    else
        return 0.0;
}

double position_compatibility(const Edge& edge1_, const Edge& edge2_)
{
    double lavg = (edge1_.vec().length() + edge2_.vec().length()) / 2.0;
    //防止出现除以0的情况
    if (lavg > EPSILON)
    {
        const vec mid1 = center(edge1_._start, edge1_._end);
        const vec mid2 = center(edge2_._start, edge2_._end);
        return lavg / (lavg + (mid1 - mid2).length());
    }
    else
        return 0.0;
}

double visibility_compability(const Edge& edge1_, const Edge& edge2_)
{
    return min(edge_visibility(edge1_, edge2_), edge_visibility(edge2_, edge1_));
}

struct Graph {
    vector<Edge>_edges;
    double _K;                                  // 全局弹簧刚度系数
    int _I;                                     // 当前迭代次数的循环回数
    int _iter;                                  // 当前迭代次数剩余的循环回数
    int _cycles;                                // 剩余的迭代次数
    double _compatibilityThreshold;

    double _edgeDistance;
    double _S;

    void build_compatibility_lists();
    int iterate();
    int update_cycle();
    void add_subvisions();
};


void Graph::build_compatibility_lists()
{

    int edgesNum = (int)_edges.size();
    double comp = 0.0;
    for (int i = 0; i < edgesNum; i++)
    {
        for (int j = i + 1; j < edgesNum; j++)
        {
            comp = angle_compatilibity(_edges[i], _edges[j]) * scale_compatibility(_edges[i], _edges[j]) * position_compatibility(_edges[i], _edges[j]) * visibility_compability(_edges[i], _edges[j]);
            if (comp >= _compatibilityThreshold)
            {
                _edges[i]._compatibleEdges.push_back(j);
                _edges[j]._compatibleEdges.push_back(i);
            }
        }
    }
}


int Graph::iterate()
{
    int edgesNum = (int)_edges.size();
    cerr << "iterator() edgesNum: " << edgesNum << endl;
    vector<vector<point> > forces(edgesNum, vector<point>((int)_edges[0]._subdivs.size(), point(0.0, 0.0)));

    // spring forces
    for (int i = 0; i < edgesNum; i++) {
        //cerr << "iterator() add_spring_forces " << i << endl;
        _edges[i].add_spring_forces(forces[i], _K);
    }
       

    // electrostatic forces
    for (int i = 0; i < edgesNum; i++)
    {
        int compatibleEdgesNum = (int)_edges[i]._compatibleEdges.size();
        for (int j = 0; j < compatibleEdgesNum; j++) {
            //cerr << "iterator() electrostatic forces " << i <<" "<<j << endl;
            _edges[i].add_electrostatic_forces(forces[i], _edges[_edges[i]._compatibleEdges[j]], _edgeDistance);
        }
    }

    // update edges
    for (int i = 0; i < edgesNum; i++)
        _edges[i].update(forces[i], _S);

    _iter--;
    return _iter;
}

int Graph::update_cycle()
{
    _S *= 0.5;
    _I = 2 * _I / 3;
    _iter = _I;
    _cycles--;
    return _cycles;
}

void Graph::add_subvisions()
{
    int edgesNum = (int)_edges.size();
    for (int i = 0; i < edgesNum; i++) {
        //cerr << "add_subdivisionsstart " << i << endl;
        _edges[i].add_subdivisions();
       // cerr << "add_subdivisions end " << i <<" subdivsNum "<<_edges[i]._subdivs.size() << endl;
    }
}


void OutputEdge(const point& start, const point& end) {
    //cout << "plt.plot([" << start.x << "," << end.x << "],[" << start.y << "," << end.y << "])" << '\n';
    cout << start.x << " " << end.x << " " << start.y << " " << end.y << '\n';
}

signed main() {

    freopen("C:\\Users\\pengfei\\Desktop\\FDEB\\edge.txt", "r", stdin);
    //freopen("C:\\Users\\pengfei\\Desktop\\FDEB\\log.txt", "w", stderr);
    freopen("C:\\Users\\pengfei\\Desktop\\result.txt", "w", stdout);

    Graph G;
    for (int _ = 0; _ < EdgeNumber; ++_) {
        point _start, _end;
        cin >> _start.x >> _start.y >> _end.x >> _end.y;
        G._edges.push_back(Edge(_start, _end));
    }
    G._cycles = 5;
    G._K = 1;
    G._S = 0.1;
    G._I = 50;
    G._compatibilityThreshold = 0.1;
    G._edgeDistance = 1e-4;

    G.build_compatibility_lists();
    G.add_subvisions();
    do
    {
        //cerr << "update_cycle end" << endl;
        while (G.iterate() > 0);
       // cerr << "iterate end" << endl;
        G.add_subvisions();
        //cerr << "add_subvisions end" << endl;
    } while (G.update_cycle() > 0);


    
    //cout << "import matplotlib.pyplot as plt" << endl;
    for (auto& edge : G._edges) {
        OutputEdge(edge._start, edge._subdivs[0]);
        for (int j = 0; j + 1 < int(edge._subdivs.size()); ++j) {
            OutputEdge(edge._subdivs[j], edge._subdivs[j + 1]);
        }
        OutputEdge(edge._subdivs.back(), edge._end);
    }
    //cout << "plt.show()" << endl;


    return 0;
}