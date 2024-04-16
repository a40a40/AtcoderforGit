//
#include <algorithm>
#include <bits/stdc++.h>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <stdio.h>
#include <string.h>
#include <utility>
#include <vector>
// #include <atcoder/all>

using namespace std;
#define ll long long
#define ull unsigned long long
#define rep(repindex, s, n)                                                    \
    for (ll(repindex) = (s); (repindex) < n; (repindex)++)
#define rrep(repindex, s, n)                                                   \
    for (ll(repindex) = (s); (repindex) >= n; (repindex)--)
#define bit(n, k) ((n >> k) & 1) /*nのk bit目*/

bool lrdef(ll x, ll l, ll r) { return (x >= l && x < r); }

// 凡例 : hwdef(h, w, 0, H, 0, W)
bool hwdef(ll h, ll w, ll hu, ll hd, ll wl, ll wr) {
    return lrdef(h, hu, hd) && lrdef(w, wl, wr);
}

const vector<vector<ll>> neib4 = {{0, -1}, {0, 1}, {-1, 0}, {1, 0}};
const vector<vector<ll>> neib8 = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0},
                                  {1, 0},   {-1, 1}, {0, 1},  {1, 1}};
struct Edge {
    ll num;
    ll from;
    ll to;
    ll cost;
};
void chmin(ll &a, ll b) {
    if (b < a)
        a = b;
}
void chmin2(ull &a, ull b) {
    if (b < a)
        a = b;
}
void chmax(ll &a, ll b) {
    if (b > a)
        a = b;
}
ll INF = 8e18;
ll mod = 998244353;
ll mod1e7 = 1000000007LL;
ll euc(ll a, ll b) {
    if (b > a)
        return euc(b, a);
    if (b == 0) {
        return a;
    }
    return euc(b, a % b);
}
pair<ll, ll> extgcd(ll a, ll b) {
    if (b == 0)
        return {1, 0};
    auto [bX, bY] = extgcd(b, a % b);
    return {bY, bX - bY * (a / b)};
}
ll lcm(ll a, ll b) { return a * b / euc(a, b); }
ll lcm(vector<ll> args) {
    ll ret = args[0];
    for (int i = 1; i < args.size(); i++) {
        ret = lcm(ret, args[i]);
    }
    return ret;
}

template <class T = long long> // formatCounter
ll LCS_base(vector<T> &vl, vector<T> &vr, vector<T> &ve) {
    if (vl.size() == 0 || vr.size() == 0)
        return 0;

    vector<vector<ll>> vt(vl.size() + 1, vector<ll>(vr.size() + 1, 0));
    rep(i, 1, vl.size() + 1) {
        rep(j, 1, vr.size() + 1) {
            if (vl[i - 1] == vr[j - 1])
                chmax(vt[i][j], vt[i - 1][j - 1] + 1);

            chmax(vt[i][j], vt[i][j - 1]);
            chmax(vt[i][j], vt[i - 1][j]);
        }
    }
    ll i = vl.size();
    ll j = vr.size();
    while (i > 0 && j > 0) {
        if (vl[i - 1] == vr[j - 1]) {
            ve.push_back(vl[i - 1]);
            i--;
            j--;
        } else if (vt[i][j] == vt[i][j - 1]) {
            j--;
        } else if (vt[i][j] == vt[i - 1][j]) {
            i--;
        }
    }

    reverse(ve.begin(), ve.end());

    return vt[vl.size()][vr.size()];
}

string LCS(string S, string T) {
    vector<char> vl, vr;
    for (char t : S)
        vl.push_back(t);
    for (char t : T)
        vr.push_back(t);

    vector<char> ve;
    LCS_base<char>(vl, vr, ve);
    string ret = "";
    for (auto t : ve)
        ret += t;

    return ret;
}

ull power(ull a, ull b, ull MOD = mod) {
    ull ans = 1;
    while (b) {
        if (b % 2) {
            ans = ans * a;
            if (MOD)
                ans = ans % MOD;
        }
        a = a * a;
        if (MOD) {
            a = a % MOD;
        }
        b /= 2;
    }
    return ans;
}
vector<pair<char, ll>> lanlen(string S) {
    char bef = S[0];
    bef--;
    vector<pair<char, ll>> ret;
    rep(i, 0, (ll)S.size()) {
        if (S[i] != bef) {
            ret.push_back({S[i], 1});
        } else {
            ret[ret.size() - 1].second++;
        }
        bef = S[i];
    }
    return ret;
}
// グラフを与える
//
vector<ll> topological_sort(vector<vector<ll>> &G) {
    ll N = G.size();
    vector<ll> deg(N, 0);
    rep(i, 0, N) {
        for (ll t : G[i])
            deg[t]++;
    }
    queue<ll> q;
    rep(i, 0, N) {
        if (deg[i] == 0)
            q.push(i);
    }
    vector<ll> ret;
    while (q.size()) {
        ll now = q.front();
        q.pop();
        ret.push_back(now);
        for (ll t : G[now]) {
            deg[t]--;
            if (deg[t] == 0)
                q.push(t);
        }
    }
    return ret;
}

template <class T = long long> struct UnionFindTemplate {
    vector<T> par; // par[i]:iの親の番号　(例) par[3] = 2 : 3の親が2
    vector<T> vsize;
    T Gc;

    UnionFindTemplate() {}

    UnionFindTemplate(T N)
        : par(N), vsize(N) { // 最初は全てが根であるとして初期化
        for (T i = 0; i < N; i++)
            par[i] = i;
        for (T i = 0; i < N; i++)
            vsize[i] = 1;
        Gc = N;
    }

    T root(T x) { // データxが属する木の根を再帰で得る：root(x) = {xの木の根}
        if (par[x] == x)
            return x;
        return par[x] = root(par[x]);
    }

    void unite(T x, T y) { // xとyの木を併合
        T rx = root(x);    // xの根をrx
        T ry = root(y);    // yの根をry
        if (rx == ry)
            return; // xとyの根が同じ(=同じ木にある)時はそのまま
        Gc--;
        vsize[ry] += vsize[rx];
        par[rx] =
            ry; // xとyの根が同じでない(=同じ木にない)時：xの根rxをyの根ryにつける
    }

    bool same(T x, T y) { // 2つのデータx, yが属する木が同じならtrueを返す
        T rx = root(x);
        T ry = root(y);
        return rx == ry;
    }
    T size(T x) { return vsize[root(x)]; }
};
using UnionFind = UnionFindTemplate<ll>;

vector<ll> era(ll n) {
    vector<ll> ret(n, 0);
    ret[1] = 1;
    rep(i, 1, n) {
        if (!ret[i]) {
            for (int j = 1; i * j <= n; j++) {
                ret[i * j] = i;
            }
        }
    }
    return ret;
}
vector<vector<ll>> extera(ll n) {
    vector<ll> v = era(n);
    vector<vector<ll>> ret(n);
    rep(i, 2, n) {
        ll now = i;
        while (now > 1) {
            ret[i].push_back(v[now]);
            now /= v[now];
        }
    }

    return ret;
}

// ヒストグラフ内最大長方形
long long histglaphMaxRect(vector<ll> v) {
    for (ll t : v)
        assert(t >= 0);
    ll ret = 0;
    v.push_back(-1);

    ll N = v.size();
    stack<tuple<ll, ll, ll>> s;
    s.push({-1, -1, -1});
    rep(i, 0, N) {
        ll lastleftindex = i;
        while (get<0>(s.top()) > v[i]) {
            auto [s1, s2, s3] = s.top();
            s3 = i;
            if (ret < s1 * (i - s2)) {
                ret = s1 * (i - s2);
            }
            lastleftindex = s2;
            s.pop();
        }
        s.push({v[i], lastleftindex, -1});
    }

    v.erase(--(v.end()));
    return ret;
};

template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S),
          F (*composition)(F, F), F (*id)()>
class mylazy_segtree {
  private:
    bool debugflg = false;
    ll n;
    ll size;
    ll log;
    vector<S> d;
    vector<F> ld;
    void update(ll k) {
        if (k < size) {
            d[k] = op(d[2 * k], d[2 * k + 1]);
        }
    }
    void all_apply(ll k, F f) {
        d[k] = mapping(f, d[k]);
        if (k < size)
            ld[k] = composition(f, ld[k]);
    }
    void push(ll k) {
        if (k < size) {
            all_apply(2 * k, ld[k]);
            all_apply(2 * k + 1, ld[k]);
            ld[k] = id();
        }
    }
    S prod(ll l, ll r, ll a, ll b, ll k) {
        assert(0 <= l && l <= r && r <= n);
        if (b <= l || r <= a)
            return e();
        push(k);
        if (l <= a && b <= r) {
            return d[k];
        }
        S t1 = prod(l, r, a, (a + b) / 2, 2 * k);
        S t2 = prod(l, r, (a + b) / 2, b, 2 * k + 1);
        return op(t1, t2);
    }
    void apply(ll l, ll r, F f, ll a, ll b, ll k) {
        if (b <= l || r <= a)
            return;
        push(k);
        if (l <= a && b <= r) {
            all_apply(k, f);
            // push(k);
        } else {
            apply(l, r, f, a, (a + b) / 2, 2 * k);
            apply(l, r, f, (a + b) / 2, b, 2 * k + 1);
            update(k);
        }
        return;
    }

  public:
    mylazy_segtree() {}
    mylazy_segtree(ll _n) : n(_n) {
        size = 1;
        log = 1;
        while (size < n) {
            size <<= 1LL;
            log++;
        }
        d.resize(2 * size, e());
        ld.resize(size, id());
    }
    mylazy_segtree(vector<S> v) : mylazy_segtree(v.size()) {
        rep(i, 0, n) set(i, v[i]);
        rrep(i, size - 1, 1) update(i);
    }

    void set(ll i, S p) {
        i += size;
        for (ll j = log; j >= 1; j--)
            push(i >> j);
        d[i] = p;
        for (ll j = 1; j <= log; j++)
            update(i >> j);
    }
    S get(ll i) {
        i += size;
        for (ll j = log; j >= 1; j--)
            push(i >> j);
        return d[i];
    }
    S prod(ll l, ll r) {
        S tmp = prod(l, r, 0, size, 1);
        return tmp;
    }
    S all_prod() { return d[1]; }
    void apply(ll l, ll r, F f) {
        // cout << "apply start" << endl;
        apply(l, r, f, 0, size, 1);
        // cout << "apply end" << endl;
    }
};

struct MyCombination {
  private:
  public:
    ll n;
    vector<ll> memokai;
    vector<ll> memokaiinv;

    MyCombination() : n(0) {}
    MyCombination(ll _n) : n(_n) {
        if (_n <= 0)
            return;
        memokai.resize(n + 1);
        memokai[1] = 1;
        rep(i, 2, n + 1) { memokai[i] = (memokai[i - 1] * i) % mod1e7; }
        memokaiinv.resize(n + 1);
        rep(i, 1, n + 1) {
            memokaiinv[i] = power(memokai[i], mod1e7 - 2LL, mod1e7);
        }
    }
    ll comb(ll a, ll b) {
        if (b == 0 || a == b)
            return 1;
        if (b == 1)
            return a;
        ll ret = 1;
        ret = (memokai[a] * memokaiinv[a - b]) % mod1e7;
        ret = (ret * memokaiinv[b]) % mod1e7;
        return ret;
    }
};

// #include <atcoder/all>
// using namespace atcoder;
/*
atcoder::lazy_segtree<S, op, e, F, mapping, composition, id> seg(N);
*/
// #include <atcoder/all>
// using mint = modint998244353;
vector<ll> ten(211000), one(211000);
using S = ll;
using F = ll;
const F ID = INF;
S op(S a, S b) { return a; }
S e() { return 0; }
S mapping(F f, S x) { return 0; }
F composition(F f, F g) { return 0; }
F id() { return INF; }

string AtoB(ll K, ll keta) {
    if (K == 0)
        return "0";
    string ans = "";
    while (K > 0) {
        ans = (char)('0' + K % keta) + ans;
        K /= keta;
    }
    return ans;
}
string AtoB(string K, ll keta) { return AtoB(atoll(&K[0]), keta); }

/*
 *   ABC121
 *   XORWorld
 */
ll XORsum(ll in) {
    assert(0 <= in);
    if (in % 2) {
        if ((in / 2) % 2 == 0)
            return 1;
        else
            return 0;
    } else {
        if ((in / 2) % 2 == 0)
            return in;
        else
            return in + 1;
    }
}

// priority_queue
using PQ = pair<ll, ll>;
// 昇順
// auto compare = [](PQ a, PQ b) {
//     // 比較を定義
//     return a > b;
// };
// 降順
// auto compare = [](PQ a, PQ b) {
//     // 比較を定義
//     return a < b;
// };
// auto compare = [](PQ a, PQ b) {
//     // 比較を定義
//     return (a.first == b.first) ? a.second > b.second : a.first > b.first;
// };
auto compare = [](PQ a, PQ b) {
    // 比較を定義
    return (a.first == b.first) ? a.second < b.second : a.first < b.first;
};

priority_queue<PQ, vector<PQ>,
               decltype(compare) // 比較関数オブジェクトを指定
               >
    pq{compare};

// 全方位木dp 実績
// { "AOJ GRL_5_A", "AOJ 1595 Traffic Tree", "EDPC V" , "限界集落", "ABC 222 F"}
struct edgeinfo {
    ll from;
    ll to;
    edgeinfo(ll par, ll in) {
        from = par;
        to = in;
    }
};

template <
    class reroot_S, reroot_S (*reroot_merge)(reroot_S, reroot_S, edgeinfo, ll),
    reroot_S (*reroot_sum_merge)(reroot_S, reroot_S), reroot_S (*reroot_id)()>
class reroot {
  private:
  public:
    ll N;
    vector<vector<pair<ll, ll>>> G;
    vector<reroot_S> dp;
    vector<vector<reroot_S>> leftsum;
    vector<vector<reroot_S>> rightsum;
    vector<map<ll, ll>> child_num;
    vector<reroot_S> vans;

    // カスタマイズ
    ll M;
    // カスタマイズ

    reroot(ll n) {
        N = n;
        G.resize(N);
        dp.resize(N, reroot_id()); // 単位元
        leftsum.resize(N);
        rightsum.resize(N);
        child_num.resize(N);
        vans.resize(N);
    }

    void set(ll a, ll b, ll w1 = 1, ll w2 = 1) {
        G[a].push_back({b, w1});
        G[b].push_back({a, w1});
    }

    reroot_S dfs1(ll in, ll par) {
        rep(i, 0, G[in].size()) {
            auto [to, w] = G[in][i];
            // ★1
            // child_num[in][to] = current_size++;
            if (to == par)
                continue;

            dp[in] = reroot_merge(dp[in], dfs1(to, in), edgeinfo(in, to), w);
        }
        return dp[in];
    }

    reroot_S dfs2(ll in, ll par) {
        leftsum[in].resize(G[in].size() + 1, reroot_id());
        rightsum[in].resize(G[in].size() + 1, reroot_id());

        rep(i, 0, G[in].size()) {
            auto [to, w] = G[in][i];
            leftsum[in][i + 1] =
                reroot_merge(leftsum[in][i], dp[to], edgeinfo(in, to), w);
        }
        rrep(i, G[in].size() - 1, 0) {
            auto [to, w] = G[in][i];
            rightsum[in][i] =
                reroot_merge(rightsum[in][i + 1], dp[to], edgeinfo(in, to), w);
        }

        reroot_S ans = reroot_id();
        for (auto [to, w] : G[in]) {
            if (to == par) {
                ans = reroot_merge(ans, dp[par], edgeinfo(in, to), w);
            } else {
                ans = reroot_merge(ans, dp[to], edgeinfo(in, to), w);
            }
        }
        vans[in] = ans;
        // ★1
        rep(i, 0, G[in].size()) {
            auto [to, w] = G[in][i];
            reroot_S ret = reroot_id();
            if (to == par) {
                ret = reroot_merge(ret, dp[par], edgeinfo(in, to), w);
            } else {
                // merge 2
                dp[in] = reroot_sum_merge(leftsum[in][i], rightsum[in][i + 1]);

                ret = dfs2(to, in);
                ret = reroot_merge(ret, dp[to], edgeinfo(in, to), w);
            }
        }
        return dp[in];
    }

    void build() {
        dfs1(0, -1);
        dfs2(0, -1);
    }

    reroot_S get(ll in) { return vans[in]; }
};

using reroot_S = ll;
reroot_S reroot_merge(reroot_S a, reroot_S b, edgeinfo etoe, ll w) {
    // ---
    return max(a, b + w);
}
reroot_S reroot_sum_merge(reroot_S a, reroot_S b) {
    // ------
    return max(a, b);
}
reroot_S reroot_id() {
    // ------
    return 0;
}
// reroot<reroot_S,reroot_merge,reroot_sum_merge,reroot_id>

class Main {
  private:
  public:
    void solve() {
        ll N;
        cin >> N;
        vector<ll> a(N);
        rep(i, 0, N) cin >> a[i];

        vector<vector<vector<ll>>> dp(
            N + 1, vector<vector<ll>>(N + 1, vector<ll>(N + 1, 0)));

        dp[0][0][0] = 1;
        rep(i, 0, N) {
            rep(j, 0, N) {
                rep(k, 0, N) {
                    if (dp[i][j][k] == 0)
                        continue;
                    dp[i + 1][j + 1][(k + a[i]) % (j + 1)] += dp[i][j][k];
                    dp[i + 1][j + 1][(k + a[i]) % (j + 1)] %= mod;
                    dp[i + 1][j][k] += dp[i][j][k];
                    dp[i + 1][j][k] %= mod;
                }
            }
        }

        ll ans = 0;
        rep(i, 1, N + 1) ans = (ans + dp[N][i][0]) % mod;

        cout << ans << endl;

        return;
    };
};

#include <random>
void test() {
    cout << 1000000000 << " " << 1000000000 << endl;
    ll N = 200000;
    cout << N << endl;
    ll p = 1;
    ll q = 1;
    rep(i, 0, N) {
        cout << p << " " << q << endl;
        p += 2;
        q += 2;
    }
    ll A = 200000;
    ll B = 200000;
    cout << A << endl;
    ll a = 0;
    ll b = 0;
    rep(i, 0, A) {
        cout << a << " ";
        a += 2;
    }
    cout << endl;
    cout << B << endl;
    rep(i, 0, B) {
        cout << b << " ";
        b += 2;
    }
    cout << endl;
}
void wraptest(ll n) {
    rep(i, 0, n) { test(); }
}

//
// テストケース作成 _DEBB1
// ランダムテストケース実施 _DEBB2
// その他 _EXEC
#define _EXEC 1

int submain() {

    Main m;
    ll t;
    t = 1;
#if defined _DEBB1 || defined _DEBB2
    cin >> t;
#endif

    rep(i, 0, t) {
#ifdef _DEBB2
        ll cn;
        cin >> cn;
        cout << cn << " : " << endl;
#endif
#ifndef _DEBB1
        m.solve();
#endif
    }

    {
#ifdef _DEBB1
        cout << t << endl;
        rep(i, 0, t) {
            cout << i << endl;
            test();
        }
#endif
    }

    return 0;
}

#ifndef INCDEF
int main() { submain(); };
#endif