#include <bits/stdc++.h>

using namespace std;
#define ll long long
#define ull unsigned long long
#define rep(repindex, s, n)                                                    \
    for (ll(repindex) = (s); (repindex) < n; (repindex)++)
#define rrep(repindex, s, n)                                                   \
    for (ll(repindex) = (s); (repindex) >= n; (repindex)--)
#define bit(n, k) ((n >> k) & 1) /*nのk bit目*/

bool lrdef(ll x, ll l, ll r) { return (x >= l && x < r); }

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

using reroot_S = long long;
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
    void setM(ll m) { M = m; };

    reroot_S reroot_merge(reroot_S a, reroot_S b, ll w) {
        return max(a, b + 1);
    }
    reroot_S reroot_sum_merge(reroot_S a, reroot_S b) { return max(a, b); }
    reroot_S reroot_id() { return 0; }
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

            dp[in] = reroot_merge(dp[in], dfs1(to, in), w);
        }
        return dp[in];
    }

    reroot_S dfs2(ll in, ll par) {
        leftsum[in].resize(G[in].size() + 1, reroot_id());
        rightsum[in].resize(G[in].size() + 1, reroot_id());

        rep(i, 0, G[in].size()) {
            auto [to, w] = G[in][i];
            leftsum[in][i + 1] = reroot_merge(leftsum[in][i], dp[to], w);
        }
        rrep(i, G[in].size() - 1, 0) {
            auto [to, w] = G[in][i];
            rightsum[in][i] = reroot_merge(rightsum[in][i + 1], dp[to], w);
        }

        reroot_S ans = reroot_id();
        for (auto [to, w] : G[in]) {
            if (to == par) {
                ans = reroot_merge(ans, dp[par], w);
            } else {
                ans = reroot_merge(ans, dp[to], w);
            }
        }
        vans[in] = ans;
        // ★1
        rep(i, 0, G[in].size()) {
            auto [to, w] = G[in][i];
            reroot_S ret = reroot_id();
            if (to == par) {
                ret = reroot_merge(ret, dp[par], w);
            } else {
                // merge 2
                dp[in] = reroot_sum_merge(leftsum[in][i], rightsum[in][i + 1]);

                ret = dfs2(to, in);
                ret = reroot_merge(ret, dp[to], w);
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

int main(void) {
    ll N;
    cin >> N;
    reroot rr(N);
    rep(i, 0, N - 1) {
        ll s, t, w;
        // cin >> s >> t >> w;
        cin >> s >> t;
        w = 1;
        s--;
        t--;
        rr.set(s, t);
    }

    rr.build();
    ll ans = 0;

    rep(i, 0, N) cout << 2 * (N - 1) - rr.get(i) << endl;

    return 0;
}
