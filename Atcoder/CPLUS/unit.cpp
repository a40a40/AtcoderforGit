#define INCDEF 1
#include "test.cpp"

// ------------------------------------------------------------------------------------
vector<ll> D_ABC222F;
using reroot_S_ABC222F = ll;
reroot_S_ABC222F reroot_merge_ABC222F(reroot_S_ABC222F a, reroot_S_ABC222F b,
                                      edgeinfo etoe, ll w) {
    // ---
    return max(max(a, D_ABC222F[etoe.to] + w), b + w);
}
reroot_S_ABC222F reroot_sum_merge_ABC222F(reroot_S_ABC222F a,
                                          reroot_S_ABC222F b) {
    // ------
    return max(a, b);
}
reroot_S_ABC222F reroot_id_ABC222F() {
    // ------
    return 0;
}

int solve_ABC222F() {
    ll N;
    cin >> N;

    D_ABC222F.resize(N);
    reroot<reroot_S_ABC222F, reroot_merge_ABC222F, reroot_sum_merge_ABC222F,
           reroot_id_ABC222F>
        rr1(N);
    rep(i, 0, N - 1) {
        ll a, b, c;
        cin >> a >> b >> c;
        a--;
        b--;
        rr1.set(a, b, c);
    }
    rep(i, 0, N) cin >> D_ABC222F[i];

    rr1.build();

    rep(i, 0, N) cout << rr1.get(i) << endl;
    return 0;
}
// ------------------------------------------------------------------------------------

int main() {

    cout << "Hello World" << endl;
    solve_ABC222F();

    return 0;
}
