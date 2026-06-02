#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

typedef long double ld;
typedef long long ll;
typedef pair<ll, ld> pif;

ld Pow(ld x, ll i) {
    if (i == 0) return 1;
    ld tmp = Pow(x,i/2);
    tmp *= tmp;

    if(i%2) tmp *= x;
    return tmp;
}

const ld zb = 0.000000001;

bool is_zero(const ld &c, const ld zero = zb) {
    return abs(c) < zero;
}
bool is_equal(const pif& a, const pif& b) {
    return is_zero(a.second-b.second) && (a.first == b.first);
}

class Poly {
public:
    Poly() : n(1), v({{0,0}}) {}
    explicit Poly(const ld& c) : n(1), v({{0,c}}) {}
    Poly(const ll& exp, const ld& c) : n(exp+1), v({{exp,c}}) {
        if(is_zero(c) && exp != 0) {
            n = 1;
            v[0] = {0,0};
        }
    }
    explicit Poly(const vector<ld>& V) : n(1) {
        for (int i = 0; i < (int)V.size(); ++i) {
            if (!is_zero(V[i])) {
                v.emplace_back(i, V[i]);
                n = i + 1;
            }
        }
        if (v.empty()) v.emplace_back(0,0);
        sort(v.begin(), v.end());
    }
    Poly(const vector<pif>& V) : n(1) {
        for (auto &ec : V) {
            if (ec.second != 0) {
                v.push_back(ec);
                n = ec.first + 1; // only sorted V
            }
        }
        if (v.empty()) v.emplace_back(0,0);
        sort(v.begin(), v.end());
    }
    Poly(const Poly& p) = default;
    ~Poly() = default;

    friend bool operator==(const Poly& p, const Poly& q);
    friend bool operator==(const Poly& p, const ld &c);
    friend Poly operator-(const Poly& p, const Poly& q);
    friend Poly operator+(const Poly& p, const Poly& q);
    friend Poly operator*(const Poly& p, const Poly& q);
    friend Poly operator*(const ld &c, const Poly& q);
    friend pair<Poly,Poly> operator/(const Poly& p, const Poly& q);

    friend ostream& operator<<(ostream& os, const Poly& p);

    ld eval(ld x) {
        ld ans = 0;
        for (auto &ec : v)
            ans += Pow(x,ec.first)*ec.second;

        return ans;
    }

private:
    ll n; // deg+1
    vector<pif> v; // c, exp
};

bool operator==(const Poly& p, const Poly& q) {
    if (p.v.size() != q.v.size()) return false;

    for (int i = 0; i < (int)p.v.size(); ++i)
        if (!is_equal(p.v[i], q.v[i]))
            return false;

    return true;
}
bool operator==(const Poly& p, const ld &c) {
    Poly C(c);
    return (p == C);
}
Poly operator-(const Poly& p, const Poly& q) {
    vector<pif> v;

    int i = 0, j = 0;
    while(i < p.v.size() || j < q.v.size()) {
        if (j >= q.v.size()) v.push_back(p.v[i++]);
        else if (i >= p.v.size()) {
            v.emplace_back(q.v[j].first, -q.v[j].second);
            ++j;
        } else {
            if (p.v[i].first == q.v[j].first) {
                ld val = p.v[i].second - q.v[j].second;
                if (!is_zero(val))
                    v.emplace_back(p.v[i].first, val);

                ++i;
                ++j;
            } else if (p.v[i].first < q.v[j].first) {
                v.push_back(p.v[i++]);
            } else {
                v.emplace_back(q.v[j].first, -q.v[j].second);
                ++j;
            }
        }
    }
    if(v.empty()) return {};
    return {v};
}
Poly operator+(const Poly& p, const Poly& q) {
    vector<pif> v;
    int i = 0, j = 0;
    while(i < p.v.size() || j < q.v.size()) {
        if (j >= q.v.size()) v.push_back(p.v[i++]);
        else if (i >= p.v.size()) v.push_back(q.v[j++]);
        else {
            if (p.v[i].first == q.v[j].first) {
                ld val = p.v[i].second+q.v[j].second;
                if (!is_zero(val))
                    v.emplace_back(p.v[i].first, val);

                ++i;
                ++j;
            } else if (p.v[i].first < q.v[j].first)
                v.push_back(p.v[i++]);
            else
                v.push_back(q.v[j++]);
        }
    }

    if(v.empty()) return {};
    return {v};
}

Poly operator*(const Poly& p, const Poly& q) {
    vector<pif> raw;

    raw.reserve(p.v.size() * q.v.size());

    for (const auto& ecp : p.v)
        for (const auto& ecq : q.v)
            raw.emplace_back(ecp.first + ecq.first, ecp.second * ecq.second);

    sort(raw.begin(), raw.end());

    vector<pif> merged;
    for (auto& ec : raw) {
        if (!merged.empty() && merged.back().first == ec.first) {
            merged.back().second += ec.second;
        } else {
            if (!merged.empty() && is_zero(merged.back().second)) {
                merged.pop_back();
            }

            merged.push_back(ec);
        }
    }
    if (!merged.empty() && is_zero(merged.back().second)) {
        merged.pop_back();
    }

    if (merged.empty()) merged.emplace_back(0, 0);

    return {merged};
}
Poly operator*(const ld &c, const Poly& q) {
    Poly C(c);
    return C*q;
}
pair<Poly,Poly> operator/(const Poly& p, const Poly& q) {
    if(p == 0) return {Poly(0),Poly(0)};
    if(q == 0) exit(0);
    if (q.n > p.n) return {Poly(0), p};

    Poly quot(0);
    Poly rem = p;

    while (rem.n >= q.n && !(rem == 0)) {
        ld c = rem.v.back().second / q.v.back().second;
        ll exp = rem.n - q.n;

        Poly mul(exp, c);


        if (mul == 0) {
            rem.v.pop_back();
            if (rem.v.empty()) {
                rem = Poly(0);
            } else {
                rem.n = rem.v.back().first + 1;
            }
            continue;
        }

        quot = quot + mul;
        rem = rem - q * mul;
    }

    return {quot, rem};
}
ostream& operator<<(ostream& os, const Poly& p) {
    bool b = false;
    for(int i = (int)p.v.size() - 1; i >= 0; --i) {
        if(is_zero(p.v[i].second, 0.0000001)) continue;
        b = true;
        if (i == (int)p.v.size() - 1 && p.v[i].second < 0) os << '-';
        else if (i != (int)p.v.size() - 1) os << (p.v[i].second >= 0 ? '+' : '-') << ' ';

        if (abs(p.v[i].second) != 1 || i == 0) os << abs(p.v[i].second);
        if (p.v[i].first > 0) os << 'x';
        if (p.v[i].first > 1) os << '^' << p.v[i].first;
        os << ' ';
    }
    if(!b) os << '0';

    return os;
}

int main() {
    Poly p({{4,21.11},{10,50.3}});
    Poly q({1,2,3.5});
    cout << "p(x) = " << p << '\n';
    cout << "q(x) = " << q << '\n';
    auto div = p/q;
    cout << "p(x) = q(x)*m(x) + r(x)" << '\n';
    cout << "m(x) = " << div.first << '\n';
    cout << "r(x) = " << div.second << '\n';

    cout << '\n';
    cout << q*div.first + div.second << "\n\n";

    ld x = 0.9;
    cout << "p(" << x << ") = " << p.eval(x) << '\n';
    cout << "q(" << x << ")*m(" << x << ") + r(" << x << ") = " <<
        q.eval(x) * div.first.eval(x) + div.second.eval(x) << '\n';

    return 0;
}
