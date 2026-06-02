//
// Created by Yonatan Rappoport on 02/06/2026.
//

#ifndef SPARSE_POLY_SPARSEPOLY_H
#define SPARSE_POLY_SPARSEPOLY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "Utils.h"
using namespace std;


class SparsePoly {
public:
    SparsePoly();
    explicit SparsePoly(const ld& c);
    SparsePoly(const ll& exp, const ld& c);
    explicit SparsePoly(const vector<ld>& V);
    SparsePoly(vector<pif> V);
    SparsePoly(const SparsePoly& p) = default;
    ~SparsePoly() = default;

    friend bool operator==(const SparsePoly& p, const SparsePoly& q);
    friend bool operator==(const SparsePoly& p, const ld &c);
    friend SparsePoly operator-(const SparsePoly& p, const SparsePoly& q);
    friend SparsePoly operator+(const SparsePoly& p, const SparsePoly& q);
    friend SparsePoly operator*(const SparsePoly& p, const SparsePoly& q);
    friend SparsePoly operator*(const ld &c, const SparsePoly& q);
    friend pair<SparsePoly,SparsePoly> operator/(const SparsePoly& p, const SparsePoly& q);

    friend ostream& operator<<(ostream& os, const SparsePoly& p);

    ld eval(ld x) const;

private:
    ll n; // deg+1
    vector<pif> v; // c, exp
};

SparsePoly::SparsePoly() : SparsePoly(0) {}
SparsePoly::SparsePoly(const ld& c) : n(1), v({{0,c}}) {}
SparsePoly::SparsePoly(const ll& exp, const ld& c) : n(exp+1), v({{exp,c}}) {
    if(is_zero(c) && exp != 0) {
        n = 1;
        v[0] = {0,0};
    }
}
SparsePoly::SparsePoly(const vector<ld>& V) : n(1) {
    for (int i = 0; i < (int)V.size(); ++i) {
        if (!is_zero(V[i])) {
            v.emplace_back(i, V[i]);
            n = i + 1;
        }
    }
    if (v.empty()) v.emplace_back(0,0);
    sort(v.begin(), v.end());
}
SparsePoly::SparsePoly(vector<pif> V) : n(1) {
    sort(V.begin(), V.end());
    for (auto &ec : V) {
        if (ec.second != 0) {
            if(!v.empty() && v.back().first == ec.first) {
                v.back().second += ec.second;
            } else {
                v.push_back(ec);
                n = ec.first + 1;
            }
        }
    }
    if (v.empty()) v.emplace_back(0,0);

}

bool operator==(const SparsePoly& p, const SparsePoly& q) {
    if (p.v.size() != q.v.size()) return false;

    for (int i = 0; i < (int)p.v.size(); ++i)
        if (!is_equal(p.v[i], q.v[i]))
            return false;

    return true;
}
bool operator==(const SparsePoly& p, const ld &c) {
    SparsePoly C(c);
    return (p == C);
}

SparsePoly operator-(const SparsePoly& p, const SparsePoly& q) {
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
SparsePoly operator+(const SparsePoly& p, const SparsePoly& q) {
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

SparsePoly operator*(const SparsePoly& p, const SparsePoly& q) {
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
SparsePoly operator*(const ld &c, const SparsePoly& q) {
    SparsePoly C(c);
    return C*q;
}

pair<SparsePoly,SparsePoly> operator/(const SparsePoly& p, const SparsePoly& q) {
    if(p == 0) return {SparsePoly(0),SparsePoly(0)};
    if(q == 0) exit(0);
    if (q.n > p.n) return {SparsePoly(0), p};

    SparsePoly quot(0);
    SparsePoly rem = p;

    while (rem.n >= q.n && !(rem == 0)) {
        ld c = rem.v.back().second / q.v.back().second;
        ll exp = rem.n - q.n;

        SparsePoly mul(exp, c);


        if (mul == 0) {
            rem.v.pop_back();
            if (rem.v.empty()) {
                rem = SparsePoly(0);
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

ostream& operator<<(ostream& os, const SparsePoly& p) {
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

ld SparsePoly::eval(ld x) const {
    ld ans = 0;
    for (auto &ec : v)
        ans += Pow(x,ec.first)*ec.second;
    return ans;
}

#endif //SPARSE_POLY_SPARSEPOLY_H
