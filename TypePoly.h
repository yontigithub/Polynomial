//
// Created by Yonatan Rappoport on 02/06/2026.
//

#ifndef SPARSE_POLY_TYPEPOLY_H
#define SPARSE_POLY_TYPEPOLY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "Utils.h"
using namespace std;

template <typename T>
class TypePoly {
public:
    TypePoly();
    explicit TypePoly(const T& c);
    TypePoly(const ll& exp, const T& c);
    explicit TypePoly(const vector<T>& V);
    TypePoly(vector<pair<ll,T>> V);
    TypePoly(const TypePoly<T>& p) = default;
    ~TypePoly() = default;

    template <typename U>
    friend bool operator==(const TypePoly<U>& p, const TypePoly<U>& q);
    template <typename U>
    friend bool operator==(const TypePoly<U>& p, const U &c);
    template <typename U>
    friend TypePoly<U> operator-(const TypePoly<U>& p, const TypePoly<U>& q);
    template <typename U>
    friend TypePoly<U> operator+(const TypePoly<U>& p, const TypePoly<U>& q);
    template <typename U>
    friend TypePoly<U> operator*(const TypePoly<U>& p, const TypePoly<U>& q);
    template <typename U>
    friend TypePoly<U> operator*(const U &c, const TypePoly<U>& q);
    template <typename U>
    friend pair<TypePoly<U>,TypePoly<U>> operator/(const TypePoly<U>& p, const TypePoly<U>& q);
    template <typename U>
    friend ostream& operator<<(ostream& os, const TypePoly<U>& p);

    T eval(T x) const;

private:
    ll n; // deg+1
    vector<pair<ll, T>> v; // c, exp
};

template <typename T>
TypePoly<T>::TypePoly() : TypePoly<T>(T(0)) {}
template <typename T>
TypePoly<T>::TypePoly(const T& c) : n(1), v({{0,c}}) {}
template <typename T>
TypePoly<T>::TypePoly(const ll& exp, const T& c) : n(exp+1), v({{exp,c}}) {
    if(is_zero(c) && exp != 0) {
        n = 1;
        v[0] = {0,T(0)};
    }
}
template <typename T>
TypePoly<T>::TypePoly(const vector<T>& V) : n(1) {
    for (int i = 0; i < (int)V.size(); ++i) {
        if (!is_zero(V[i])) {
            v.emplace_back(i, V[i]);
            n = i + 1;
        }
    }
    if (v.empty()) v.emplace_back(0,T(0));
    sort(v.begin(), v.end());
}
template <typename T>
TypePoly<T>::TypePoly(vector<pair<ll, T>> V) : n(1) {
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
    if (v.empty()) v.emplace_back(0,T(0));

}

template <typename T>
bool operator==(const TypePoly<T>& p, const TypePoly<T>& q) {
    if (p.v.size() != q.v.size()) return false;

    for (int i = 0; i < (int)p.v.size(); ++i)
        if (!is_equal(p.v[i], q.v[i]))
            return false;

    return true;
}
template <typename T>
bool operator==(const TypePoly<T>& p, const T &c) {
    TypePoly<T> C(c);
    return (p == C);
}

template <typename T>
TypePoly<T> operator-(const TypePoly<T>& p, const TypePoly<T>& q) {
    vector<pair<ll, T>> v;

    int i = 0, j = 0;
    while(i < p.v.size() || j < q.v.size()) {
        if (j >= q.v.size()) v.push_back(p.v[i++]);
        else if (i >= p.v.size()) {
            v.emplace_back(q.v[j].first, -q.v[j].second);
            ++j;
        } else {
            if (p.v[i].first == q.v[j].first) {
                T val = p.v[i].second - q.v[j].second;
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
template <typename T>
TypePoly<T> operator+(const TypePoly<T>& p, const TypePoly<T>& q) {
    vector<pair<ll, T>> v;
    int i = 0, j = 0;
    while(i < p.v.size() || j < q.v.size()) {
        if (j >= q.v.size()) v.push_back(p.v[i++]);
        else if (i >= p.v.size()) v.push_back(q.v[j++]);
        else {
            if (p.v[i].first == q.v[j].first) {
                T val = p.v[i].second+q.v[j].second;
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

template <typename T>
TypePoly<T> operator*(const TypePoly<T>& p, const TypePoly<T>& q) {
    vector<pair<ll, T>> raw;

    raw.reserve(p.v.size() * q.v.size());

    for (const auto& ecp : p.v)
        for (const auto& ecq : q.v)
            raw.emplace_back(ecp.first + ecq.first, ecp.second * ecq.second);

    sort(raw.begin(), raw.end());

    vector<pair<ll, T>> merged;
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

    if (merged.empty()) merged.emplace_back(0, T(0));

    return {merged};
}
template <typename T>
TypePoly<T> operator*(const T &c, const TypePoly<T>& q) {
    TypePoly<T> C(c);
    return C*q;
}

template <typename T>
pair<TypePoly<T>,TypePoly<T>> operator/(const TypePoly<T>& p, const TypePoly<T>& q) {
    if(p == T(0)) return {TypePoly<T>(T(0)),TypePoly<T>(T(0))};
    if(q == T(0)) exit(0);
    if (q.n > p.n) return {TypePoly<T>(T(0)), p};

    TypePoly<T> quot(T(0));
    TypePoly<T> rem = p;

    while (rem.n >= q.n && !(rem == T(0))) {
        T c = rem.v.back().second / q.v.back().second;
        ll exp = rem.n - q.n;

        TypePoly<T> mul(exp, c);


        if (mul == T(0)) {
            rem.v.pop_back();
            if (rem.v.empty()) {
                rem = TypePoly<T>(T(0));
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

template <typename T>
ostream& operator<<(ostream& os, const TypePoly<T>& p) {
    bool b = false;
    for(int i = (int)p.v.size() - 1; i >= 0; --i) {
        if(is_zero(p.v[i].second, 0.0000001)) continue;
        b = true;
        if (i == (int)p.v.size() - 1 && p.v[i].second < T(0)) os << '-';
        else if (i != (int)p.v.size() - 1) os << (p.v[i].second >= T(0) ? '+' : '-') << ' ';

        if (abs(p.v[i].second) != 1 || i == 0) os << abs(p.v[i].second);
        if (p.v[i].first > 0) os << 'x';
        if (p.v[i].first > 1) os << '^' << p.v[i].first;
        os << ' ';
    }
    if(!b) os << '0';

    return os;
}

template <typename T>
T TypePoly<T>::eval(T x) const {
    T ans = 0;
    for (auto &ec : v)
        ans += Pow(x,ec.first)*ec.second;
    return ans;
}

#endif //SPARSE_POLY_TYPEPOLY_H
