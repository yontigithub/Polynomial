//
// Created by Yonatan Rappoport on 02/06/2026.
//

#ifndef SPARSE_POLY_DENSEPOLY_H
#define SPARSE_POLY_DENSEPOLY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "Utils.h"
using namespace std;

class DensePoly  {
public:
    DensePoly () : n(1), v({0}) {}
    DensePoly (const ld& c) : n(1), v({c}) {}
    DensePoly (const ld& c, const int& exp) : n(exp+1) {
        v.resize(exp+1);
        for(int i = 0; i < exp; ++i) v[i] = 0;
        v[exp] = c;
    }
    DensePoly (const vector<ld>& V) : n((int)V.size()), v(V) {
        if(V.empty()) {
            ++n;
            v.push_back(0);
        }
    }
    DensePoly (const DensePoly & p) = default;
    ~DensePoly () = default;

    friend bool operator==(DensePoly & p, DensePoly & q);
    friend bool operator==(DensePoly & p, const ld &c);
    friend DensePoly  operator-(const DensePoly & p, const DensePoly & q);
    friend DensePoly  operator+(const DensePoly & p, const DensePoly & q);
    friend DensePoly  operator*(const DensePoly & p, const DensePoly & q);
    friend DensePoly  operator*(const ld &c, const DensePoly & q);
    friend ostream& operator<<(ostream& os, DensePoly & p);

    ld eval(ld x) {
        ld ans = 0;
        for (int i = 0; i < n; ++i) {
            ans += Pow(x,i)*v[i];
        }
        return ans;
    }
    void simplify() {
        int t = 1;
        for (int i = n - 1; i >= 0; --i) {
            if (!is_zero(v[i])) {
                t=i+1;
                break;
            }
        }
        if (t == n) return;
        v.resize(t);
        n = t;
    }
    bool is_const() {
        simplify();
        return (n == 1);
    }
    int deg() const {
        return n - 1;
    }

    pair<DensePoly ,DensePoly > div(DensePoly & p, DensePoly & q) { // p/q

        if(p == 0) return {DensePoly (0),DensePoly (0)}; // simplifies
        if(q == 0) return {NULL,NULL};// simplifies
        if (q.n > p.n) return {DensePoly (0), p};

        ld a = p.v[p.n - 1]/q.v[q.n - 1];
        int exp = p.n - q.n;

        DensePoly  mul(a,exp);
        DensePoly  sub = q*mul;
        DensePoly  new_p = p - sub;
        auto qr = div(new_p, q);
        qr.first = qr.first + mul;
        return qr;
    }

    int n; // deg+1
    vector<ld> v;
};

bool operator==(DensePoly & p, DensePoly & q) {
    p.simplify(); q.simplify();
    if (p.n != q.n) return false;

    for (int i = 0; i < p.n; ++i)
        if (p.v[i] != q.v[i])
            return false;

    return true;
}
bool operator==(DensePoly & p, const ld &c) {
    DensePoly  C(c);
    return (p == C);
}
DensePoly  operator-(const DensePoly & p, const DensePoly & q) {
    vector<ld> v(max(p.n, q.n) );
    for (int i = 0; i < max(p.n, q.n); ++i) {
        if (i < min(p.n,q.n)) v[i] = p.v[i] - q.v[i];
        else if (i < p.n) {v[i] = p.v[i];}
        else v[i] = -q.v[i];
    }
    return {v};
}
DensePoly  operator+(const DensePoly & p, const DensePoly & q) {
    vector<ld> v(max(p.n, q.n));
    for (int i = 0; i < max(p.n, q.n); ++i) {
        if (i < min(p.n,q.n)) v[i] = p.v[i] + q.v[i];
        else if (i < p.n) {v[i] = p.v[i];}
        else v[i] = q.v[i];
    }
    return {v};
}

DensePoly  operator*(const DensePoly & p, const DensePoly & q) {
    vector<ld> v(q.n*p.n);

    for (int i = 0; i < q.n*p.n; ++i) {
        v[i] = 0;
        for (int j = 0; j <= min(p.n-1,i); ++j) {
            if (i - j >= q.n) continue;
            v[i] += p.v[j]*q.v[i - j];
        }
    }
    return {v};
}
DensePoly  operator*(const ld &c, const DensePoly & q) {
    DensePoly  C(c);
    return C*q;
}
ostream& operator<<(ostream& os, DensePoly & p) {
    p.simplify();

    for(int i = p.n - 1; i > 0; --i) {
        if (is_zero(p.v[i])) continue;
        if (i == p.n - 1 && p.v[p.n - 1] < 0) os << '-';
        else if (i != p.n - 1) os << (p.v[i] >= 0 ? '+' : '-') << ' ';

        if (abs(p.v[i]) != 1) os << abs(p.v[i]);
        os << 'x';
        if (i > 1) os << '^' << i;
        os << ' ';
    }

    if(!is_zero(p.v[0]) || p.n == 1) os << (p.v[0] >= 0 ? '+' : '-') << ' ' << abs(p.v[0]);
    return os;
}


#endif //SPARSE_POLY_DENSEPOLY_H
