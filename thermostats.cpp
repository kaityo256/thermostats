//------------------------------------------------------------------------
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
//------------------------------------------------------------------------
const int LOOP = 10000000;
const int O_LOOP = 1000;
const double dt = 0.001;
const double T = 1.0;
//------------------------------------------------------------------------
struct Vars {
  double p;
  double q;
  double zeta;
  double eta;
  double theta;
  Vars (double _p, double _q, double _zeta, double _eta, double _theta) {
    p = _p;
    q = _q;
    zeta = _zeta;
    eta = _eta;
    theta = _theta;
  }
  Vars () {
    p = 0.0;
    q = 0.0;
    zeta = 0.0;
    eta = 0.0;
    theta = 0.0;
  }
  Vars operator+(const Vars &v) {
    return Vars(p + v.p, q + v.q, zeta + v.zeta, eta + v.eta, theta + v.theta);
  }
  Vars operator*(const double f) {
    return Vars(p * f, q * f, zeta * f, eta * f, theta * f);
  }
};
//------------------------------------------------------------------------
class NoseHoover {
private:
  double Q;
public:
  NoseHoover() {
    Q = 4.0;
  }
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta / Q;
    dv.q = v.p;
    dv.zeta = v.p * v.p - T;
    dv.eta = v.zeta * T / Q;
    return dv;
  }
  double H(Vars v) {
    double e = 0;
    e += (v.p * v.p * 0.5);
    e += (v.q * v.q * 0.5);
    e += (v.zeta * v.zeta * 0.5 / Q);
    e += v.eta;
    return e;
  }
};
//------------------------------------------------------------------------
class KineticMoments {
private:
  double Qzeta, Qeta;
public:
  KineticMoments() {
    Qzeta = 4.0;
    Qeta = 6.0;
  }
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta / Qzeta - v.p * v.p * v.p * v.eta / Qeta;
    dv.q = v.p;
    dv.zeta = v.p * v.p - T;
    dv.eta = v.p * v.p * v.p * v.p - 3.0 * T * v.p * v.p;
    dv.theta = v.zeta * T / Qzeta + 3.0 * T * v.p * v.p * v.eta / Qeta;
    return dv;
  }
  double H(Vars v) {
    double e = 0;
    e += (v.p * v.p * 0.5);
    e += (v.q * v.q * 0.5);
    e += (v.zeta * v.zeta * 0.5 / Qzeta);
    e += (v.eta * v.eta * 0.5 / Qeta);
    e += v.theta;
    return e;
  }
};
//------------------------------------------------------------------------
class NoseHooverChain {
private:
  double Qzeta;
  double Qeta;
public:
  NoseHooverChain() {
    Qzeta = 2.0;
    Qeta = 5.0;
  }
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta / Qzeta;
    dv.q = v.p;
    dv.zeta = v.p * v.p - T - v.eta * v.zeta / Qeta;
    dv.eta = v.zeta * v.zeta / Qzeta - T;
    dv.theta = T * (v.zeta / Qzeta + v.eta / Qeta);
    return dv;
  }
  double H(Vars v) {
    double e = 0;
    e += (v.p * v.p * 0.5);
    e += (v.q * v.q * 0.5);
    e += (v.zeta * v.zeta * 0.5 / Qzeta);
    e += (v.eta * v.eta * 0.5 / Qeta);
    e += v.theta;
    return e;
  }
};

//------------------------------------------------------------------------
class Langevin {
public:
  std::mt19937 mt;
  const double gamma;
  const double D;
  std::normal_distribution<double> nd;
  Langevin() : mt(1), gamma(1.0), D(sqrt(2.0 * gamma * T / dt)), nd(0.0, D) {
  }
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - gamma * v.p + nd(mt);
    dv.q = v.p;
    return dv;
  }
  double H(Vars ) {
    return 0.0;
  }
};
//------------------------------------------------------------------------
template <class DIFF>
class RungeKutta {
private:
  const double dt;
  const double hdt;
  DIFF diff;
public:
  RungeKutta(double _dt) : dt(_dt), hdt(_dt * 0.5) {};
  Vars operator()(Vars v) {
    Vars k1 = diff(v);
    Vars k2 = diff(v + k1 * hdt);
    Vars k3 = diff(v + k2 * hdt);
    Vars k4 = diff(v + k3 * dt);
    v = v + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
    return v;
  }
  double H(Vars v) {
    return diff.H(v);
  }
};
//------------------------------------------------------------------------
template <class DIFF>
class Euler {
private:
  const double dt;
  DIFF diff;
public:
  Euler(double _dt) : dt(_dt) {};
  Vars operator()(Vars v) {
    v = v + diff(v) * dt;
    return v;
  }
  double H(Vars v) {
    return diff.H(v);
  }
};
//------------------------------------------------------------------------
template <class Method>
void integrate(Method f, std::string name) {
  Vars v;
  v.q = 1.0;
  double t = 0;
  std::string pfile = name + "_ps.dat";
  std::ofstream ofs_p(pfile.c_str());
  std::cout << name << std::endl;
  std::vector<double> data_e;
  for (int i = 0; i < LOOP; i++) {
    t += dt;
    v = f(v);
    if (i % O_LOOP == 0) {
      ofs_p << t << " ";
      ofs_p << v.p << " ";
      ofs_p << v.q << " ";
      ofs_p << f.H(v) << std::endl;
      data_e.push_back(v.p * v.p * 0.5 + v.q * v.q * 0.5);
    }
  }
  std::sort(data_e.begin(), data_e.end());
  const int n = data_e.size();
  const double nd = static_cast<double>(n);
  std::string efile = name + ".dat";
  std::ofstream ofs(efile.c_str());
  for (int i = 0; i < n; i++) {
    ofs << data_e[i] << " " << static_cast<double>(i) / nd << std::endl;
  }
}
//------------------------------------------------------------------------
int
main(void) {
  integrate(RungeKutta<NoseHoover>(dt), "nh_rk");
  integrate(Euler<NoseHoover>(dt), "nh_euler");

  integrate(RungeKutta<KineticMoments>(dt), "km_rk");
  integrate(Euler<KineticMoments>(dt), "km_euler");

  integrate(RungeKutta<NoseHooverChain>(dt), "nhc_rk");
  integrate(Euler<NoseHooverChain>(dt), "nhc_euler");

  integrate(Euler<Langevin>(dt), "langevin");
}
//------------------------------------------------------------------------
