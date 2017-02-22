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
const double Q = 1.0;
const double T = 1.0;
//------------------------------------------------------------------------
struct Vars {
  double p;
  double q;
  double zeta;
  double eta;
  Vars (double _p, double _q, double _zeta, double _eta) {
    p = _p;
    q = _q;
    zeta = _zeta;
    eta = _eta;
  }
  Vars () {
    p = 0.0;
    q = 0.0;
    zeta = 0.0;
    eta = 0.0;
  }
  Vars operator+(const Vars &v) {
    return Vars(p + v.p, q + v.q, zeta + v.zeta, eta + v.eta);
  }
  Vars operator*(const double f) {
    return Vars(p * f, q * f, zeta * f, eta * f);
  }
};
//------------------------------------------------------------------------
class NoseHoover {
public:
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta;
    dv.q = v.p;
    dv.zeta = (v.p * v.p - T) / Q;
    dv.eta = T * v.zeta;
    return dv;
  }
};
//------------------------------------------------------------------------
class KineticMoments {
public:
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta - v.p * v.p * v.p * v.eta;
    dv.q = v.p;
    dv.zeta = (v.p * v.p - T) / Q;
    dv.eta = (v.p * v.p * v.p * v.p - 3.0 * T * v.p * v.p) / Q;
    return dv;
  }
};
//------------------------------------------------------------------------
class NoseHooverChain {
public:
  Vars operator()(Vars v) {
    Vars dv;
    dv.p = -v.q - v.p * v.zeta;
    dv.q = v.p;
    dv.zeta = (v.p * v.p - T) / Q - v.eta * v.zeta;
    dv.eta = (v.zeta * v.zeta - T) / Q;
    return dv;
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
};
//------------------------------------------------------------------------
template <class Method>
void integrate(Method f, std::string name) {
  Vars v(0.0, 1.0, 0.0, 0.0);
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
      ofs_p << v.q << std::endl;
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
  integrate(RungeKutta<NoseHoover>(dt), "nose_hoover");
  integrate(RungeKutta<KineticMoments>(dt), "kinetic_moments");
  integrate(RungeKutta<NoseHooverChain>(dt), "nose_hoover_chain");
  integrate(Euler<Langevin>(dt), "langevin");
}
//------------------------------------------------------------------------
