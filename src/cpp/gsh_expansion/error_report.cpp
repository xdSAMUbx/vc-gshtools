// error_report.cpp: implementation of error_report.hpp (consolidated numerical, omission and commission errors).
#include "error_report.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace gsh
{
  namespace
  {
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double kTolerance = 1e-12;
    constexpr double kDigitCap = 16.0;

    // Fills the E_rel part of the numerical error from its maximum value.
    void fill_relative(NumericalError &ne, double e)
    {
      ne.available = true;
      ne.max_relative_error = e;
      if (std::isnan(e))
      {
        ne.within_1e12 = false;
        ne.digits = kNaN;
        return;
      }
      ne.within_1e12 = e <= kTolerance;
      ne.digits = e <= 0.0 ? kDigitCap : std::min(kDigitCap, -std::log10(e));
    }

    // Number formatting for the text report: "n/a" when not finite.
    std::string fmt(double v, int precision, bool scientific)
    {
      if (!std::isfinite(v))
        return "n/a";
      std::ostringstream os;
      os << (scientific ? std::scientific : std::fixed) << std::setprecision(precision) << v;
      return os.str();
    }

    // Number formatting for JSON: null when not finite.
    std::string json_number(double v)
    {
      if (!std::isfinite(v))
        return "null";
      std::ostringstream os;
      os << std::setprecision(17) << v;
      return os.str();
    }

    std::string json_bool(bool available, bool value)
    {
      if (!available)
        return "null";
      return value ? "true" : "false";
    }

    std::string json_escape(const std::string &s)
    {
      std::string out;
      out.reserve(s.size());
      for (char ch : s)
      {
        const unsigned char c = static_cast<unsigned char>(ch);
        switch (c)
        {
        case '"':
          out += "\\\"";
          break;
        case '\\':
          out += "\\\\";
          break;
        case '\n':
          out += "\\n";
          break;
        case '\r':
          out += "\\r";
          break;
        case '\t':
          out += "\\t";
          break;
        default:
          if (c < 0x20)
          {
            char buf[8];
            std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned>(c));
            out += buf;
          }
          else
            out += ch;
        }
      }
      return out;
    }
  }

  ErrorReport make_error_report(const Model &sigma_model, int nmax, double r, double gamma, const std::vector<std::pair<double, double>> &calc_ref,
                                const Points *pts, const ReportOptions &opt)
  {
    if (nmax < 2)
      throw std::invalid_argument("make_error_report: nmax must be >= 2");
    if (sigma_model.nmax < 2)
      throw std::invalid_argument("make_error_report: the sigma model needs degree >= 2");
    if (opt.n_tail <= nmax)
      throw std::invalid_argument("make_error_report: n_tail must exceed nmax");

    ErrorReport rep;
    rep.model_name = sigma_model.name;
    rep.GM = sigma_model.GM;
    rep.R = sigma_model.R;
    rep.model_nmax = sigma_model.nmax;
    rep.nmax = nmax;
    rep.r = r;
    rep.gamma = gamma;

    // Numerical error: the digit-loss bound depends only on the degree and colatitude.
    NumericalError &ne = rep.numerical;
    ne.colatitude = opt.colatitude;
    ne.digit_loss_bound = digit_loss_bound(nmax, opt.colatitude);
    if (!ne.available)
    {
      ne.max_relative_error = kNaN;
      ne.digits = kNaN;
    }
    if (!calc_ref.empty())
    {
      double emax = 0.0;
      for (const auto &p : calc_ref)
      {
        const double e = relative_error(p.first, p.second);
        if (std::isnan(e))
        {
          emax = e;
          break;
        }
        emax = std::max(emax, e);
      }
      fill_relative(ne, emax);
    }

    // Omission error of the undulation N, on the model sphere R (no upward attenuation).
    OmissionError &oe = rep.omission;
    oe.nmax = nmax;
    oe.n_tail = opt.n_tail;
    oe.kaula = omission_kaula_undulation(sigma_model.R, nmax);
    oe.tr = omission_tscherning_rapp_undulation(nmax, opt.n_tail, gamma, opt.tr);
    oe.tr_model = kNaN;
    if (opt.coefficients && opt.coefficients->nmax > nmax && opt.n_tail > opt.coefficients->nmax)
    {
      const Model &coef = *opt.coefficients;
      const std::vector<double> sig2 = degree_variances_with_tail(coef, opt.n_tail, opt.tr);
      oe.tr_model = omission_error(sig2, nmax, ErrorFunctional::ZetaEll, coef.GM, coef.R, coef.R, gamma);
      oe.has_tr_model = true;
      oe.coefficient_nmax = coef.nmax;
    }

    // Commission error: the sigma model cannot give degrees above its own maximum.
    const int nce = std::min(nmax, sigma_model.nmax);
    rep.commission.nmax = nce;
    rep.commission.global_N = commission_global(sigma_model, ErrorFunctional::ZetaEll, r, gamma, 2, nce);
    if (pts && pts->size() > 0)
    {
      ErrorOptions eo;
      eo.nmin = 2;
      eo.nmax = nce;
      const std::vector<double> ce = commission_error(sigma_model, *pts, ErrorFunctional::ZetaEll, eo);
      CommissionError &c = rep.commission;
      c.has_points = true;
      c.count = ce.size();
      c.min = ce.front();
      c.max = ce.front();
      double sum = 0.0;
      for (double x : ce)
      {
        c.min = std::min(c.min, x);
        c.max = std::max(c.max, x);
        sum += x;
      }
      c.mean = sum / static_cast<double>(ce.size());
    }
    return rep;
  }

  void set_relative_error(ErrorReport &report, double max_relative_error)
  {
    fill_relative(report.numerical, max_relative_error);
  }

  std::string format_error_report(const ErrorReport &rep)
  {
    std::ostringstream o;
    const NumericalError &ne = rep.numerical;
    const OmissionError &oe = rep.omission;
    const CommissionError &ce = rep.commission;

    o << "Error report: " << (rep.model_name.empty() ? std::string("(unnamed model)") : rep.model_name) << "\n";
    o << "  model degree " << rep.model_nmax << ", run nmax " << rep.nmax << ", R = " << fmt(rep.R, 1, false)
      << " m, GM = " << fmt(rep.GM, 6, true) << " m^3/s^2, gamma = " << fmt(rep.gamma, 4, false) << " m/s^2\n";

    o << "Numerical error (NE):\n";
    if (ne.available)
      o << "  E_rel = " << fmt(ne.max_relative_error, 3, true) << ", E_rel <= 1e-12: " << (ne.within_1e12 ? "yes" : "no")
        << ", digits retained = " << fmt(ne.digits, 2, false) << " (cap 16)\n";
    else
      o << "  E_rel: not computed (no calc/ref values given)\n";
    o << "  digit-loss bound (M = " << rep.nmax << ", colatitude " << fmt(rad2deg(ne.colatitude), 1, false)
      << " deg): " << fmt(ne.digit_loss_bound, 2, false) << " digits\n";

    o << "Omission error (undulation, degrees " << rep.nmax + 1 << ".." << oe.n_tail << ", m):\n";
    o << "  Kaula rule                 : " << fmt(oe.kaula, 4, false) << "\n";
    o << "  Tscherning-Rapp            : " << fmt(oe.tr, 4, false) << "\n";
    if (oe.has_tr_model)
      o << "  Tscherning-Rapp (model to " << oe.coefficient_nmax << ") : " << fmt(oe.tr_model, 4, false) << "\n";
    else
      o << "  Tscherning-Rapp (model)    : not computed (needs a coefficient model with degree > nmax)\n";

    o << "Commission error (undulation, degrees 2.." << ce.nmax << ", m):\n";
    o << "  CE(N) global               : " << fmt(ce.global_N, 4, false) << "\n";
    if (ce.has_points)
      o << "  points (" << ce.count << ")               : min " << fmt(ce.min, 4, false) << ", mean " << fmt(ce.mean, 4, false) << ", max "
        << fmt(ce.max, 4, false) << "\n";
    else
      o << "  points                     : not computed\n";
    return o.str();
  }

  std::string error_report_json(const ErrorReport &rep)
  {
    const NumericalError &ne = rep.numerical;
    const OmissionError &oe = rep.omission;
    const CommissionError &ce = rep.commission;
    std::ostringstream o;
    o << "{";
    o << "\"model\":{\"name\":\"" << json_escape(rep.model_name) << "\",\"GM_m3_s2\":" << json_number(rep.GM) << ",\"R_m\":" << json_number(rep.R)
      << ",\"model_nmax\":" << rep.model_nmax << "},";
    o << "\"run\":{\"nmax\":" << rep.nmax << ",\"r_m\":" << json_number(rep.r) << ",\"gamma_m_s2\":" << json_number(rep.gamma) << "},";
    o << "\"numerical\":{\"available\":" << (ne.available ? "true" : "false") << ",\"max_relative_error\":" << json_number(ne.max_relative_error)
      << ",\"within_1e-12\":" << json_bool(ne.available, ne.within_1e12) << ",\"digits\":" << json_number(ne.digits)
      << ",\"digit_loss_bound_digits\":" << json_number(ne.digit_loss_bound) << ",\"colatitude_rad\":" << json_number(ne.colatitude) << "},";
    o << "\"omission\":{\"nmax\":" << oe.nmax << ",\"n_tail\":" << oe.n_tail << ",\"kaula_m\":" << json_number(oe.kaula)
      << ",\"tr_m\":" << json_number(oe.tr) << ",\"tr_model_available\":" << (oe.has_tr_model ? "true" : "false")
      << ",\"tr_model_m\":" << json_number(oe.has_tr_model ? oe.tr_model : kNaN) << ",\"coefficient_nmax\":"
      << (oe.has_tr_model ? std::to_string(oe.coefficient_nmax) : std::string("null")) << "},";
    o << "\"commission\":{\"nmax\":" << ce.nmax << ",\"global_CE_N_m\":" << json_number(ce.global_N) << ",\"points\":";
    if (ce.has_points)
      o << "{\"count\":" << ce.count << ",\"min_m\":" << json_number(ce.min) << ",\"mean_m\":" << json_number(ce.mean)
        << ",\"max_m\":" << json_number(ce.max) << "}";
    else
      o << "null";
    o << "}}";
    return o.str();
  }
}
