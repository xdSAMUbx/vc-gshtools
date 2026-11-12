// errors.cpp: commission and omission errors (errors.hpp) and the consolidated error report (error_report.hpp).
//
// Units: the undulation and height anomaly values are in metres; gravity functionals in m/s^2; gamma in m/s^2;
// radii in metres; degree variances in (dimensionless coefficient)^2 units; Tscherning-Rapp amplitude A in mGal^2.
#include <cstdint>
#include <utility>
#include <vector>
#include "common.hpp"
#include "error_report.hpp"
#include "errors.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace gshpy
{
  void bind_errors(py::module_ &m)
  {
    using namespace gsh;

    // ---- enums and options -------------------------------------------------------------------------------------
    py::enum_<ErrorFunctional>(m, "ErrorFunctional", "Functional whose error is evaluated.")
        .value("T", ErrorFunctional::T)
        .value("ZetaEll", ErrorFunctional::ZetaEll)
        .value("GravityAnomaly", ErrorFunctional::GravityAnomaly)
        .value("GravityDisturbance", ErrorFunctional::GravityDisturbance);
    m.def("error_functional_name", [](ErrorFunctional f) { return std::string(error_functional_name(f)); }, py::arg("functional"));

    py::class_<ErrorOptions>(m, "ErrorOptions", "Degree range, ellipsoid and threads for commission_error.")
        .def(py::init<>())
        .def_readwrite("nmin", &ErrorOptions::nmin)
        .def_readwrite("nmax", &ErrorOptions::nmax, "< 0 = degree of the sigma model")
        .def_readwrite("ellipsoid", &ErrorOptions::ellipsoid)
        .def_readwrite("threads", &ErrorOptions::threads);

    py::class_<TscherningRapp>(m, "TscherningRapp", "Tscherning-Rapp degree-variance model parameters.")
        .def(py::init<>())
        .def_readwrite("A", &TscherningRapp::A, "amplitude [mGal^2]")
        .def_readwrite("B", &TscherningRapp::B)
        .def_readwrite("R_B", &TscherningRapp::R_B, "Bjerhammar sphere radius [m]")
        .def_readwrite("R_E", &TscherningRapp::R_E, "Earth radius [m]")
        .def("s", &TscherningRapp::s, "s = R_B^2 / R_E^2");

    // ---- commission and omission -------------------------------------------------------------------------------
    m.def("commission_error", [](const Model &sigma, const Points &pts, ErrorFunctional f, const ErrorOptions &opts) {
           std::vector<double> ce;
           {
             py::gil_scoped_release g;
             ce = commission_error(sigma, pts, f, opts);
           }
           return to_numpy(std::move(ce));
         },
          py::arg("sigma"), py::arg("pts"), py::arg("functional"), py::arg("opts") = ErrorOptions(),
          "One-sigma commission error per point, in the unit of the functional (metres for ZetaEll, m^2/s^2 for T).");

    m.def("commission_per_degree", [](const Model &sigma, ErrorFunctional f, double r, double gamma, int nmin, int nmax) {
           return to_numpy(commission_per_degree(sigma, f, r, gamma, nmin, nmax));
         },
          py::arg("sigma"), py::arg("functional"), py::arg("r"), py::arg("gamma"), py::arg("nmin"), py::arg("nmax"),
          "Area-mean commission contribution per degree (nmax + 1 entries).");

    m.def("commission_global", &commission_global, py::arg("sigma"), py::arg("functional"), py::arg("r"), py::arg("gamma"),
          py::arg("nmin"), py::arg("nmax"), "Global one-sigma commission error (same unit as the functional).");

    m.def("degree_variances", [](const Model &model) { return to_numpy(degree_variances(model)); }, py::arg("model"),
          "Degree variances sigma_n^2 = sum_m (C^2 + S^2), n = 0..model.nmax.");

    m.def("degree_variances_with_tail", [](const Model &model, int n_end, const TscherningRapp &p) {
           return to_numpy(degree_variances_with_tail(model, n_end, p));
         },
          py::arg("model"), py::arg("n_end"), py::arg("tr") = TscherningRapp(),
          "Degree variances of model, extended with the Tscherning-Rapp model up to n_end.");

    m.def("tscherning_rapp_variances", [](double GM, double R, int n_first, int n_last, const TscherningRapp &p) {
           return to_numpy(tscherning_rapp_variances(GM, R, n_first, n_last, p));
         },
          py::arg("GM"), py::arg("R"), py::arg("n_first"), py::arg("n_last"), py::arg("tr") = TscherningRapp(),
          "Degree variances implied by the Tscherning-Rapp model on the sphere of radius R.");

    m.def("omission_kaula_undulation", &omission_kaula_undulation, py::arg("R"), py::arg("nmax"),
          "Kaula-rule omission error of N [m]: R 1e-5 / (nmax + 1).");

    m.def("omission_tscherning_rapp_undulation", [](int nmax, int n_tail, double gamma, const TscherningRapp &p) {
           return omission_tscherning_rapp_undulation(nmax, n_tail, gamma, p);
         },
          py::arg("nmax"), py::arg("n_tail"), py::arg("gamma"), py::arg("tr") = TscherningRapp(),
          "Tscherning-Rapp omission error of N [m] for degrees nmax+1 .. n_tail.");

    m.def("omission_error", [](const DArr &sigma2, int nmax, ErrorFunctional f, double GM, double R, double r, double gamma) {
           return omission_error(copy_in<std::vector<double>>(sigma2, "sigma2"), nmax, f, GM, R, r, gamma);
         },
          py::arg("sigma2"), py::arg("nmax"), py::arg("functional"), py::arg("GM"), py::arg("R"), py::arg("r"),
          py::arg("gamma"), "Omission error sqrt(sum_{n>nmax} w_n^2 sigma_n^2) from degree variances sigma2.");

    m.def("relative_error", &relative_error, py::arg("calc"), py::arg("ref"), "|calc - ref| / |ref|");
    m.def("digit_loss_bound", &digit_loss_bound, py::arg("M"), py::arg("theta"),
          "Decimal digits lost by the Legendre recurrence at degree M and colatitude theta [rad].");

    // ---- error report ------------------------------------------------------------------------------------------
    py::class_<NumericalError>(m, "NumericalError", "Numerical error part of the report.")
        .def_readonly("available", &NumericalError::available)
        .def_readonly("max_relative_error", &NumericalError::max_relative_error)
        .def_readonly("within_1e12", &NumericalError::within_1e12)
        .def_readonly("digits", &NumericalError::digits)
        .def_readonly("colatitude", &NumericalError::colatitude)
        .def_readonly("digit_loss_bound", &NumericalError::digit_loss_bound);

    py::class_<OmissionError>(m, "OmissionError", "Omission error part of the report (undulation N, metres).")
        .def_readonly("nmax", &OmissionError::nmax)
        .def_readonly("n_tail", &OmissionError::n_tail)
        .def_readonly("kaula", &OmissionError::kaula)
        .def_readonly("tr", &OmissionError::tr)
        .def_readonly("has_tr_model", &OmissionError::has_tr_model)
        .def_readonly("tr_model", &OmissionError::tr_model)
        .def_readonly("coefficient_nmax", &OmissionError::coefficient_nmax);

    py::class_<CommissionError>(m, "CommissionError", "Commission error part of the report.")
        .def_readonly("nmax", &CommissionError::nmax)
        .def_readonly("global_N", &CommissionError::global_N)
        .def_readonly("has_points", &CommissionError::has_points)
        .def_readonly("count", &CommissionError::count)
        .def_readonly("min", &CommissionError::min)
        .def_readonly("mean", &CommissionError::mean)
        .def_readonly("max", &CommissionError::max);

    py::class_<ErrorReport>(m, "ErrorReport", "Consolidated error report of a run.")
        .def_readonly("model_name", &ErrorReport::model_name)
        .def_readonly("GM", &ErrorReport::GM)
        .def_readonly("R", &ErrorReport::R)
        .def_readonly("model_nmax", &ErrorReport::model_nmax)
        .def_readonly("nmax", &ErrorReport::nmax)
        .def_readonly("r", &ErrorReport::r)
        .def_readonly("gamma", &ErrorReport::gamma)
        .def_readonly("numerical", &ErrorReport::numerical)
        .def_readonly("omission", &ErrorReport::omission)
        .def_readonly("commission", &ErrorReport::commission);

    py::class_<ReportOptions>(m, "ReportOptions", "Optional settings of make_error_report (coefficients: see make_error_report).")
        .def(py::init<>())
        .def_readwrite("colatitude", &ReportOptions::colatitude, "colatitude [rad] for the digit-loss bound (default pi/2)")
        .def_readwrite("n_tail", &ReportOptions::n_tail, "upper degree of the Tscherning-Rapp sums (default 10800)")
        .def_readwrite("tr", &ReportOptions::tr);

    m.def("make_error_report", [](const Model &sigma_model, int nmax, double r, double gamma,
                                  const std::vector<std::pair<double, double>> &calc_ref, py::object pts_obj,
                                  ReportOptions opt, py::object coefficients_obj) {
           const Points *pts = nullptr;
           if (!pts_obj.is_none())
             pts = pts_obj.cast<const Points *>();
           const Model *coef = nullptr;
           if (!coefficients_obj.is_none())
             coef = coefficients_obj.cast<const Model *>();
           opt.coefficients = coef;
           ErrorReport rep;
           {
             py::gil_scoped_release g;
             rep = make_error_report(sigma_model, nmax, r, gamma, calc_ref, pts, opt);
           }
           return rep;
         },
          py::arg("sigma_model"), py::arg("nmax"), py::arg("r"), py::arg("gamma") = 9.81,
          py::arg("calc_ref") = std::vector<std::pair<double, double>>(), py::arg("pts") = py::none(),
          py::arg("opt") = ReportOptions(), py::arg("coefficients") = py::none(),
          "Builds the error report. sigma_model: load_gfc_sigma. r [m]: radius of the global commission value.\n"
          "calc_ref: optional (calc, ref) pairs for the numerical error. pts: optional points for point statistics.\n"
          "coefficients: optional load_gfc model (not sigma) for the Tscherning-Rapp model tail.");

    m.def("set_relative_error", &set_relative_error, py::arg("report"), py::arg("max_relative_error"),
          "Sets the numerical part from a precomputed maximum relative error.");
    m.def("format_error_report", &format_error_report, py::arg("report"), "Human-readable text (several lines).");
    m.def("error_report_json", &error_report_json, py::arg("report"), "Single-line JSON object.");
  }
}
