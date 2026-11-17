// SCIP custom branching-rule demo (SCIP 10.1.0, C API used from C++).
//
// Implements a custom branching rule "owi" (Objective-Weighted Infeasibility):
//   among the LP-fractional branching candidates, select the variable that
//   maximises  (|obj| + 1) * min(frac, 1 - frac).
// This is a deliberately custom score, different from SCIP's built-in
// mostinf / fullstrong / relpscost rules.
//
// The instance is a small 3-constraint 0-1 ILP whose LP relaxation has
// several fractional variables and an integrality gap much larger than 1.
// That matters: SCIP detects an all-integral objective and rounds the dual
// bound down to an integer, so a single-knapsack gap below 1 is closed at
// the root and no branching ever happens. With this instance branching is
// unavoidable and the custom rule has real candidates to choose from.
//
// The demo solves the instance twice:
//   (A) with the custom rule at high priority,
//   (B) with the custom rule suppressed (SCIP defaults),
// prints each branching decision, verifies both solutions against an
// independent brute-force optimum, and reports the node counts so the
// effect of the custom rule is visible.
//
// Build (in CI only, via the GitHub Actions workflow):
//   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
//   ./build/scip_branching_demo

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>

namespace {

struct Item {
   double profit;
   double w[3];
};

// Three "cyclic" resource types. Adding the three capacity constraints gives
// 4 * sum(x) <= 18, i.e. sum(x) <= 4.5, while an integer solution can contain
// at most 4 items. The LP relaxation therefore stays fractional (three
// fractional basic variables) and its bound is far above the integer optimum.
const std::vector<Item> kItems = {
   {20, {2, 1, 1}}, {21, {2, 1, 1}}, {22, {2, 1, 1}}, {23, {2, 1, 1}},
   {20, {1, 2, 1}}, {21, {1, 2, 1}}, {22, {1, 2, 1}}, {23, {1, 2, 1}},
   {20, {1, 1, 2}}, {21, {1, 1, 2}}, {22, {1, 1, 2}}, {23, {1, 1, 2}},
};
const double kCapacity[3] = {6.0, 6.0, 6.0};

constexpr const char* kRuleName    = "owi";
constexpr int         kPriorityOn  = 1000000;
constexpr int         kPriorityOff = -1000000;
constexpr int         kTraceLimit  = 40;

// Statistics of the custom rule, collected during a solve.
struct RuleStats {
   long long calls    = 0;
   long long branched = 0;
};
RuleStats g_stats;
bool      g_verbose = true;

// ---------------------------------------------------------------------------
// Custom branching rule callbacks
// ---------------------------------------------------------------------------

static SCIP_DECL_BRANCHEXECLP(branchExecLpOwi)
{
   SCIP_VAR** cands     = nullptr;
   SCIP_Real* candsfrac = nullptr;
   int        ncands    = 0;

   SCIP_CALL( SCIPgetLPBranchCands(scip, &cands, nullptr, &candsfrac,
                                   nullptr, &ncands, nullptr) );

   ++g_stats.calls;

   double bestscore = -1.0;
   int    bestcand  = -1;
   double bestfrac  = 0.0;
   double bestobj   = 0.0;

   for( int i = 0; i < ncands; ++i )
   {
      const double frac   = candsfrac[i];                   // in [0,1)
      const double infeas = std::min(frac, 1.0 - frac);    // "most fractional" part
      const double obj    = std::fabs(SCIPvarGetObj(cands[i]));
      const double score  = (obj + 1.0) * infeas;          // <-- custom score

      if( score > bestscore )
      {
         bestscore = score;
         bestcand  = i;
         bestfrac  = frac;
         bestobj   = obj;
      }
   }

   if( bestcand < 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( g_verbose && g_stats.calls <= kTraceLimit )
   {
      std::cout << "  [owi] node=" << SCIPgetNNodes(scip)
                << " ncands=" << ncands
                << " -> branch on <" << SCIPvarGetName(cands[bestcand]) << ">"
                << " frac=" << bestfrac
                << " |obj|=" << bestobj
                << " score=" << bestscore << '\n';
   }

   SCIP_CALL( SCIPbranchVar(scip, cands[bestcand], nullptr, nullptr, nullptr) );
   ++g_stats.branched;
   *result = SCIP_BRANCHED;
   return SCIP_OKAY;
}

static SCIP_RETCODE includeCustomRule(SCIP* scip, int priority)
{
   SCIP_BRANCHRULE* rule = nullptr;
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &rule, kRuleName,
      "custom objective-weighted infeasibility branching rule",
      priority, -1, 1.0, nullptr) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, rule, branchExecLpOwi) );
   return SCIP_OKAY;
}

// ---------------------------------------------------------------------------
// Model building / solving
// ---------------------------------------------------------------------------

static void trySetInt(SCIP* scip, const char* name, int value)
{
   const SCIP_RETCODE rc = SCIPsetIntParam(scip, name, value);
   if( rc != SCIP_OKAY )
      std::cout << "  [warn] cannot set " << name << " (rc=" << rc << ")\n";
}

// Independent brute-force optimum (verification source of truth).
static double bruteForceOptimum(const std::vector<Item>& items,
                                const double cap[3])
{
   const std::size_t   n     = items.size();
   const std::uint64_t total = std::uint64_t(1) << n;
   double best = 0.0;
   for( std::uint64_t mask = 0; mask < total; ++mask )
   {
      double w[3] = {0.0, 0.0, 0.0};
      double p    = 0.0;
      for( std::size_t i = 0; i < n; ++i )
      {
         if( mask & (std::uint64_t(1) << i) )
         {
            w[0] += items[i].w[0];
            w[1] += items[i].w[1];
            w[2] += items[i].w[2];
            p    += items[i].profit;
         }
      }
      if( w[0] <= cap[0] && w[1] <= cap[1] && w[2] <= cap[2] && p > best )
         best = p;
   }
   return best;
}

struct SolveResult {
   SCIP_STATUS status       = SCIP_STATUS_UNKNOWN;
   double      objective    = 0.0;
   long long   nodes        = 0;
   long long   ruleCalls    = 0;
   long long   ruleBranched = 0;
   std::vector<double> solution;
};

static SCIP_RETCODE runSolve(int priority, bool verbose, SolveResult& out)
{
   SCIP* scip = nullptr;
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   // Keep the LP relaxation gap so the branching rule is actually exercised,
   // and keep SCIP quiet so our own trace is readable.
   trySetInt(scip, "display/verblevel", 0);
   trySetInt(scip, "presolving/maxrounds", 0);
   trySetInt(scip, "separating/maxroundsroot", 0);
   trySetInt(scip, "separating/maxrounds", 0);
   trySetInt(scip, "propagating/maxroundsroot", 0);
   trySetInt(scip, "propagating/maxrounds", 0);

   g_stats   = RuleStats{};
   g_verbose = verbose;
   SCIP_CALL( includeCustomRule(scip, priority) );

   SCIP_CALL( SCIPcreateProbBasic(scip, "custom_branching_mip") );
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   const int n = static_cast<int>(kItems.size());
   std::vector<SCIP_VAR*> vars(n, nullptr);

   for( int i = 0; i < n; ++i )
   {
      const std::string name = "x" + std::to_string(i + 1);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name.c_str(),
                                    0.0, 1.0, kItems[i].profit,
                                    SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   for( int k = 0; k < 3; ++k )
   {
      std::vector<SCIP_Real> coefs(n, 0.0);
      for( int i = 0; i < n; ++i )
         coefs[i] = kItems[i].w[k];

      const std::string cname = "cap" + std::to_string(k + 1);
      SCIP_CONS* cons = nullptr;
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, cname.c_str(), n,
                                           vars.data(), coefs.data(),
                                           -SCIPinfinity(scip), kCapacity[k]) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIP_CALL( SCIPsolve(scip) );

   out.status       = SCIPgetStatus(scip);
   out.nodes        = SCIPgetNNodes(scip);
   out.ruleCalls    = g_stats.calls;
   out.ruleBranched = g_stats.branched;

   SCIP_SOL* sol = SCIPgetBestSol(scip);
   if( sol != nullptr )
   {
      out.objective = SCIPgetSolOrigObj(scip, sol);
      out.solution.resize(n);
      for( int i = 0; i < n; ++i )
         out.solution[i] = SCIPgetSolVal(scip, sol, vars[i]);
   }

   for( int i = 0; i < n; ++i )
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   SCIP_CALL( SCIPfree(&scip) );
   return SCIP_OKAY;
}

} // namespace

int main()
{
   const double expected = bruteForceOptimum(kItems, kCapacity);

   std::cout << "SCIP custom branching-rule demo\n";
   std::cout << "  rule      : 'owi' = (|obj|+1) * min(frac, 1-frac)\n";
   std::cout << "  instance  : 3-constraint 0-1 ILP, n=" << kItems.size()
             << ", capacities=(" << kCapacity[0] << ", " << kCapacity[1]
             << ", " << kCapacity[2] << ")\n";
   std::cout << "  brute force optimum = " << expected << "\n\n";

   std::cout << "=== (A) custom rule enabled (priority " << kPriorityOn << ") ===\n";
   SolveResult a;
   if( runSolve(kPriorityOn, true, a) != SCIP_OKAY )
   {
      std::cerr << "SCIP error in solve A\n";
      return 1;
   }
   std::cout << "  ... (trace capped at " << kTraceLimit << " decisions)\n";
   std::cout << "  status=" << a.status << " objective=" << a.objective
             << " nodes=" << a.nodes
             << " rule_calls=" << a.ruleCalls
             << " rule_branched=" << a.ruleBranched << "\n\n";

   std::cout << "=== (B) custom rule suppressed (SCIP defaults) ===\n";
   SolveResult b;
   if( runSolve(kPriorityOff, false, b) != SCIP_OKAY )
   {
      std::cerr << "SCIP error in solve B\n";
      return 1;
   }
   std::cout << "  status=" << b.status << " objective=" << b.objective
             << " nodes=" << b.nodes << "\n\n";

   bool ok = true;
   if( a.status != SCIP_STATUS_OPTIMAL || b.status != SCIP_STATUS_OPTIMAL )
   {
      std::cerr << "ERROR: a solve did not reach optimality\n";
      ok = false;
   }
   if( std::fabs(a.objective - expected) > 1e-6 ||
       std::fabs(b.objective - expected) > 1e-6 )
   {
      std::cerr << "ERROR: objective does not match the brute-force optimum\n";
      ok = false;
   }

   double w[3] = {0.0, 0.0, 0.0};
   double p    = 0.0;
   for( std::size_t i = 0; i < a.solution.size(); ++i )
   {
      if( a.solution[i] > 0.5 )
      {
         w[0] += kItems[i].w[0];
         w[1] += kItems[i].w[1];
         w[2] += kItems[i].w[2];
         p    += kItems[i].profit;
      }
   }
   std::cout << "  solution resource usage = (" << w[0] << ", " << w[1]
             << ", " << w[2] << ") / (" << kCapacity[0] << ", "
             << kCapacity[1] << ", " << kCapacity[2] << ")  profit=" << p << '\n';
   if( w[0] > kCapacity[0] + 1e-6 ||
       w[1] > kCapacity[1] + 1e-6 ||
       w[2] > kCapacity[2] + 1e-6 ||
       std::fabs(p - expected) > 1e-6 )
   {
      std::cerr << "ERROR: reported solution is infeasible or suboptimal\n";
      ok = false;
   }

   std::cout << "\n  nodes  custom=" << a.nodes << "  default=" << b.nodes
             << "  (difference shows the custom rule's effect)\n";
   if( a.ruleBranched == 0 )
   {
      std::cerr << "ERROR: the custom branching rule was never executed\n";
      ok = false;
   }

   if( !ok )
   {
      std::cerr << "FAILED\n";
      return 1;
   }

   std::cout << "\nBUILD_AND_RUN_OK\n";
   return 0;
}
