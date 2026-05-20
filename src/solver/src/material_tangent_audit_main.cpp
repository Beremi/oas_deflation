#include "element.h"
#include "material_csl.h"
#include "material_ldpm.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std :: filesystem;

fs :: path GlobPaths :: BASEDIR;
Model *masterModel = nullptr;

namespace {

struct AuditRow {
    std :: string material;
    std :: string statusName;
    std :: string caseName;
    std :: string evaluation;
    std :: string branch;
    bool supported = true;
    double eps = std :: numeric_limits< double > :: quiet_NaN();
    double relerr = std :: numeric_limits< double > :: quiet_NaN();
    double cosine = std :: numeric_limits< double > :: quiet_NaN();
    double kpNorm = std :: numeric_limits< double > :: quiet_NaN();
    double fdNorm = std :: numeric_limits< double > :: quiet_NaN();
    double mismatch = std :: numeric_limits< double > :: quiet_NaN();
    std :: string note;
};

struct StepState {
    Vector strain;
    double volumetricStrain = 0.;
};

struct AuditCase {
    std :: string name;
    bool frozen = false;
    std :: vector< StepState > acceptedHistory;
    StepState trial;
    Vector direction;
};

struct MaterialSpec {
    std :: string name;
    std :: vector< std :: string > supportedBranches;
    std :: vector< std :: string > unsupportedBranches;
    std :: vector< AuditCase > cases;
    std :: unique_ptr< Material > ( *makeMaterial )();
};

Vector vec3(double a, double b, double c) {
    Vector v = Vector :: Zero(3);
    v [ 0 ] = a;
    v [ 1 ] = b;
    v [ 2 ] = c;
    return v;
}

Vector normalized(Vector v) {
    const double n = v.norm();
    if ( n > 0. ) {
        v /= n;
    }
    return v;
}

std :: unique_ptr< Material > makeCSLMaterial() {
    auto mat = std :: make_unique< CSLMaterial >(3);
    std :: istringstream iss("E0 3.000000e10 alpha 0.300000 density 2200.000000 ft 2.000000e6 Gt 500.000000");
    mat->readFromLine(iss);
    mat->setId(0);
    mat->init(nullptr);
    return mat;
}

std :: unique_ptr< Material > makeLDPMMaterial() {
    auto mat = std :: make_unique< LDPMMaterial >(3);
    std :: istringstream iss("E0 3.100000e10 alpha 0.300000 density 2200.000000 ft 2.000000e7 Gt 5000.000000 fc0 2.400000e9 fs0 3.000000e9 damage_residuum 0");
    mat->readFromLine(iss);
    mat->setId(0);
    mat->init(nullptr);
    return mat;
}

class AuditStatus {
public:
    explicit AuditStatus(std :: unique_ptr< Material > inputMaterial)
        : material(std :: move(inputMaterial)), element(3) {
        status.reset(material->giveNewMaterialStatus(&element, 0) );
        status->initializeStressAndStrainVector();
        status->init();
    }

    MaterialStatus &giveStatus() {
        return *status;
    }

private:
    std :: unique_ptr< Material > material;
    MaterialTestElement element;
    std :: unique_ptr< MaterialStatus > status;
};

void evaluate(MaterialStatus &status, const StepState &state, bool frozen, bool update) {
    status.setParameterValue("volumetric_strain", state.volumetricStrain);
    status.setTotalTempStrain(state.strain);
    if ( frozen ) {
        status.computeStressWithFrozenIntVars(1.);
    } else {
        status.computeStress(1.);
    }
    if ( update ) {
        status.update();
    }
}

void applyAcceptedHistory(MaterialStatus &status, const std :: vector< StepState > &history) {
    for ( const StepState &state : history ) {
        evaluate(status, state, false, true);
    }
    status.resetTemporaryVariables();
}

AuditRow runBranchAudit(const MaterialSpec &spec, const AuditCase &testCase, const std :: string &branch, double eps) {
    AuditStatus audit(spec.makeMaterial() );
    MaterialStatus &status = audit.giveStatus();
    applyAcceptedHistory(status, testCase.acceptedHistory);

    const auto acceptedSnapshot = status.cloneState();
    const std :: uint64_t acceptedHash = status.stateHash();

    evaluate(status, testCase.trial, testCase.frozen, false);
    const Vector baseStress = status.giveTempStress();
    const Matrix tangent = status.giveStiffnessTensor(branch);
    const Vector kp = tangent * testCase.direction;

    status.restoreStateFrom( *acceptedSnapshot );
    StepState perturbed = testCase.trial;
    perturbed.strain += eps * testCase.direction;
    evaluate(status, perturbed, testCase.frozen, false);
    const Vector fd = ( status.giveTempStress() - baseStress ) / eps;

    status.restoreStateFrom( *acceptedSnapshot );
    if ( status.stateHash() != acceptedHash ) {
        std :: cerr << "material snapshot hash mismatch after restore in " << spec.name << " / " << testCase.name << '\n';
        exit(EXIT_FAILURE);
    }

    AuditRow row;
    row.material = spec.name;
    row.statusName = status.giveName();
    row.caseName = testCase.name;
    row.evaluation = testCase.frozen ? "frozen" : "actual";
    row.branch = branch;
    row.eps = eps;
    row.kpNorm = kp.norm();
    row.fdNorm = fd.norm();
    row.mismatch = ( kp - fd ).norm();
    row.relerr = row.mismatch / std :: max(row.fdNorm, 1e-30);
    const double cosineDenom = row.kpNorm * row.fdNorm;
    row.cosine = cosineDenom > 0. ? kp.dot(fd) / cosineDenom : 0.;
    return row;
}

AuditRow unsupportedRow(const MaterialSpec &spec, const AuditCase &testCase, const std :: string &branch) {
    AuditStatus audit(spec.makeMaterial() );
    AuditRow row;
    row.material = spec.name;
    row.statusName = audit.giveStatus().giveName();
    row.caseName = testCase.name;
    row.evaluation = testCase.frozen ? "frozen" : "actual";
    row.branch = branch;
    row.supported = false;
    row.note = "unsupported";
    return row;
}

std :: vector< AuditCase > cslCases() {
    const Vector d = normalized(vec3(1., 0.25, -0.10) );
    return {
        { "elastic_loading", false, {}, { vec3(2.0e-5, 3.0e-6, -2.0e-6), 2.0e-5 }, d },
        { "damage_growth", false, {}, { vec3(9.0e-5, 1.5e-5, -8.0e-6), 9.0e-5 }, d },
        { "damage_growth_frozen", true, {}, { vec3(9.0e-5, 1.5e-5, -8.0e-6), 9.0e-5 }, d },
        { "damaged_unloading", false, { { vec3(1.2e-4, 2.0e-5, -1.0e-5), 1.2e-4 } }, { vec3(5.0e-5, 7.0e-6, -3.0e-6), 5.0e-5 }, d },
        { "damaged_unloading_frozen", true, { { vec3(1.2e-4, 2.0e-5, -1.0e-5), 1.2e-4 } }, { vec3(5.0e-5, 7.0e-6, -3.0e-6), 5.0e-5 }, d },
    };
}

std :: vector< AuditCase > ldpmCases() {
    const Vector d = normalized(vec3(1., 0.20, -0.08) );
    return {
        { "elastic_loading", false, {}, { vec3(1.0e-4, 1.0e-5, -5.0e-6), 1.0e-4 }, d },
        { "damage_growth", false, {}, { vec3(9.0e-4, 1.2e-4, -5.0e-5), 9.0e-4 }, d },
        { "damage_growth_frozen", true, {}, { vec3(9.0e-4, 1.2e-4, -5.0e-5), 9.0e-4 }, d },
        { "damaged_unloading", false, { { vec3(1.1e-3, 1.4e-4, -6.0e-5), 1.1e-3 } }, { vec3(4.0e-4, 4.0e-5, -2.0e-5), 4.0e-4 }, d },
        { "damaged_unloading_frozen", true, { { vec3(1.1e-3, 1.4e-4, -6.0e-5), 1.1e-3 } }, { vec3(4.0e-4, 4.0e-5, -2.0e-5), 4.0e-4 }, d },
    };
}

void writeRows(std :: ostream &out, const std :: vector< AuditRow > &rows) {
    out << "material\tstatus\tcase\tevaluation\tbranch\tsupported\teps\trelerr\tcosine\tkp_norm\tfd_norm\tmismatch\tnote\n";
    out << std :: scientific << std :: setprecision(10);
    for ( const AuditRow &row : rows ) {
        out << row.material << '\t'
            << row.statusName << '\t'
            << row.caseName << '\t'
            << row.evaluation << '\t'
            << row.branch << '\t'
            << ( row.supported ? 1 : 0 ) << '\t';
        if ( row.supported ) {
            out << row.eps << '\t'
                << row.relerr << '\t'
                << row.cosine << '\t'
                << row.kpNorm << '\t'
                << row.fdNorm << '\t'
                << row.mismatch << '\t';
        } else {
            out << "nan\tnan\tnan\tnan\tnan\tnan\t";
        }
        out << ( row.note.empty() ? "ok" : row.note ) << '\n';
    }
}

bool rowMatches(const AuditRow &row, const std :: string &material, const std :: string &caseName, const std :: string &branch) {
    return row.supported && row.material == material && row.caseName == caseName && row.branch == branch;
}

bool expectedLinearPathFailed(const AuditRow &row) {
    constexpr double tolerance = 5e-4;
    if ( rowMatches(row, "CSLMaterial", "elastic_loading", "secant") ||
         rowMatches(row, "CSLMaterial", "elastic_loading", "tangent") ||
         rowMatches(row, "CSLMaterial", "elastic_loading", "archived_csl_damage_tangent") ||
         rowMatches(row, "CSLMaterial", "elastic_loading", "consistent") ||
         rowMatches(row, "CSLMaterial", "damage_growth", "archived_csl_damage_tangent") ||
         rowMatches(row, "CSLMaterial", "damage_growth", "consistent") ||
         rowMatches(row, "CSLMaterial", "damage_growth_frozen", "secant") ||
         rowMatches(row, "CSLMaterial", "damage_growth_frozen", "tangent") ||
         rowMatches(row, "CSLMaterial", "damage_growth_frozen", "archived_csl_damage_tangent") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading", "secant") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading", "tangent") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading", "archived_csl_damage_tangent") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading", "consistent") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading_frozen", "secant") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading_frozen", "tangent") ||
         rowMatches(row, "CSLMaterial", "damaged_unloading_frozen", "archived_csl_damage_tangent") ||
         rowMatches(row, "LDPMMaterial", "elastic_loading", "consistent") ||
         rowMatches(row, "LDPMMaterial", "damage_growth", "consistent") ||
         rowMatches(row, "LDPMMaterial", "damaged_unloading", "consistent") ||
         rowMatches(row, "LDPMMaterial", "damage_growth_frozen", "secant") ||
         rowMatches(row, "LDPMMaterial", "damaged_unloading_frozen", "secant") ) {
        return row.relerr > tolerance || !std :: isfinite(row.relerr);
    }
    return false;
}

} // namespace

int main(int argc, char **argv) {
    fs :: path outputPath;
    for ( int i = 1; i < argc; i++ ) {
        const std :: string arg(argv [ i ]);
        if ( arg == "--output" && i + 1 < argc ) {
            outputPath = argv [ ++i ];
        } else {
            std :: cerr << "Usage: " << argv [ 0 ] << " [--output audit.tsv]\n";
            return EXIT_FAILURE;
        }
    }

    const double eps = 1e-8;
    const std :: vector< MaterialSpec > specs = {
        { "CSLMaterial", { "elastic", "secant", "tangent", "archived_csl_damage_tangent", "consistent" }, {}, cslCases(), makeCSLMaterial },
        { "LDPMMaterial", { "elastic", "secant", "consistent" }, { "tangent" }, ldpmCases(), makeLDPMMaterial },
    };

    std :: vector< AuditRow > rows;
    for ( const MaterialSpec &spec : specs ) {
        for ( const AuditCase &testCase : spec.cases ) {
            for ( const std :: string &branch : spec.supportedBranches ) {
                if ( branch == "consistent" && testCase.frozen ) {
                    AuditRow row = unsupportedRow(spec, testCase, branch);
                    row.note = "not_applicable_to_frozen_evaluation";
                    rows.push_back(row);
                } else {
                    rows.push_back(runBranchAudit(spec, testCase, branch, eps) );
                }
            }
            for ( const std :: string &branch : spec.unsupportedBranches ) {
                rows.push_back(unsupportedRow(spec, testCase, branch) );
            }
        }
    }

    writeRows(std :: cout, rows);
    if ( !outputPath.empty() ) {
        if ( outputPath.has_parent_path() ) {
            fs :: create_directories(outputPath.parent_path() );
        }
        std :: ofstream output(outputPath);
        if ( !output.is_open() ) {
            std :: cerr << "Cannot open audit output " << outputPath << '\n';
            return EXIT_FAILURE;
        }
        writeRows(output, rows);
    }

    bool failed = false;
    for ( const AuditRow &row : rows ) {
        if ( expectedLinearPathFailed(row) ) {
            std :: cerr << "Unexpected tangent FD error in linear/frozen path: "
                       << row.material << ' ' << row.caseName << ' ' << row.branch
                       << " relerr=" << row.relerr << '\n';
            failed = true;
        }
    }
    if ( failed ) {
        return EXIT_FAILURE;
    }

    std :: cerr << "material tangent audit passed expected linear/frozen checks\n";
    return EXIT_SUCCESS;
}
