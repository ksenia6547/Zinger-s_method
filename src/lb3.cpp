#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>
#include <stdexcept>
#include <vector>

static const double PI = std::acos(-1.0);
static const double DD2R = PI / 180.0;
static const double DR2D = 180.0 / PI;
static const double DH2R = PI / 12.0;
static const double DR2H = 12.0 / PI;

struct Params {
    double phi;
    double ae;
    double de;
    double aw;
    double dw;
};

struct Bracket {
    double a;
    double b;
    bool found;
};

double F(double s, const Params& p)
{
    return std::sin(p.phi) * (std::sin(p.dw) - std::sin(p.de))
         + std::cos(p.phi) * (std::cos(p.dw) * std::cos(s - p.aw)
                               - std::cos(p.de) * std::cos(s - p.ae));
}

double dF(double s, const Params& p)
{
    return std::cos(p.phi) * (- std::cos(p.dw) * std::sin(s - p.aw)
                               + std::cos(p.de) * std::sin(s - p.ae));
}

double zenithDist(double s0, double alpha, double delta, double phi)
{
    double cz = std::sin(phi) * std::sin(delta)
              + std::cos(phi) * std::cos(delta) * std::cos(s0 - alpha);
    cz = std::max(-1.0, std::min(1.0, cz));
    return std::acos(cz) * DR2D;
}

std::vector<Bracket> findAllBrackets(const Params& p, int steps = 5000)
{
    std::vector<Bracket> result;
    double h = 2.0 * PI / steps;
    double x0 = 0.0;
    double f0 = F(x0, p);

    for (int i = 1; i <= steps; ++i) {
        double x1 = i * h;
        double f1 = F(x1, p);
        if (std::signbit(f0) != std::signbit(f1))
            result.push_back({x0, x1, true});
        x0 = x1;
        f0 = f1;
    }
    return result;
}

double bisection(double a, double b, const Params& p)
{
    const double EPS = std::numeric_limits<double>::epsilon();
    const int MAXIT = 300;
    double fa = F(a, p);

    std::cout << std::left
              << std::setw(5) << "Iter"
              << std::setw(26) << "  a"
              << std::setw(26) << "  b"
              << std::setw(26) << "  mid"
              << std::setw(20) << "  f(mid)"
              << "  длина отрезка\n";

    double mid = a;
    for (int i = 1; i <= MAXIT; ++i) {
        mid = 0.5 * (a + b);
        double fm = F(mid, p);
        double len = b - a;

        std::cout << std::left << std::fixed
                  << std::setw(5) << i
                  << "  " << std::setw(24) << std::setprecision(17) << a
                  << "  " << std::setw(24) << b
                  << "  " << std::setw(24) << mid
                  << "  " << std::setw(18)
                  << std::scientific << std::setprecision(6) << fm
                  << "  " << len
                  << "\n";

        if (std::fabs(fm) == 0.0 || len < 2.0 * EPS * std::fabs(mid)) {
            std::cout << "\n Бисекция сошлась за " << i << " итераций.\n";
            return mid;
        }

        if (std::signbit(fa) != std::signbit(fm)) { b = mid; }
        else { a = mid; fa = fm; }
    }
    std::cout << "\n Достигнут лимит итераций.\n";
    return mid;
}

double newton(double x0, const Params& p)
{
    const double EPS = std::numeric_limits<double>::epsilon();
    const int MAXIT = 50;

    std::cout << std::left
              << std::setw(5) << "Iter"
              << std::setw(26) << "  x"
              << std::setw(26) << "  f(x)"
              << std::setw(26) << "  f'(x)"
              << "  шаг Δx = f/f'\n";

    double x = x0;
    for (int i = 1; i <= MAXIT; ++i) {
        double fx = F(x, p);
        double dfx = dF(x, p);

        if (std::fabs(dfx) < 1e-30)
            throw std::runtime_error(
                "Метод Ньютона расходится: производная близка к нулю.");

        double step = fx / dfx;
        double xnew = x - step;

        std::cout << std::left
                  << std::setw(5) << i
                  << "  " << std::setw(24)
                  << std::setprecision(17) << std::fixed    << x
                  << "  " << std::setw(24)
                  << std::scientific << std::setprecision(6) << fx
                  << "  " << std::setw(24) << dfx
                  << "  " << step
                  << "\n";

        if (std::fabs(step) <= EPS * std::fabs(xnew) + EPS) {
            std::cout << "\n Метод Ньютона сошёлся за " << i << " итераций.\n";
            return xnew;
        }
        x = xnew;
    }
    std::cout << "\n Достигнут лимит итераций.\n";
    return x;
}

void printResult(const std::string& method, double s0, const Params& p)
{
    s0 = std::fmod(s0, 2.0 * PI);
    if (s0 < 0.0) s0 += 2.0 * PI;

    double s0_hours = s0 * DR2H;
    int hh = static_cast<int>(s0_hours);
    int mm = static_cast<int>((s0_hours - hh) * 60.0);
    double ss = ((s0_hours - hh) * 60.0 - mm) * 60.0;

    double ze = zenithDist(s0, p.ae, p.de, p.phi);
    double zw = zenithDist(s0, p.aw, p.dw, p.phi);

    std::cout << "\n" << method << "\n";
    std::cout << "    s0      = "
              << std::fixed << std::setprecision(17) << s0 << " рад\n";
    std::cout << "    s0      = "
              << std::setprecision(10) << s0 * DR2D << " °\n";
    std::cout << "    s0      = "
              << hh << "ч " << mm << "м "
              << std::setprecision(6) << ss << "с\n";
    std::cout << "    ze      = "
              << std::setprecision(10) << ze << " ° (восточная)\n";
    std::cout << "    zw      = "
              << zw << " ° (западная)\n";
    std::cout << "    |ze−zw| = "
              << std::scientific << std::setprecision(3)
              << std::fabs(ze - zw) << " °\n";
    std::cout << "    f(s0)   = " << F(s0, p) << "\n";
}

void hmsToRad(double h, double m, double s, double& out)
{
    out = (h + m / 60.0 + s / 3600.0) * DH2R;
}

void dmsToRad(double d, double m, double s, double& out)
{
    double sign = (d < 0.0) ? -1.0 : 1.0;
    out = sign * (std::fabs(d) + m / 60.0 + s / 3600.0) * DD2R;
}

Params readParams()
{
    std::cout << "  0 — Альтаир, Вега (φ=20°)\n"
              << "  1 — W392, E713 (φ=60°31'57\")\n"
              << "  2 — W452, E655 (φ=60°31'57\")\n"
              << "  3 — W2547, E3488 (φ=60°31'57\")\n"
              << "  4 — ввод вручную\n"
              << "выбор: ";

    int choice;
    std::cin >> choice;

    Params p;
    double phi60 = (60.0 + 31.0/60.0 + 57.0/3600.0) * DD2R;

    if (choice == 0) {
        p.phi = 20.0 * DD2R;
        hmsToRad(19, 50, 47.0, p.ae);
        dmsToRad(8, 52, 0.0, p.de);
        hmsToRad(18, 36, 56.0, p.aw);
        dmsToRad(38, 47, 0.0, p.dw);
        std::cout << "\nВосточная: Альтаир α=19h50m47s δ=+8°52'\n"
                  << "Западная: Вега α=18h36m56s δ=+38°47'\n"
                  << "φ = 20°\n\n";
    } else if (choice == 1) {
        p.phi = phi60;
        hmsToRad(19, 45, 10.660, p.ae);
        dmsToRad(37, 24, 33.812, p.de);
        hmsToRad(10, 40, 9.742, p.aw);
        dmsToRad(31, 50, 41.176, p.dw);
        std::cout << "\nВосточная: E713 α=19h45m10.660s δ=+37°24'33.812\"\n"
                  << "Западная:  W392 α=10h40m09.742s δ=+31°50'41.176\"\n"
                  << "φ = 60°31'57\"\n\n";
    } else if (choice == 2) {
        p.phi = phi60;
        hmsToRad(17, 53, 58.541, p.ae);
        dmsToRad(56, 51, 40.293, p.de);
        hmsToRad(12, 16, 43.174, p.aw);
        dmsToRad(56, 53, 29.747, p.dw);
        std::cout << "\nВосточная: E655   α=17h53m58.541s  δ=+56°51'40.293\"\n"
                  << "Западная:  W452   α=12h16m43.174s  δ=+56°53'29.747\"\n"
                  << "φ = 60°31'57\"\n\n";
    } else if (choice == 3) {
        p.phi = phi60;
        hmsToRad(19, 25, 6.972, p.ae);
        dmsToRad(29, 39, 54.648, p.de);
        hmsToRad(12, 27, 41.393, p.aw);
        dmsToRad(27, 7, 35.226, p.dw);
        std::cout << "\nВосточная: E3488  α=19h25m06.972s  δ=+29°39'54.648\"\n"
                  << "Западная:  W2547  α=12h27m41.393s  δ=+27°07'35.226\"\n"
                  << "φ = 60°31'57\"\n\n";
    } else {
        double d, m, s;
        std::cout << "\nШирота φ: ";
        std::cin >> d >> m >> s;
        dmsToRad(d, m, s, p.phi);

        std::cout << "Восточная звезда:\n";
        std::cout << "  α (ч м с): ";
        std::cin >> d >> m >> s;
        hmsToRad(d, m, s, p.ae);
        std::cout << "  δ (° ' \"): ";
        std::cin >> d >> m >> s;
        dmsToRad(d, m, s, p.de);

        std::cout << "Западная звезда:\n";
        std::cout << "  α (ч м с): ";
        std::cin >> d >> m >> s;
        hmsToRad(d, m, s, p.aw);
        std::cout << "  δ (° ' \"): ";
        std::cin >> d >> m >> s;
        dmsToRad(d, m, s, p.dw);
    }

    return p;
}

int main()
{
    try {
        Params p = readParams();

        std::vector<Bracket> brackets = findAllBrackets(p);

        if (brackets.empty()) {
            std::cerr << "\nОшибка: нет момента, когда зенитные расстояния совпадают.\n";
            return 1;
        }

        Bracket chosen = brackets[0];
        double best_ze = 1e9;

        for (auto& br : brackets) {
            double mid = 0.5 * (br.a + br.b);
            double ze = zenithDist(mid, p.ae, p.de, p.phi);
            if (ze < best_ze) {
                best_ze = ze;
                chosen = br;
            }
        }

        if (best_ze >= 90.0) {
            std::cerr << "\nВсе найденные решения — звёзды под горизонтом\n";
        }

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Выбран отрезок: ["
                  << chosen.a << ", " << chosen.b << "] рад"
                  << "  (ze ~ " << std::setprecision(2) << best_ze << "°)\n";
        std::cout << "f(a) = " << std::scientific << F(chosen.a, p)
                  << ",  f(b) = " << F(chosen.b, p) << "\n";

        double x0 = 0.5 * (chosen.a + chosen.b);
        std::cout << "Начальное приближение x0 = "
                  << std::fixed << std::setprecision(8) << x0 << " рад\n";

        std::cout << "\n Метод бисекции:\n";
        double s0_bis = bisection(chosen.a, chosen.b, p);

        std::cout << "\n Метод Ньютона:\n";
        double s0_newt = newton(x0, p);

        std::cout << "\n Результат\n";
        printResult("Бисекция", s0_bis, p);
        printResult("Метод Ньютона", s0_newt, p);

        double diff = std::fabs(s0_bis - s0_newt);
        std::cout << "\nРасхождение между методами: "
                  << std::scientific << std::setprecision(3) << diff
                  << " рад (" << diff * DR2D << " °)\n";

        std::cout << "\nЗенитные расстояния двух звёзд равны.\n";

    } catch (const std::exception& e) {
        std::cerr << "\nИСКЛЮЧЕНИЕ: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
