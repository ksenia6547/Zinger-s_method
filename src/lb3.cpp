#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>
#include <stdexcept>

static const double PI = std::acos(-1.0);
static const double DEG = PI/180.0;
static const double RAD = 180.0/PI;
static const double HOUR = PI/12.0;

struct Params {
    double phi;
    double ae;
    double de;
    double aw;
    double dw;
};

double F(double s, const Params& p)
{
    return std::sin(p.phi) * (std::sin(p.dw) - std::sin(p.de)) + std::cos(p.phi) * (  std::cos(p.dw) * std::cos(s - p.aw) - std::cos(p.de) * std::cos(s - p.ae));
}

double dF(double s, const Params& p)
{
    return std::cos(p.phi) * (- std::cos(p.dw) * std::sin(s - p.aw) + std::cos(p.de) * std::sin(s - p.ae));
}

struct Bracket { double a, b; bool found; };

Bracket findBracket(const Params& p, int steps = 5000)
{
    double h = 2.0 * PI / steps;
    double x0 = 0.0;
    double f0 = F(x0, p);

    for (int i = 1; i <= steps; ++i) {
        double x1 = i * h;
        double f1 = F(x1, p);
        if (f0 * f1 < 0.0)
            return {x0, x1, true};
        x0 = x1;
        f0 = f1;
    }
    return {0, 0, false};
}

double bisection(double a, double b, const Params& p)
{
    const double EPS = std::numeric_limits<double>::epsilon();
    const int MAXIT = 300;
    double fa = F(a, p);

    std::cout << std::left
              << std::setw(5)  << "Iter"
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

        std::cout << std::left  << std::fixed
                  << std::setw(5)  << i
                  << "  " << std::setw(24) << std::setprecision(17) << a
                  << "  " << std::setw(24) << b
                  << "  " << std::setw(24) << mid
                  << "  " << std::setw(18) << std::scientific << std::setprecision(6) << fm
                  << "  " << len
                  << "\n";

        if (std::abs(fm) == 0.0 || len < 2.0 * EPS * std::abs(mid)) {
            std::cout << "\n Бисекция сошлась за " << i << " итераций.\n";
            return mid;
        }

        if (fa * fm < 0.0) { b = mid;}
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
              << "  шаг  Δx = f/f'\n";

    double x = x0;
    for (int i = 1; i <= MAXIT; ++i) {
        double fx = F(x, p);
        double dfx = dF(x, p);

        if (std::abs(dfx) < 1e-30)
            throw std::runtime_error("Метод Ньютона расходится, т.к. производная близка к нулю.");

        double step = fx / dfx;
        double xnew = x - step;

        std::cout << std::left
                  << std::setw(5)  << i
                  << "  " << std::setw(24) << std::setprecision(17) << std::fixed    << x
                  << "  " << std::setw(24) << std::scientific << std::setprecision(6) << fx
                  << "  " << std::setw(24) << dfx
                  << "  " << step
                  << "\n";

        if (std::abs(step) <= EPS * std::abs(xnew) + EPS) {
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
    if (s0 < 0) s0 += 2.0 * PI;

    double s0_hours = s0 / HOUR;
    int hh = (int)s0_hours;
    int mm = (int)((s0_hours - hh) * 60.0);
    double ss = ((s0_hours - hh) * 60.0 - mm) * 60.0;

    auto zd = [&](double s0_, double alpha, double delta) -> double {
        double cz = std::sin(p.phi) * std::sin(delta) + std::cos(p.phi) * std::cos(delta) * std::cos(s0_ - alpha);
        cz = std::max(-1.0, std::min(1.0, cz));
        return std::acos(cz) * RAD;
    };
    double ze = zd(s0, p.ae, p.de);
    double zw = zd(s0, p.aw, p.dw);

    std::cout << "\n" << method << "\n";
    std::cout << "    s0  = " << std::fixed << std::setprecision(17) << s0 << " рад\n";
    std::cout << "    s0  = " << std::setprecision(10) << s0 * RAD << " °\n";
    std::cout << "    s0  = " << hh << "ч " << mm << "м "
              << std::setprecision(6) << ss << "с\n";
    std::cout << "    ze  = " << std::setprecision(10) << ze << " ° (восточная)\n";
    std::cout << "    zw  = " << zw << " ° (западная)\n";
    std::cout << "    |ze−zw| = " << std::scientific << std::setprecision(3)
              << std::abs(ze - zw) << " °\n";
    std::cout << "    f(s0)   = " << F(s0, p) << "\n";
}

Params readParams()
{
    Params p;
    p.phi= 20.0 * DEG;
    p.ae = (19.0 + 50.0/60.0 + 47.0/3600.0) * HOUR;
    p.de = (8.0 + 52.0/60.0) * DEG;
    p.aw= (18.0 + 36.0/60.0 + 56.0/3600.0) * HOUR;
    p.dw = (38.0 + 47.0/60.0) * DEG;
    
    std::cout << "Данные:\n"
              << "Восточная звезда Альтаир\n"
              << "Прямое восхождение α = 19h 50m 47s\n"
              << "Склонение δ = +8° 52'\n"
              << "Западная звезда Вега\n"
              << "Прямое восхождение α = 18h 36m 56s\n"
              << "Склонение δ = +38° 47'\n"
              << "Широта наблюдлателя φ = 20°\n\n";
    
    return p;
}

int main()
{
    try {
        Params p = readParams();
        Bracket br = findBracket(p);

        if (!br.found) {
            std::cerr <<
                "\n oшибка, нет момента, когда зенитные расстояния совпадают.\n";
            return 1;
        }

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Отрезок: [" << br.a << ", " << br.b << "] рад\n";
        std::cout << "f(a) = " << std::scientific << F(br.a, p)
                  << ", f(b) = " << F(br.b, p) << "\n";

        double x0 = 0.5 * (br.a + br.b);
        std::cout << "Начальное приближение x0 = "
                  << std::fixed << x0 << " рад\n";

        std::cout << "\n Метод бисекции:\n";
        double s0_bis = bisection(br.a, br.b, p);

        std::cout << "\n Метод Ньютона:\n";
        double s0_newt = newton(x0, p);

        std::cout << "\n Результат\n";

        printResult("Бисекция", s0_bis, p);
        printResult("Метод Ньютона", s0_newt, p);

        double diff = std::abs(s0_bis - s0_newt);
        std::cout << "\nРасхождение между методами: "
                  << std::scientific << std::setprecision(3) << diff
                  << " рад  (" << diff * RAD << " °)\n";

        std::cout <<
            "\nЗенитные расстояния двух звезд равны\n";

    } catch (const std::exception& e) {
        std::cerr << "\nИСКЛЮЧЕНИЕ: " << e.what() << "\n";
        return 1;
    }
    return 0;
}