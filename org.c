#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <stdbool.h>
#include <assert.h>
/////////////////////////////////////////////////
// Code review robione przez Wojciech Sobiński //
/////////////////////////////////////////////////
bool is_zero(double x) {
    return fabs(x) < 1e-10;
}

bool is_equiv(double x, double y) {
    return is_zero(x - y);
}

bool is_gt(double x, double y) {
    return x > y && !is_equiv(x, y);
}

bool is_lt(double x, double y) {
    return x < y && !is_equiv(x, y);
}

// czy x >= y
bool is_geq(double x, double y) {
    return !is_lt(x, y);
}

// czy x <= y
bool is_leq(double x, double y) {
    return !is_gt(x, y);
}

double squared(double x) {
    return x * x;
}

typedef struct Point {
    double Re, Im;
} Point_t;

Point_t pt_from_cart(double x, double y) {
    return (Point_t) { .Re = x, .Im = y };
}

Point_t pt_from_real(double x) {
    return pt_from_cart(x, 0);
}

double Re(Point_t p) {
    return p.Re;
}

double Im(Point_t p) {
    return p.Im;
}

Point_t pt_conj(Point_t p) {
    return (Point_t) { .Re = Re(p), .Im = -Im(p) };
}

double norm_sq(Point_t p) {
    return squared(Re(p)) + squared(Im(p));
}

Point_t add(Point_t a, Point_t b) {
    return (Point_t) { .Re = Re(a) + Re(b), .Im = Im(a) + Im(b) };
}

Point_t mult(Point_t a, Point_t b) {
    return (Point_t) { .Re = Re(a) * Re(b) - Im(a) * Im(b), .Im = Re(a) * Im(b) + Im(a) * Re(b) };
}

Point_t inv(Point_t a) {
    Point_t x = pt_conj(a);
    double norm = norm_sq(a);
    return (Point_t) { .Re = Re(x) / norm, .Im = Im(x) / norm };
}

// a / b
Point_t pt_div(Point_t a, Point_t b) {
    return mult(a, inv(b));
}

// zwraca pkt odpowiadajacy wektorowi ab
Point_t vec_betw(Point_t a, Point_t b) {
    return (Point_t) { .Re = Re(b) - Re(a), .Im = Im(b) - Im(a) };
}

// traktuje a i b jako wektory na plaszczyznie euklidesowej i liczy ich iloczyn wektorowy
double pt_vec_mult(Point_t a, Point_t b) {
    return Im(a) * Re(b) - Re(a) * Im(b);
}

// prostokat, kolo, zgiecie kartki
typedef enum {
    Rect,
    Circ,
    Comp
} Paper_Type_t;

typedef struct Paper {
    Paper_Type_t type;
} Paper_t;

typedef struct Rect_Paper {
    Paper_t paper;
    // xl xr to wspolrzedne pionowych prostych, yb yt poziomych
    double xl, xr;
    double yb, yt;
} Rect_Paper_t;

typedef struct Circ_Paper {
    Paper_t paper;
    double radius;
    Point_t origin;
} Circ_Paper_t;

typedef struct Comp_Paper {
    Paper_t paper;
    // co zginano
    Paper_t* base;
    // a i b to pkty wzdluz ktorych zginano
    Point_t a, b;
} Comp_Paper_t;

bool in_paper(Point_t poi, const Paper_t* pap) {
    if(pap->type == Rect) {
        const Rect_Paper_t* rect = (Rect_Paper_t*) pap;
        return  is_leq(rect->xl, poi.Re) && 
                is_leq(poi.Re, rect->xr) && 
                is_leq(rect->yb, poi.Im) && 
                is_leq(poi.Im, rect->yt);
    }

    if(pap->type == Circ) {
        const Circ_Paper_t* circ = (Circ_Paper_t*) pap;
        return is_leq(norm_sq(vec_betw(circ->origin, poi)), squared(circ->radius));
    }

    assert(pap->type != Comp);

    return false;
}

// odbija t wzgledem odcinka ab
// uwaga a != b !
Point_t reflect(Point_t t, Point_t a, Point_t b) {
    t = vec_betw(a, t);
    b = vec_betw(a, b);

    // (liczba zespolona z okr jedn o argumencie b)^2
    Point_t x_sq = pt_div(mult(b, b), pt_from_real(norm_sq(b)));

    // liczy conj(t) * x^2 + a, co jest uproszczonym wzorem conj(t * conj(x)) * x + a, ktory trywialnie jest wzorem na odbicie t
    return add(mult(pt_conj(t), x_sq), a);
}

int cnt_points(Paper_t* pap, Point_t point) {
    if(pap->type != Comp) {
        return in_paper(point, pap);
    }

    Comp_Paper_t* comp = (Comp_Paper_t*) pap;
    // wektor odpowiadajacy prostej wzdluz ktorej zginalismy
    Point_t ab = vec_betw(comp->a, comp->b);
    // sprawdzamy po ktorej stronie zgiecia lezal
    double side = pt_vec_mult(ab, vec_betw(comp->a, point));

    // bez zmiennych pomocniczych aby byla tail rekurencja
    if(is_lt(side, 0)) {
        return cnt_points(comp->base, reflect(point, comp->a, comp->b)) + cnt_points(comp->base, point);
    } else if(is_leq(side, 0)) {
        return cnt_points(comp->base, point);
    }

    return 0;
}

void init_rect(Rect_Paper_t* rect, double xl, double yb, double xr, double yt) {
    rect->paper.type = Rect;
    rect->xl = xl;
    rect->xr = xr;
    rect->yb = yb;
    rect->yt = yt;
}

void init_circ(Circ_Paper_t* circ, double x, double y, double r) {
    circ->paper.type = Circ;
    circ->origin = pt_from_cart(x, y);
    circ->radius = r;
}

void init_comp(Comp_Paper_t* comp, Paper_t* base, double ax, double ay, double bx, double by) {
    comp->paper.type = Comp;
    comp->base = base;
    comp->a = pt_from_cart(ax, ay);
    comp->b = pt_from_cart(bx, by);
}

int main() {
    int n, q;
    scanf("%d %d", &n, &q);
    Paper_t** papers = malloc((size_t) n * sizeof(Paper_t*));

    for(int i = 0; i < n; ++i) {
        char type;
        // ignorujemy leading whitespace stad " %c"
        scanf(" %c", &type);

        switch(type) {
            case 'P': {
                double xl, yb, xr, yt;
                scanf("%lf %lf %lf %lf", &xl, &yb, &xr, &yt);
                papers[i] = malloc(sizeof(Rect_Paper_t));
                Rect_Paper_t* rect = (Rect_Paper_t*) papers[i];
                init_rect(rect, xl, yb, xr, yt);
                break;
            }
            case 'K': {
                double x, y, r;
                scanf("%lf %lf %lf", &x, &y, &r);
                papers[i] = malloc(sizeof(Circ_Paper_t));
                Circ_Paper_t* circ = (Circ_Paper_t*) papers[i];
                init_circ(circ, x, y, r);
                break;
            }
            case 'Z': {
                int k;
                double ax, ay, bx, by;
                scanf("%d %lf %lf %lf %lf", &k, &ax, &ay, &bx, &by);
                papers[i] = malloc(sizeof(Comp_Paper_t));
                Comp_Paper_t* comp = (Comp_Paper_t*) papers[i];
                init_comp(comp, papers[k - 1], ax, ay, bx, by);
                break;
            }
        };
    }

    for(int i = 0; i < q; ++i) {
        int k;
        double x, y;
        scanf("%d %lf %lf", &k, &x, &y);

        printf("%d\n", cnt_points(papers[k - 1], pt_from_cart(x, y)));
    }

    for(int i = 0; i < n; ++i) {
        free(papers[i]);
    }
    free(papers);
}
/////////////////////////////////////////////////
// Code review robione przez Wojciech Sobiński //
/////////////////////////////////////////////////