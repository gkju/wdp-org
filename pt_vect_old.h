typedef struct Pt_Vector {
    int sz;
    int capacity;
    Point_t* tab;
} Pt_Vector_t;

Pt_Vector_t get_pt_vec(int sz) {
    int capacity = sz | 1;
    return (Pt_Vector_t) { .sz = sz, .capacity = capacity, .tab = malloc(sizeof(Point_t) * capacity) };
}

int get_sz(Pt_Vector_t vec) {
    return vec.sz;
}

void set_sz(Pt_Vector_t* vec, int sz) {
    vec->sz = sz;
}

void clear(Pt_Vector_t* vec) {
    vec->sz = 0;
}

void free_pt_vec(Pt_Vector_t vec) {
    return free(vec.tab);
}

void push_to_pt_vec(Pt_Vector_t* vec, Point_t p) {
    if(vec->capacity < vec->sz + 1) {
        vec->tab = realloc(vec->tab, 2 * vec->capacity * sizeof(Point_t));
        vec->capacity *= 2;
    }

    vec->tab[vec->sz++] = p;
}

Point_t get_from_pt_vec(Pt_Vector_t vec, int i) {
    return vec.tab[i];
}

Point_t set_pt_vec_at(Pt_Vector_t vec, int i, Point_t p) {
    assert(i < vec.sz);
    return vec.tab[i] = p;
}