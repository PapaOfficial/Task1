#include <sys/types.h>
#include <stdio.h>
#include <gmp.h>

//Структура "точка" содержит координаты точки
struct point{mpz_t x; mpz_t y; mpz_t z;};

//Функция сложения точек для формы Гессе (add2)
struct point addition(struct point p1, struct point p2, mpz_t p){
    struct point p3;
    mpz_init(p3.x);
    mpz_init(p3.y);
    mpz_init(p3.z);
    mpz_t mul;
    mpz_init(mul);
    
    mpz_mul(mul, p1.y, p2.x);
    mpz_mul(mul, mul, p1.y);
    mpz_mul(p3.x, mul, p2.z);
    mpz_mul(mul, p1.z, p2.y);
    mpz_mul(mul, mul, p1.x);
    mpz_mul(mul, mul, p2.y);
    mpz_sub(p3.x, p3.x, mul);
    
    mpz_mul(mul, p1.x, p2.z);
    mpz_mul(mul, mul, p1.x);
    mpz_mul(p3.y, mul, p2.y);
    mpz_mul(mul, p1.y, p2.x);
    mpz_mul(mul, mul, p1.z);
    mpz_mul(mul, mul, p2.x);
    mpz_sub(p3.y, p3.y, mul);
    
    mpz_mul(mul, p1.z, p2.y);
    mpz_mul(mul, mul, p1.z);
    mpz_mul(p3.z, mul, p2.x);
    mpz_mul(mul, p1.x, p2.z);
    mpz_mul(mul, mul, p1.y);
    mpz_mul(mul, mul, p2.z);
    mpz_sub(p3.z, p3.z, mul);
    
    mpz_fdiv_r(p3.x, p3.x, p);
    mpz_fdiv_r(p3.y, p3.y, p);
    mpz_fdiv_r(p3.z, p3.z, p);
    
    mpz_clear(mul);
        
    return p3;
}


//Функция удвоения точки для формы Гессе (dbl2)
struct point doubling(struct point p1, mpz_t p){
    struct point p2;
    mpz_init(p2.x);
    mpz_init(p2.y);
    mpz_init(p2.z);
    mpz_t xxx, yyy, zzz;
    mpz_init(xxx);
    mpz_init(yyy);
    mpz_init(zzz);
    
    mpz_mul(xxx, p1.x, p1.x);
    mpz_mul(xxx, xxx, p1.x);
    
    mpz_mul(yyy, p1.y, p1.y);
    mpz_mul(yyy, yyy, p1.y);
    
    mpz_mul(zzz, p1.z, p1.z);
    mpz_mul(zzz, zzz, p1.z);
    
    mpz_sub(p2.x, xxx, zzz);
    mpz_mul(p2.x, p2.x, p1.y);
    
    mpz_sub(p2.y, zzz, yyy);
    mpz_mul(p2.y, p2.y, p1.x);
    
    mpz_sub(p2.z, yyy, xxx);
    mpz_mul(p2.z, p2.z, p1.z);
    
    mpz_clear(xxx);
    mpz_clear(yyy);
    mpz_clear(zzz);
    
    mpz_fdiv_r(p2.x, p2.x, p);
    mpz_fdiv_r(p2.y, p2.y, p);
    mpz_fdiv_r(p2.z, p2.z, p);
    
    return p2;
}




//Функция, вычисляющая кратную точку по алгоритму "Лесенка Монтгомери",
//функция принимает число k, модуль p, точки P и O и возвращает кратную точку
struct point montgomery_ladder(struct point P, struct point O, mpz_t k, mpz_t p){
    mp_bitcnt_t popcount = mpz_popcount(k);
    long int i = 0;
    while(popcount > 0){
        if (mpz_tstbit(k, i) == 1){
            popcount = popcount - 1;
        }
        i = i + 1;
    }
    i = i - 1;
    //i - максимальный индекс единицы в двоичном представлении k
    
    struct point Q;
    mpz_init_set(Q.x, O.x);
    mpz_init_set(Q.y, O.y);
    mpz_init_set(Q.z, O.z);
    struct point R;
    mpz_init_set(R.x, P.x);
    mpz_init_set(R.y, P.y);
    mpz_init_set(R.z, P.z);
    
    while(i >= 0){
        if (mpz_tstbit(k, i) == 0){
            R = addition(R, Q, p);
            Q = doubling(Q, p);
        } else {
            Q = addition(Q, R, p);
            R = doubling(R, p);
        }
        i = i - 1;
    }
    return Q;
}



//Функция, вычисляющая кратную точку по бинарному алгоритму,
//функция принимает число k, модуль p, точки P и O и возвращает кратную точку
struct point binary_algorithm(struct point P, struct point O, mpz_t k, mpz_t p){
    mp_bitcnt_t popcount = mpz_popcount(k);
    long int i = 0;
    while(popcount > 0){
        if (mpz_tstbit(k, i) == 1){
            popcount = popcount - 1;
        }
        i = i + 1;
    }
    i = i - 1;
    //i - максимальный индекс единицы в двоичном представлении k
    
    struct point Q;
    mpz_init_set(Q.x, O.x);
    mpz_init_set(Q.y, O.y);
    mpz_init_set(Q.z, O.z);

    while(i >= 0){
        Q = doubling(Q, p);
        if (mpz_tstbit(k, i) == 1){
            Q = addition(Q, P, p);
        }
        i = i - 1;
    }
    return Q;
}









int main(void) {
    //Сначала определяем D, p, q, начальную точку P, (a, b для формы Вейерштрасса)
    //и точку бесконечности O
    mpz_t A, B, C, a, b, p, D, q, k;
    mpz_init_set_str(a,"-2835", 10);
    mpz_init_set_str(b,"9774", 10);
    mpz_init_set_str(p,"115792089237316195423570985008687907853269984665640564039457584007913111864739",10);
    mpz_init_set_str(D, "3", 10);
    mpz_init_set_str(q,"115792089237316195423570985008687907852907080286716537199505774922796921406320",10);
    mpz_init(A);
    mpz_init(B);
    mpz_init(C);
    mpz_init(k);
    
    struct point P;
    mpz_init_set_str(P.x, "93528601524582384978654134272795468222405403055717890280271688132874849008326", 10);
    mpz_init_set_str(P.y, "14443324612566128911211262381388707474030458136470034119105598903952521080679", 10);
    mpz_init_set_str(P.z, "1", 10);
    
    struct point O;
    mpz_init_set_ui(O.x, 1);
    mpz_init_set_si(O.y, -1);
    mpz_init_set_ui(O.z, 0);
    
    
    
    
    //генерируем псевдорандомное k
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_urandomm(k, state, p);
    
    
    //вычисляем кратную точку
    struct point Q;
    Q = montgomery_ladder(P, O, k, p);
    
    
    //тест 4.1
    //проверяем, что точка Q лежит на кривой
    mpz_mul(A, Q.x, Q.x);
    mpz_mul(A, A, Q.x);
    mpz_mul(B, Q.y, Q.y);
    mpz_mul(B, B, Q.y);
    mpz_mul(C, Q.z, Q.z);
    mpz_mul(C, C, Q.z);
    mpz_add(A, A, B);
    mpz_add(A, A, C);
    mpz_fdiv_r(A, A, p);
    printf("%s\n", "Тест 4.1:");
    printf("%s", "k=");
    gmp_printf("%Zd\n", k);
    printf("%s\n", "Левая часть уравнения:");
    gmp_printf("%Zd\n", A);

    mpz_mul_ui(B, D, 3);
    mpz_mul(B, B, Q.x);
    mpz_mul(B, B, Q.y);
    mpz_mul(B, B, Q.z);
    mpz_fdiv_r(B, B, p);
    printf("%s\n", "Правая часть уравнения:");
    gmp_printf("%Zd\n", B);
    
    
    //тест 4.2
    //q - порядок группы точек
    mpz_set(k, q);
    struct point Q1;
    Q1 = montgomery_ladder(P, O, k, p);
    
    printf("\n\n%s\n", "Тест 4.2:");
    printf("%s\n", "Координаты полученной точки:");
    printf("%s", "x=");
    gmp_printf("%Zd\n", Q1.x);
    printf("%s", "y=");
    gmp_printf("%Zd\n", Q1.y);
    printf("%s", "z=");
    gmp_printf("%Zd\n", Q1.z);
    printf("%s\n", "так как z=0, то это точка бесконечности");
    
    
    //тест 4.3
    //q - порядок группы точек
    mpz_set(k, q);
    mpz_add_ui(k, k, 1); //k = q + 1
    struct point Q2;
    Q2 = montgomery_ladder(P, O, k, p);
    
    //так как мы в проективном пространстве, то Q2 = P, если их координаты пропорциональны
    printf("\n\n%s\n", "Тест 4.3:");
    printf("%s\n", "Первое уравнение:");
    printf("%s\n", "Полученная точка равна точке P, если их координаты пропорциональны");
    
    mpz_mul(A, Q2.x, P.y);
    mpz_mul(B, Q2.y, P.x);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q2.x * P.y");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q2.y * P.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q2.x, P.z);
    mpz_mul(B, Q2.z, P.x);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q2.x * P.z");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q2.z * P.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q2.z, P.y);
    mpz_mul(B, Q2.y, P.z);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q2.z * P.y");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q2.y * P.z");
    gmp_printf("%Zd\n", B);
    
    
    printf("\n%s\n", "Второе уравнение:");
    printf("%s\n", "Если P=(x,y,z), то -P=(y,x,z), P1=-P");
    printf("%s\n", "Аналогично смотрим, чтобы координаты были пропорциональны: ");
    
    mpz_set(k, q);
    mpz_sub_ui(k, k, 1); //k = q - 1
    struct point Q3;
    Q3 = montgomery_ladder(P, O, k, p);
    
    //Q3 = -P, если их координаты пропорциональны
    mpz_mul(A, Q3.x, P.x);
    mpz_mul(B, Q3.y, P.y);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q3.x * P1.y");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q3.y * P1.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q3.x, P.z);
    mpz_mul(B, Q3.z, P.y);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q3.x * P1.z");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q3.z * P1.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q3.y, P.z);
    mpz_mul(B, Q3.z, P.x);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q3.y * P1.z");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q3.z * P1.y");
    gmp_printf("%Zd\n", B);
    
    
    //Тест 4.4
    //генерируем псевдорандомное k1
    mpz_t k1;
    mpz_init(k1);
    mpz_urandomm(k1, state, p);
    
    //вычисляем кратную точку от k1
    struct point Q4;
    Q4 = montgomery_ladder(P, O, k1, p);
    
    //генерируем псевдорандомное k2
    mpz_t k2;
    mpz_init(k2);
    mpz_urandomm(k2, state, p);
    
    //вычисляем кратную точку от k2
    struct point Q5;
    Q5 = montgomery_ladder(P, O, k2, p);
    
    //в Q4 пишем сумму точек Q4 и Q5
    Q4 = addition(Q4, Q5, p);
    
    //в k1 записываем сумму k1 и k2
    mpz_add(k1, k1, k2);
    mpz_fdiv_r(k1, k1, p);
    
    //в Q5 помещаем кратную точку от (k1 + k2)
    Q5 = montgomery_ladder(P, O, k1, p);
    
    //координаты Q4 и координаты Q5 должны быть пропорциональны
    printf("\n\n%s\n", "Тест 4.4:");
    printf("%s", "k1=");
    gmp_printf("%Zd\n", k1);
    printf("%s", "k2=");
    gmp_printf("%Zd\n", k2);
    printf("%s\n", "Q4 - точка слева знака равно, Q5 - точка справа");
    
    mpz_mul(A, Q4.x, Q5.y);
    mpz_mul(B, Q4.y, Q5.x);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q4.x * Q5.y");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q4.y * Q5.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q4.x, Q5.z);
    mpz_mul(B, Q4.z, Q5.x);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q4.x * Q5.z");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q4.z * Q5.x");
    gmp_printf("%Zd\n", B);
    
    mpz_mul(A, Q4.z, Q5.y);
    mpz_mul(B, Q4.y, Q5.z);
    mpz_fdiv_r(A, A, p);
    mpz_fdiv_r(B, B, p);
    
    printf("%s\n", "Q4.z * Q5.y");
    gmp_printf("%Zd\n", A);
    printf("%s\n", "Q4.y * Q5.z");
    gmp_printf("%Zd\n", B);
    
    //чистим за собой
    mpz_clears(A,B,C,a,b,p,D,q,k,P.x,P.y,P.z,O.x,O.y,O.z,Q.x,Q.y,Q.z,Q1.x,Q1.y,Q1.z,Q2.x,Q2.y,Q2.z,Q3.x,Q3.y,Q3.z,Q4.x,Q4.y,Q4.z,Q5.x,Q5.y,Q5.z, NULL);
    
    return 0;
}
